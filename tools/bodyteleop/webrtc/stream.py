import asyncio
import aiortc
from aiortc.contrib.media import MediaRelay

import abc
import dataclasses
import logging
from typing import Callable, Awaitable, List, Any, Optional


@dataclasses.dataclass
class StreamingOffer:
  sdp: str
  video: List[str]


ConnectionProvider = Callable[[StreamingOffer], Awaitable[aiortc.RTCSessionDescription]]


class WebRTCBaseStream(abc.ABC):
  def __init__(self,
               consumed_camera_types: List[str],
               consume_audio: bool,
               video_producer_tracks: List[aiortc.MediaStreamTrack],
               audio_producer_tracks: List[aiortc.MediaStreamTrack],
               add_data_channel: bool):
    self.peer_connection = aiortc.RTCPeerConnection()
    self.media_relay = MediaRelay()
    self.expected_incoming_camera_types = consumed_camera_types
    self.expected_incoming_audio = consume_audio
    self.expected_number_of_incoming_media = None

    self.incoming_camera_tracks = dict()
    self.incoming_audio_tracks = []
    self.outgoing_video_tracks = video_producer_tracks
    self.outgoing_audio_tracks = audio_producer_tracks

    self.add_data_channel = add_data_channel
    self.messaging_channel = None
    self.incoming_message_handlers = []

    self.incoming_media_ready_event = asyncio.Event()
    self.connection_attempted_event = asyncio.Event()

    self.peer_connection.on("connectionstatechange", self._on_connectionstatechange)
    self.peer_connection.on("datachannel", self._on_incoming_datachannel)
    self.peer_connection.on("track", self._on_incoming_track)

    self.logger = logging.getLogger("WebRTCStream")
    self.logger.setLevel(logging.DEBUG)

  def _log_debug(self, msg: Any, *args):
    self.logger.debug(f"{type(self)}() {msg}", *args)

  @property
  def _number_of_incoming_media(self):
    media = len(self.incoming_camera_tracks) + len(self.incoming_audio_tracks)
    media += int(self.messaging_channel is not None) if not self.add_data_channel else 0
    return media

  def _add_consumer_transceivers(self):
    for _ in self.expected_incoming_camera_types:
      self.peer_connection.addTransceiver("video", direction="recvonly")
    if self.expected_incoming_audio:
      self.peer_connection.addTransceiver("audio", direction="recvonly")

  def _add_producer_tracks(self):
    for track in self.outgoing_video_tracks:
      sender = self.peer_connection.addTrack(track)
      if hasattr(track, "codec_preference") and track.codec_preference() is not None:
        transceiver = next(t for t in self.peer_connection.getTransceivers() if t.sender == sender)
        self._force_codec(transceiver, track.codec_preference(), "video")
    for track in self.outgoing_audio_tracks:
      self.peer_connection.addTrack(track)

  def _add_messaging_channel(self, channel: Optional[aiortc.RTCDataChannel] = None):
    if not channel:
      channel = self.peer_connection.createDataChannel("data", ordered=True)

    for handler in self.incoming_message_handlers:
      channel.on("message", self._create_channel_handler_wrapper(channel, handler))
    self.messaging_channel = channel

  def _create_channel_handler_wrapper(self, channel, message_handler):
    async def handler_wrapper(message):
      await message_handler(channel, message)
    return handler_wrapper

  def _force_codec(self, transceiver, codec, stream_type):
    codec_mime = f"{stream_type}/{codec.upper()}"
    codecs = aiortc.RTCRtpSender.getCapabilities(stream_type).codecs
    codec = [codec for codec in codecs if codec.mimeType == codec_mime]
    transceiver.setCodecPreferences(codec)

  def _on_connectionstatechange(self):
    self._log_debug("connection state is", self.peer_connection.connectionState)
    if self.peer_connection.connectionState in ['connected', 'failed']:
      self.connection_attempted_event.set()

  def _on_incoming_track(self, track: aiortc.MediaStreamTrack):
    self._log_debug("got track:", track.kind, track.id)
    if track.kind == "video":
      parts = track.id.split(":") # format: "camera_type:camera_id"
      if len(parts) < 2:
        return

      camera_type = parts[0]
      if camera_type in self.expected_incoming_camera_types:
        self.incoming_camera_tracks[camera_type] = track
    elif track.kind == "audio":
      if self.expected_incoming_audio:
        self.incoming_audio_tracks.append(track)
    self._on_after_media()

  def _on_incoming_datachannel(self, channel: aiortc.RTCDataChannel):
    self._log_debug("got data channel:", channel.label)
    if channel.label == "data" and self.messaging_channel is None:
      self._add_messaging_channel(channel)
    self._on_after_media()

  def _on_after_media(self):
    if self._number_of_incoming_media == self.expected_number_of_incoming_media:
      self.incoming_media_ready_event.set()

  def _parse_incoming_streams(self, remote_sdp: str):
    desc = aiortc.sdp.SessionDescription.parse(remote_sdp)
    sending_medias = [m for m in desc.media if m.direction in ["sendonly", "sendrecv"]]
    incoming_media_count = len(sending_medias)
    if not self.add_data_channel:
      channel_medias = [m for m in desc.media if m.kind == "application"]
      incoming_media_count += len(channel_medias)
    self.expected_number_of_incoming_media = incoming_media_count

  def has_incoming_video_track(self, camera_type: str) -> bool:
    return camera_type in self.incoming_camera_tracks

  def has_incoming_audio_track(self) -> bool:
    return len(self.incoming_audio_tracks) > 0

  def has_messaging_channel(self) -> bool:
    return self.messaging_channel is not None

  def get_incoming_video_track(self, camera_type: str, buffered: bool):
    assert camera_type in self.incoming_camera_tracks, "Video tracks are not enabled on this stream"
    assert self.is_started, "Stream must be started"

    track = self.incoming_camera_tracks[camera_type]
    relay_track = self.media_relay.subscribe(track, buffered=buffered)
    return relay_track

  def get_incoming_audio_track(self, buffered: bool):
    assert len(self.incoming_audio_tracks) > 0, "Audio tracks are not enabled on this stream"
    assert self.is_started, "Stream must be started"

    track = self.incoming_audio_tracks[0]
    relay_track = self.media_relay.subscribe(track, buffered=buffered)
    return relay_track

  def get_messaging_channel(self):
    assert self.messaging_channel is not None, "Messaging channel is not enabled on this stream"
    assert self.is_started, "Stream must be started"

    return self.messaging_channel

  def set_message_handler(self, message_handler: Callable[[aiortc.RTCDataChannel, bytes], Awaitable[None]]):
    self.incoming_message_handlers.append(message_handler)
    if self.is_started:
      self.messaging_channel.on("message", self._create_channel_handler_wrapper(self.messaging_channel, message_handler))

  @property
  def is_started(self) -> bool:
    return self.peer_connection is not None and \
           self.peer_connection.localDescription is not None and \
           self.peer_connection.remoteDescription is not None

  async def wait_for_connection(self):
    await self.connection_attempted_event.wait()
    if self.peer_connection.connectionState != 'connected':
      raise ValueError("Connection failed.")
    if self.expected_number_of_incoming_media:
      await self.incoming_media_ready_event.wait()

  async def stop(self):
    await self.peer_connection.close()

  @abc.abstractmethod
  async def start(self):
    raise NotImplemented


class WebRTCOfferStream(WebRTCBaseStream):
  def __init__(self, session_provider: ConnectionProvider, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.session_provider = session_provider

  async def start(self):
    self._add_consumer_transceivers()
    if self.add_data_channel:
      self._add_messaging_channel()
    self._add_producer_tracks()

    offer = await self.peer_connection.createOffer()
    await self.peer_connection.setLocalDescription(offer)
    actual_offer = self.peer_connection.localDescription

    streaming_offer = StreamingOffer(
      sdp=actual_offer.sdp,
      video=list(self.expected_incoming_camera_types),
    )
    remote_answer = await self.session_provider(streaming_offer)
    self._parse_incoming_streams(remote_sdp=remote_answer.sdp)
    await self.peer_connection.setRemoteDescription(remote_answer)
    actual_answer = self.peer_connection.remoteDescription

    return actual_answer


class WebRTCAnswerStream(WebRTCBaseStream):
  def __init__(self, session: aiortc.RTCSessionDescription, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.session = session

  def _probe_video_codecs(self) -> List[str]:
    codecs = []
    for track in self.outgoing_video_tracks:
      if hasattr(track, "codec_preference") and track.codec_preference() is not None:
        codecs.append(track.codec_preference())

    return codecs

  def _override_incoming_video_codecs(self, remote_sdp: str, preferred_codecs: List[str]) -> str:
    desc = aiortc.sdp.SessionDescription.parse(remote_sdp)
    preferred_codec_mimes = [f"video/{c}" for c in preferred_codecs]
    for m in desc.media:
      if m.kind != "video":
        continue

      preferred_codecs = [c for c in m.rtp.codecs if c.mimeType in preferred_codec_mimes]
      if len(preferred_codecs) == 0:
        raise ValueError(f"None of {preferred_codecs} codecs is supported in remote SDP")

      m.rtp.codecs = preferred_codecs
      m.fmt = [c.payloadType for c in preferred_codecs]

    return str(desc)

  async def start(self):
    assert self.peer_connection.remoteDescription is None, "Connection already established"

    self._add_consumer_transceivers()

    # since we sent already encoded frames in some cases (e.g. livestream video tracks are in H264), we need to force aiortc to actually use it
    # we do that by overriding supported codec information on incoming sdp
    preferred_codecs = self._probe_video_codecs()
    if len(preferred_codecs) > 0:
      self.session.sdp = self._override_incoming_video_codecs(self.session.sdp, preferred_codecs)

    self._parse_incoming_streams(remote_sdp=self.session.sdp)
    await self.peer_connection.setRemoteDescription(self.session)

    self._add_producer_tracks()

    answer = await self.peer_connection.createAnswer()
    await self.peer_connection.setLocalDescription(answer)
    actual_answer = self.peer_connection.localDescription

    return actual_answer
