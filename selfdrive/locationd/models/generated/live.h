#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_35(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_3802695452116108024);
void live_err_fun(double *nom_x, double *delta_x, double *out_2043090151300958853);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_4795362139913636671);
void live_H_mod_fun(double *state, double *out_3957003700060578281);
void live_f_fun(double *state, double dt, double *out_527245089978267491);
void live_F_fun(double *state, double dt, double *out_4485576451842111367);
void live_h_4(double *state, double *unused, double *out_5083618636369239354);
void live_H_4(double *state, double *unused, double *out_5537439906076635471);
void live_h_9(double *state, double *unused, double *out_5191441572416926564);
void live_H_9(double *state, double *unused, double *out_5778629552706226116);
void live_h_10(double *state, double *unused, double *out_4694201907022710331);
void live_H_10(double *state, double *unused, double *out_3014558940665275568);
void live_h_12(double *state, double *unused, double *out_7779712826253544198);
void live_H_12(double *state, double *unused, double *out_7889847759600954350);
void live_h_35(double *state, double *unused, double *out_1653953902525287304);
void live_H_35(double *state, double *unused, double *out_8904101963449242847);
void live_h_32(double *state, double *unused, double *out_6298947203429338932);
void live_H_32(double *state, double *unused, double *out_4563097187716135959);
void live_h_13(double *state, double *unused, double *out_6115061610492001857);
void live_H_13(double *state, double *unused, double *out_5775805087570586313);
void live_h_14(double *state, double *unused, double *out_5191441572416926564);
void live_H_14(double *state, double *unused, double *out_5778629552706226116);
void live_h_33(double *state, double *unused, double *out_4623361709124221677);
void live_H_33(double *state, double *unused, double *out_6392085105621451165);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}