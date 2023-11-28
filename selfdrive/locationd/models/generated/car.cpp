#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_1881581342858382704) {
   out_1881581342858382704[0] = delta_x[0] + nom_x[0];
   out_1881581342858382704[1] = delta_x[1] + nom_x[1];
   out_1881581342858382704[2] = delta_x[2] + nom_x[2];
   out_1881581342858382704[3] = delta_x[3] + nom_x[3];
   out_1881581342858382704[4] = delta_x[4] + nom_x[4];
   out_1881581342858382704[5] = delta_x[5] + nom_x[5];
   out_1881581342858382704[6] = delta_x[6] + nom_x[6];
   out_1881581342858382704[7] = delta_x[7] + nom_x[7];
   out_1881581342858382704[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_1313133831038911976) {
   out_1313133831038911976[0] = -nom_x[0] + true_x[0];
   out_1313133831038911976[1] = -nom_x[1] + true_x[1];
   out_1313133831038911976[2] = -nom_x[2] + true_x[2];
   out_1313133831038911976[3] = -nom_x[3] + true_x[3];
   out_1313133831038911976[4] = -nom_x[4] + true_x[4];
   out_1313133831038911976[5] = -nom_x[5] + true_x[5];
   out_1313133831038911976[6] = -nom_x[6] + true_x[6];
   out_1313133831038911976[7] = -nom_x[7] + true_x[7];
   out_1313133831038911976[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_4584559026012520373) {
   out_4584559026012520373[0] = 1.0;
   out_4584559026012520373[1] = 0;
   out_4584559026012520373[2] = 0;
   out_4584559026012520373[3] = 0;
   out_4584559026012520373[4] = 0;
   out_4584559026012520373[5] = 0;
   out_4584559026012520373[6] = 0;
   out_4584559026012520373[7] = 0;
   out_4584559026012520373[8] = 0;
   out_4584559026012520373[9] = 0;
   out_4584559026012520373[10] = 1.0;
   out_4584559026012520373[11] = 0;
   out_4584559026012520373[12] = 0;
   out_4584559026012520373[13] = 0;
   out_4584559026012520373[14] = 0;
   out_4584559026012520373[15] = 0;
   out_4584559026012520373[16] = 0;
   out_4584559026012520373[17] = 0;
   out_4584559026012520373[18] = 0;
   out_4584559026012520373[19] = 0;
   out_4584559026012520373[20] = 1.0;
   out_4584559026012520373[21] = 0;
   out_4584559026012520373[22] = 0;
   out_4584559026012520373[23] = 0;
   out_4584559026012520373[24] = 0;
   out_4584559026012520373[25] = 0;
   out_4584559026012520373[26] = 0;
   out_4584559026012520373[27] = 0;
   out_4584559026012520373[28] = 0;
   out_4584559026012520373[29] = 0;
   out_4584559026012520373[30] = 1.0;
   out_4584559026012520373[31] = 0;
   out_4584559026012520373[32] = 0;
   out_4584559026012520373[33] = 0;
   out_4584559026012520373[34] = 0;
   out_4584559026012520373[35] = 0;
   out_4584559026012520373[36] = 0;
   out_4584559026012520373[37] = 0;
   out_4584559026012520373[38] = 0;
   out_4584559026012520373[39] = 0;
   out_4584559026012520373[40] = 1.0;
   out_4584559026012520373[41] = 0;
   out_4584559026012520373[42] = 0;
   out_4584559026012520373[43] = 0;
   out_4584559026012520373[44] = 0;
   out_4584559026012520373[45] = 0;
   out_4584559026012520373[46] = 0;
   out_4584559026012520373[47] = 0;
   out_4584559026012520373[48] = 0;
   out_4584559026012520373[49] = 0;
   out_4584559026012520373[50] = 1.0;
   out_4584559026012520373[51] = 0;
   out_4584559026012520373[52] = 0;
   out_4584559026012520373[53] = 0;
   out_4584559026012520373[54] = 0;
   out_4584559026012520373[55] = 0;
   out_4584559026012520373[56] = 0;
   out_4584559026012520373[57] = 0;
   out_4584559026012520373[58] = 0;
   out_4584559026012520373[59] = 0;
   out_4584559026012520373[60] = 1.0;
   out_4584559026012520373[61] = 0;
   out_4584559026012520373[62] = 0;
   out_4584559026012520373[63] = 0;
   out_4584559026012520373[64] = 0;
   out_4584559026012520373[65] = 0;
   out_4584559026012520373[66] = 0;
   out_4584559026012520373[67] = 0;
   out_4584559026012520373[68] = 0;
   out_4584559026012520373[69] = 0;
   out_4584559026012520373[70] = 1.0;
   out_4584559026012520373[71] = 0;
   out_4584559026012520373[72] = 0;
   out_4584559026012520373[73] = 0;
   out_4584559026012520373[74] = 0;
   out_4584559026012520373[75] = 0;
   out_4584559026012520373[76] = 0;
   out_4584559026012520373[77] = 0;
   out_4584559026012520373[78] = 0;
   out_4584559026012520373[79] = 0;
   out_4584559026012520373[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_5101551619635677519) {
   out_5101551619635677519[0] = state[0];
   out_5101551619635677519[1] = state[1];
   out_5101551619635677519[2] = state[2];
   out_5101551619635677519[3] = state[3];
   out_5101551619635677519[4] = state[4];
   out_5101551619635677519[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_5101551619635677519[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_5101551619635677519[7] = state[7];
   out_5101551619635677519[8] = state[8];
}
void F_fun(double *state, double dt, double *out_4910919501640013458) {
   out_4910919501640013458[0] = 1;
   out_4910919501640013458[1] = 0;
   out_4910919501640013458[2] = 0;
   out_4910919501640013458[3] = 0;
   out_4910919501640013458[4] = 0;
   out_4910919501640013458[5] = 0;
   out_4910919501640013458[6] = 0;
   out_4910919501640013458[7] = 0;
   out_4910919501640013458[8] = 0;
   out_4910919501640013458[9] = 0;
   out_4910919501640013458[10] = 1;
   out_4910919501640013458[11] = 0;
   out_4910919501640013458[12] = 0;
   out_4910919501640013458[13] = 0;
   out_4910919501640013458[14] = 0;
   out_4910919501640013458[15] = 0;
   out_4910919501640013458[16] = 0;
   out_4910919501640013458[17] = 0;
   out_4910919501640013458[18] = 0;
   out_4910919501640013458[19] = 0;
   out_4910919501640013458[20] = 1;
   out_4910919501640013458[21] = 0;
   out_4910919501640013458[22] = 0;
   out_4910919501640013458[23] = 0;
   out_4910919501640013458[24] = 0;
   out_4910919501640013458[25] = 0;
   out_4910919501640013458[26] = 0;
   out_4910919501640013458[27] = 0;
   out_4910919501640013458[28] = 0;
   out_4910919501640013458[29] = 0;
   out_4910919501640013458[30] = 1;
   out_4910919501640013458[31] = 0;
   out_4910919501640013458[32] = 0;
   out_4910919501640013458[33] = 0;
   out_4910919501640013458[34] = 0;
   out_4910919501640013458[35] = 0;
   out_4910919501640013458[36] = 0;
   out_4910919501640013458[37] = 0;
   out_4910919501640013458[38] = 0;
   out_4910919501640013458[39] = 0;
   out_4910919501640013458[40] = 1;
   out_4910919501640013458[41] = 0;
   out_4910919501640013458[42] = 0;
   out_4910919501640013458[43] = 0;
   out_4910919501640013458[44] = 0;
   out_4910919501640013458[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_4910919501640013458[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_4910919501640013458[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4910919501640013458[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4910919501640013458[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_4910919501640013458[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_4910919501640013458[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_4910919501640013458[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_4910919501640013458[53] = -9.8000000000000007*dt;
   out_4910919501640013458[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_4910919501640013458[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_4910919501640013458[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4910919501640013458[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4910919501640013458[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_4910919501640013458[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_4910919501640013458[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_4910919501640013458[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4910919501640013458[62] = 0;
   out_4910919501640013458[63] = 0;
   out_4910919501640013458[64] = 0;
   out_4910919501640013458[65] = 0;
   out_4910919501640013458[66] = 0;
   out_4910919501640013458[67] = 0;
   out_4910919501640013458[68] = 0;
   out_4910919501640013458[69] = 0;
   out_4910919501640013458[70] = 1;
   out_4910919501640013458[71] = 0;
   out_4910919501640013458[72] = 0;
   out_4910919501640013458[73] = 0;
   out_4910919501640013458[74] = 0;
   out_4910919501640013458[75] = 0;
   out_4910919501640013458[76] = 0;
   out_4910919501640013458[77] = 0;
   out_4910919501640013458[78] = 0;
   out_4910919501640013458[79] = 0;
   out_4910919501640013458[80] = 1;
}
void h_25(double *state, double *unused, double *out_2272755813087694232) {
   out_2272755813087694232[0] = state[6];
}
void H_25(double *state, double *unused, double *out_5686421704720180986) {
   out_5686421704720180986[0] = 0;
   out_5686421704720180986[1] = 0;
   out_5686421704720180986[2] = 0;
   out_5686421704720180986[3] = 0;
   out_5686421704720180986[4] = 0;
   out_5686421704720180986[5] = 0;
   out_5686421704720180986[6] = 1;
   out_5686421704720180986[7] = 0;
   out_5686421704720180986[8] = 0;
}
void h_24(double *state, double *unused, double *out_8684099650987606346) {
   out_8684099650987606346[0] = state[4];
   out_8684099650987606346[1] = state[5];
}
void H_24(double *state, double *unused, double *out_5341930735312443737) {
   out_5341930735312443737[0] = 0;
   out_5341930735312443737[1] = 0;
   out_5341930735312443737[2] = 0;
   out_5341930735312443737[3] = 0;
   out_5341930735312443737[4] = 1;
   out_5341930735312443737[5] = 0;
   out_5341930735312443737[6] = 0;
   out_5341930735312443737[7] = 0;
   out_5341930735312443737[8] = 0;
   out_5341930735312443737[9] = 0;
   out_5341930735312443737[10] = 0;
   out_5341930735312443737[11] = 0;
   out_5341930735312443737[12] = 0;
   out_5341930735312443737[13] = 0;
   out_5341930735312443737[14] = 1;
   out_5341930735312443737[15] = 0;
   out_5341930735312443737[16] = 0;
   out_5341930735312443737[17] = 0;
}
void h_30(double *state, double *unused, double *out_1667400207738112234) {
   out_1667400207738112234[0] = state[4];
}
void H_30(double *state, double *unused, double *out_1158725374592572788) {
   out_1158725374592572788[0] = 0;
   out_1158725374592572788[1] = 0;
   out_1158725374592572788[2] = 0;
   out_1158725374592572788[3] = 0;
   out_1158725374592572788[4] = 1;
   out_1158725374592572788[5] = 0;
   out_1158725374592572788[6] = 0;
   out_1158725374592572788[7] = 0;
   out_1158725374592572788[8] = 0;
}
void h_26(double *state, double *unused, double *out_7353143204702364208) {
   out_7353143204702364208[0] = state[7];
}
void H_26(double *state, double *unused, double *out_1944918385846124762) {
   out_1944918385846124762[0] = 0;
   out_1944918385846124762[1] = 0;
   out_1944918385846124762[2] = 0;
   out_1944918385846124762[3] = 0;
   out_1944918385846124762[4] = 0;
   out_1944918385846124762[5] = 0;
   out_1944918385846124762[6] = 0;
   out_1944918385846124762[7] = 1;
   out_1944918385846124762[8] = 0;
}
void h_27(double *state, double *unused, double *out_2872134250027953305) {
   out_2872134250027953305[0] = state[3];
}
void H_27(double *state, double *unused, double *out_1016037937207852123) {
   out_1016037937207852123[0] = 0;
   out_1016037937207852123[1] = 0;
   out_1016037937207852123[2] = 0;
   out_1016037937207852123[3] = 1;
   out_1016037937207852123[4] = 0;
   out_1016037937207852123[5] = 0;
   out_1016037937207852123[6] = 0;
   out_1016037937207852123[7] = 0;
   out_1016037937207852123[8] = 0;
}
void h_29(double *state, double *unused, double *out_3515493882213917574) {
   out_3515493882213917574[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1668956718906964972) {
   out_1668956718906964972[0] = 0;
   out_1668956718906964972[1] = 1;
   out_1668956718906964972[2] = 0;
   out_1668956718906964972[3] = 0;
   out_1668956718906964972[4] = 0;
   out_1668956718906964972[5] = 0;
   out_1668956718906964972[6] = 0;
   out_1668956718906964972[7] = 0;
   out_1668956718906964972[8] = 0;
}
void h_28(double *state, double *unused, double *out_1785407601671488978) {
   out_1785407601671488978[0] = state[0];
}
void H_28(double *state, double *unused, double *out_3632586990472291223) {
   out_3632586990472291223[0] = 1;
   out_3632586990472291223[1] = 0;
   out_3632586990472291223[2] = 0;
   out_3632586990472291223[3] = 0;
   out_3632586990472291223[4] = 0;
   out_3632586990472291223[5] = 0;
   out_3632586990472291223[6] = 0;
   out_3632586990472291223[7] = 0;
   out_3632586990472291223[8] = 0;
}
void h_31(double *state, double *unused, double *out_5490509443732978209) {
   out_5490509443732978209[0] = state[8];
}
void H_31(double *state, double *unused, double *out_1318710283612773286) {
   out_1318710283612773286[0] = 0;
   out_1318710283612773286[1] = 0;
   out_1318710283612773286[2] = 0;
   out_1318710283612773286[3] = 0;
   out_1318710283612773286[4] = 0;
   out_1318710283612773286[5] = 0;
   out_1318710283612773286[6] = 0;
   out_1318710283612773286[7] = 0;
   out_1318710283612773286[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_1881581342858382704) {
  err_fun(nom_x, delta_x, out_1881581342858382704);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_1313133831038911976) {
  inv_err_fun(nom_x, true_x, out_1313133831038911976);
}
void car_H_mod_fun(double *state, double *out_4584559026012520373) {
  H_mod_fun(state, out_4584559026012520373);
}
void car_f_fun(double *state, double dt, double *out_5101551619635677519) {
  f_fun(state,  dt, out_5101551619635677519);
}
void car_F_fun(double *state, double dt, double *out_4910919501640013458) {
  F_fun(state,  dt, out_4910919501640013458);
}
void car_h_25(double *state, double *unused, double *out_2272755813087694232) {
  h_25(state, unused, out_2272755813087694232);
}
void car_H_25(double *state, double *unused, double *out_5686421704720180986) {
  H_25(state, unused, out_5686421704720180986);
}
void car_h_24(double *state, double *unused, double *out_8684099650987606346) {
  h_24(state, unused, out_8684099650987606346);
}
void car_H_24(double *state, double *unused, double *out_5341930735312443737) {
  H_24(state, unused, out_5341930735312443737);
}
void car_h_30(double *state, double *unused, double *out_1667400207738112234) {
  h_30(state, unused, out_1667400207738112234);
}
void car_H_30(double *state, double *unused, double *out_1158725374592572788) {
  H_30(state, unused, out_1158725374592572788);
}
void car_h_26(double *state, double *unused, double *out_7353143204702364208) {
  h_26(state, unused, out_7353143204702364208);
}
void car_H_26(double *state, double *unused, double *out_1944918385846124762) {
  H_26(state, unused, out_1944918385846124762);
}
void car_h_27(double *state, double *unused, double *out_2872134250027953305) {
  h_27(state, unused, out_2872134250027953305);
}
void car_H_27(double *state, double *unused, double *out_1016037937207852123) {
  H_27(state, unused, out_1016037937207852123);
}
void car_h_29(double *state, double *unused, double *out_3515493882213917574) {
  h_29(state, unused, out_3515493882213917574);
}
void car_H_29(double *state, double *unused, double *out_1668956718906964972) {
  H_29(state, unused, out_1668956718906964972);
}
void car_h_28(double *state, double *unused, double *out_1785407601671488978) {
  h_28(state, unused, out_1785407601671488978);
}
void car_H_28(double *state, double *unused, double *out_3632586990472291223) {
  H_28(state, unused, out_3632586990472291223);
}
void car_h_31(double *state, double *unused, double *out_5490509443732978209) {
  h_31(state, unused, out_5490509443732978209);
}
void car_H_31(double *state, double *unused, double *out_1318710283612773286) {
  H_31(state, unused, out_1318710283612773286);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_lib_init(car)
