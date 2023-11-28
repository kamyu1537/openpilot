#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_1881581342858382704);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_1313133831038911976);
void car_H_mod_fun(double *state, double *out_4584559026012520373);
void car_f_fun(double *state, double dt, double *out_5101551619635677519);
void car_F_fun(double *state, double dt, double *out_4910919501640013458);
void car_h_25(double *state, double *unused, double *out_2272755813087694232);
void car_H_25(double *state, double *unused, double *out_5686421704720180986);
void car_h_24(double *state, double *unused, double *out_8684099650987606346);
void car_H_24(double *state, double *unused, double *out_5341930735312443737);
void car_h_30(double *state, double *unused, double *out_1667400207738112234);
void car_H_30(double *state, double *unused, double *out_1158725374592572788);
void car_h_26(double *state, double *unused, double *out_7353143204702364208);
void car_H_26(double *state, double *unused, double *out_1944918385846124762);
void car_h_27(double *state, double *unused, double *out_2872134250027953305);
void car_H_27(double *state, double *unused, double *out_1016037937207852123);
void car_h_29(double *state, double *unused, double *out_3515493882213917574);
void car_H_29(double *state, double *unused, double *out_1668956718906964972);
void car_h_28(double *state, double *unused, double *out_1785407601671488978);
void car_H_28(double *state, double *unused, double *out_3632586990472291223);
void car_h_31(double *state, double *unused, double *out_5490509443732978209);
void car_H_31(double *state, double *unused, double *out_1318710283612773286);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}