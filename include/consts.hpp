/**
 *   \file consts.hpp
 *   \brief useful constants
*/

#pragma once


#include <array>


namespace s13 {

  // math
  constexpr double gk_d2r = 0.0174533;

  // he6 run numbers
  constexpr int gk_he6_start_run = 133;
  constexpr int gk_he6_finish_run = 269;
  constexpr int gk_he6_num_of_runs = gk_he6_finish_run - gk_he6_start_run + 1;
  constexpr int gk_he6_pol_up_start_run = 133;
  constexpr int gk_he6_pol_up_finish_run = 193;
  constexpr int gk_he6_pol_down_start_run = 194;
  constexpr int gk_he6_pol_down_finish_run = 269;

  // he4 run numbers
  constexpr int gk_he4_start_run = 272;
  constexpr int gk_he4_finish_run = 316;
  constexpr int gk_he4_num_of_runs = gk_he4_finish_run - gk_he4_start_run + 1;
  constexpr int gk_he4_pol_down_start_run = 272;
  constexpr int gk_he4_pol_down_finish_run = 296;
  constexpr int gk_he4_pol_up_start_run = 297;
  constexpr int gk_he4_pol_up_finish_run = 316;

  // carbon run numbers
  constexpr int gk_carbon_start_run = 324;
  constexpr int gk_carbon_finish_run = 350;
  constexpr int gk_carbon_num_of_runs = gk_carbon_finish_run - gk_carbon_start_run + 1;

  // geometry of the setup

  constexpr double gk_dist_bdcs = 750.;
  constexpr double gk_dist_bdc_target = 1130;
  // thickness of the target [mm]
  constexpr double gk_target_z = 2.6;

  // S1 dimensions
  constexpr double gk_s1dc_width = 250 * 2;
  constexpr double gk_s1dc_height = 130 * 2;

  // FDC0 dimensions
  constexpr double gk_fdc0_width = 80 * 2;
  constexpr double gk_fdc0_height = 80 * 2;

  // BDC dimensions
  constexpr double gk_bdc_width = 80 * 2;
  constexpr double gk_bdc_height = 80 * 2;

  constexpr double gk_dist_target_fdc0 = 399.5;
  constexpr double gk_dist_fdc0_s1dc1 = 272.;
  constexpr double gk_dist_fdc0_s1dc2 = 400.;
  constexpr double gk_dist_s1dc_fdc0 = 336;
  constexpr double gk_dist_fdc0_target = 399.5;

  // ESPRIs
  constexpr int gk_espri_num_nais = 7;
  constexpr int gk_num_espris = 2;
  constexpr int gk_total_num_nais = gk_espri_num_nais * gk_num_espris;
  constexpr int gk_nai_pmts_num = gk_total_num_nais * 2;
  constexpr double gk_rdc2nai_esl_solid_angle_factor = 0.662094;
  constexpr double gk_rdc2nai_esr_solid_angle_factor = 0.627221;
  constexpr double gk_rdc_center_angle = 62.5;
  constexpr double gk_dist_rdc_target = 995 + 10; // 995.; // 1004.75?
  constexpr double gk_espri_acceptance_width = 145 * 2;
  constexpr double gk_espri_rdc_width = 220 * 2;
  constexpr double gk_espri_plastic_thickness = 4;
  constexpr double gk_espri_left_aang_a = -20.368;
  constexpr double gk_espri_left_aang_b = 1277.491;
  constexpr double gk_espri_right_aang_a = -14.646;
  constexpr double gk_espri_right_aang_b = 897.863;
  constexpr double gk_espri_computed_aang_a = -17.532;
  constexpr double gk_espri_computed_aang_b = 1095.687;
  constexpr std::array<double, gk_espri_num_nais>
  gk_espri_left_e_rel_mean_lin_fit_consts{{0, 0, 0, 0, 0, 0, 0}};
  constexpr std::array<double, gk_espri_num_nais>
  gk_espri_right_e_rel_mean_lin_fit_consts{{0, 0, 0, 0, 0, 0, 0}};

  // nais geometry
  constexpr double gk_dist_rdc_naiarray_center = 340;
  constexpr double gk_espri_naiarray_planes_offset = 65; /*!< offset of nai planes from the center */
  constexpr double gk_dist_naiarray_target = gk_dist_rdc_target + gk_dist_rdc_naiarray_center;

  // FDC2 geometry
  // TODO: clarify constant values
  constexpr double gk_fdc2_acceptance_height = 780;
  constexpr double gk_fdc2_acceptance_width = 2400;

  // HODF constants
  // threshold Q/dE value when HODF plastic is considered to be fired
  constexpr double gk_hodf_fire_threshold_q = 0.15;

  // cross-section computation

  // default setting for Y-acceptance of ESPRIs for cross-section computation
  constexpr double gk_cs_espri_y_acceptance = 132;
  // beam downscaling factor
  constexpr double gk_beam_dsc_factor = 2000;

  // total number of incident particles for he4 runs
  // R < 6 mm
  constexpr double gk_he4_total_incident = 2.855e6 * gk_beam_dsc_factor;
  // R < 10 mm
  //constexpr double gk_he4_total_incident = 5.248e6 * gk_beam_dsc_factor;

  // total number of incident particles for he6 runs
  // R < 6 mm
  constexpr double gk_he6_total_incident = 1.497e7 * gk_beam_dsc_factor;
  // R < 10 mm
  // constexpr double gk_he6_total_incident = 2.345e7 * gk_beam_dsc_factor;

  // areal target density number
  constexpr double gk_pol_tgt_areal_number_density = 1.046e22;
  constexpr double gk_carbon_tgt_areal_number_density = 9.98e21;

  // carbon background normalization factor
  constexpr double gk_carbon_bg_norm_factor = 10.706;
  constexpr double gk_carbon_to_he4_bg_norm_factor = 2.6;

  // sigma of various distributions
  constexpr double gk_phi_distribution_sigma_s1dc = 5.08;
  constexpr double gk_phi_distribution_sigma_fdc0 = 6.16;
  constexpr double gk_vertz_distribution_sigma = 16;
  constexpr double gk_rel_e_distribution_sigma = 0.2;
  constexpr double gk_theta_distribution_sigma_s1dc = 0.525;
  constexpr double gk_theta_distribution_sigma_fdc0 = 0.74;

  // constants for variable phi cut
  constexpr double gk_lin_phi_cut_fdc0_a = 0.18278;
  constexpr double gk_lin_phi_cut_fdc0_b = -6.45982;
  constexpr double gk_lin_phi_cut_s1dc_a = 0.157638;
  constexpr double gk_lin_phi_cut_s1dc_b = -5.79589;

  // constants for CS/AyPt evaluation
  constexpr int gk_cs_bin_number = 16;
  constexpr double gk_cs_range_start = 55.5;
  constexpr double gk_cs_range_end = 70;
}
