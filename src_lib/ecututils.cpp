

#include "ecututils.hpp"


std::shared_ptr<s13::ana::EDistParamsInterp>
s13::ana::make_he6_e_dist_sigma_interp() {
  std::string fname = "data/he6_e_dist_interp_sigmas.csv";
  std::shared_ptr<EDistParamsInterp> p_interp{new EDistParamsInterp(fname)};
  return p_interp;
}

std::shared_ptr<s13::ana::EDistParamsInterp>
s13::ana::make_he4_e_dist_sigma_interp() {
  std::string fname = "data/he4_e_dist_interp_sigmas.csv";
  std::shared_ptr<EDistParamsInterp> p_interp{new EDistParamsInterp(fname)};
  return p_interp;
}

std::shared_ptr<s13::ana::EDistParamsInterp>
s13::ana::make_he6_e_dist_mean_interp(s13::ana::Espri espri) {
  assert(espri != Espri::both);
  std::string espri_prefix = espri == Espri::left ? "esl" : "esr";
  std::string fname = "data/he6_" + espri_prefix + "_e_dist_interp_means.csv";
  std::shared_ptr<EDistParamsInterp> p_interp{new EDistParamsInterp(fname)};
  return p_interp;
}

std::shared_ptr<s13::ana::EDistParamsInterp>
s13::ana::make_he4_e_dist_mean_interp(s13::ana::Espri espri) {
  assert(espri != Espri::both);
  std::string espri_prefix = espri == Espri::left ? "esl" : "esr";
  std::string fname = "data/he4_" + espri_prefix + "_e_dist_interp_means.csv";
  std::shared_ptr<EDistParamsInterp> p_interp{new EDistParamsInterp(fname)};
  return p_interp;
}

std::shared_ptr<s13::ana::ERelMeanInterp>
s13::ana::make_he6_nai_erel_mean_interp() {
  std::string esl_fit_coeffs_fname = "data/he6_left_nai_e_rel_lin_fit_params.csv";
  std::string esr_fit_coeffs_fname = "data/he6_right_nai_e_rel_lin_fit_params.csv";
  auto esl_data = io::load_csv(esl_fit_coeffs_fname, 3, 2);
  auto esr_data = io::load_csv(esr_fit_coeffs_fname, 3, 2);
  auto interp = new ERelMeanInterp(esl_data, esr_data);
  std::shared_ptr<ERelMeanInterp> p_interp(interp);
  return p_interp;
}

std::shared_ptr<s13::ana::ERelMeanInterp>
s13::ana::make_he4_nai_erel_mean_interp() {
  std::string esl_fit_coeffs_fname = "data/he4_left_nai_e_rel_lin_fit_params.csv";
  std::string esr_fit_coeffs_fname = "data/he4_right_nai_e_rel_lin_fit_params.csv";
  auto esl_data = io::load_csv(esl_fit_coeffs_fname, 3, 2);
  auto esr_data = io::load_csv(esr_fit_coeffs_fname, 3, 2);
  auto interp = new ERelMeanInterp(esl_data, esr_data);
  std::shared_ptr<ERelMeanInterp> p_interp(interp);
  return p_interp;
}
