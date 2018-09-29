

#include <cassert>
#include <TMath.h>

#include "varthetacut.hpp"


static std::string
get_theta_dist_data_file_name(s13::ana::DatasetType ds_type,
                              s13::ana::Dsdcs dsdc_selector,
                              s13::ana::Espri espri_selector) {
  std::string data_file = "data/theta_sigmas";
  data_file += ds_type == s13::ana::DatasetType::he4 ? "_he4" : "_he6";
  data_file += dsdc_selector == s13::ana::Dsdcs::fdc0 ? "_fdc0" : "_s1dc";
  data_file += espri_selector == s13::ana::Espri::left ? "_esl" : "_esr";
  data_file += ".csv";
  return data_file;
}

std::shared_ptr<ROOT::Math::Interpolator>
s13::ana::make_theta_sigma_interp(DatasetType ds_type, Dsdcs dsdc_selector,
                                  Espri espri_selector) {
  assert(espri_selector != Espri::both);
  s13::misc::MessageLogger logger("make_theta_sigma_interp()");
  std::string data_file = get_theta_dist_data_file_name(ds_type, dsdc_selector, espri_selector);

  logger.info("Making interpolator for theta cut sigma "
              "parameter from data in %s...", data_file.c_str());
  io::CsvColumnsVector csv_data = s13::io::load_csv(data_file, 7, 2);
  std::vector<double>& p_angles = csv_data.at(0);
  std::vector<double>& theta_sigmas = csv_data.at(3);

  auto interp = new ROOT::Math::Interpolator{
    static_cast<unsigned int>(p_angles.size()),
    ROOT::Math::Interpolation::kLINEAR};
  interp->SetData(p_angles, theta_sigmas);

  using interp_ptr_t = std::unique_ptr<ROOT::Math::Interpolator>;
  return interp_ptr_t(interp);
}

std::shared_ptr<ROOT::Math::Interpolator>
s13::ana::make_theta_mean_interp(DatasetType ds_type, Dsdcs dsdc_selector,
                                 Espri espri_selector) {
  assert(espri_selector != Espri::both);
  s13::misc::MessageLogger logger("make_theta_sigma_interp()");
  std::string data_file = get_theta_dist_data_file_name(ds_type, dsdc_selector, espri_selector);

  logger.info("Making interpolator for theta cut mean "
              "parameter from data in %s...", data_file.c_str());
  io::CsvColumnsVector csv_data = s13::io::load_csv(data_file, 7, 2);
  std::vector<double>& p_angles = csv_data.at(0);
  std::vector<double>& theta_sigmas = csv_data.at(2);

  auto interp = new ROOT::Math::Interpolator{
    static_cast<unsigned int>(p_angles.size()),
    ROOT::Math::Interpolation::kLINEAR};
  interp->SetData(p_angles, theta_sigmas);

  using interp_ptr_t = std::unique_ptr<ROOT::Math::Interpolator>;
  return interp_ptr_t(interp);
}
