

#include <csutil.hpp>


double
s13::ana::lab2cm_solid_angle_factor(double proton_cm_theta_rad,
                                    DatasetType ds_type) {
  const double p_he4_gamma_cm = 1.12654467394258;
  const double p_he6_gamma_cm = 1.14966481938141;
  double gamma = p_he4_gamma_cm;

  if (ds_type == DatasetType::he6 || ds_type == DatasetType::carbon) {
    gamma = p_he6_gamma_cm;
  } else if (ds_type == DatasetType::he4) {
    gamma = p_he4_gamma_cm;
  } else {
    throw std::invalid_argument("Invalid dataset type!");
  }

  double theta = proton_cm_theta_rad;
  double num = gamma * std::abs(std::cos(theta) - 1);
  double denom = std::pow(std::sin(theta), 2) +
    std::pow(gamma, 2) * std::pow(1 - std::cos(theta), 2);
  denom = std::pow(denom, 3./2);
  return num / denom;
}

s13::ana::ScatteringYield
s13::ana::lab_yield_to_cm(const s13::ana::ScatteringYield& lab_yield,
                          DatasetType ds_type) {
  misc::MessageLogger log_("lab_yield_to_cm()");
  TInterpPtr angle_interp = make_lab2cm_angle_interp(ds_type);
  const double lab_theta_min =
    (lab_yield.GetBinArg(0) - lab_yield.GetBinWidth() / 2);
  const double lab_theta_max =
    (lab_yield.GetBinArg(lab_yield.GetBinsNumber() - 1)
     + lab_yield.GetBinWidth() / 2);
  // const double cm_theta_max = lab_theta_angle_to_cm(lab_theta_min);
  // const double cm_theta_min = lab_theta_angle_to_cm(lab_theta_max);
  const double cm_theta_max = angle_interp->Eval(lab_theta_min);
  const double cm_theta_min = angle_interp->Eval(lab_theta_max);

  log_.debug("CM angular range: [%.1f; %.1f]", cm_theta_min, cm_theta_max);
  ScatteringYield cm_yield(lab_yield.GetBinsNumber(), cm_theta_min, cm_theta_max);

  log_.debug("Converting yield...");
  for (int i = cm_yield.GetBinsNumber() - 1; i >= 0; --i) {
    // Double_t cm_angle = lab_theta_angle_to_cm(lab_yield.GetBinArg(i));
    Double_t cm_angle = angle_interp->Eval(lab_yield.GetBinArg(i));
    cm_yield.SetBinValue(cm_angle, lab_yield.GetBinValue(i));
    cm_yield.SetBinError(cm_angle, lab_yield.GetBinError(i));
  }

  return cm_yield;
}

s13::ana::TInterpPtr
s13::ana::make_lab2cm_angle_interp(DatasetType ds_type) {
  std::string fname = ds_type == DatasetType::he4
    ? "data/p_he4_recoil_A_lab_A_cm.csv" : "data/p_he6_recoil_A_lab_A_cm.csv";
  auto csv_data = io::load_csv(fname, 2, 2);
  assert (csv_data.size() == 2 && "Invalid number of columns in input file.");

  std::vector<double> a_cm, a_lab;
  for (size_t i = 0; i < csv_data[0].size(); i++) {
    a_lab.push_back(csv_data[0][i]);
    a_cm.push_back(csv_data[1][i]);
  }

  TInterp* interp = new TInterp{static_cast<unsigned int>(a_lab.size()),
                                ROOT::Math::Interpolation::kLINEAR};
  interp->SetData(a_lab, a_cm);
  return TInterpPtr{interp};
}

/*
  Converts CS from lab to center of mass frame.
*/
s13::ana::ScatteringYield
s13::ana::lab_cs_to_cm(const s13::ana::ScatteringYield& lab_cs,
                       DatasetType ds_type) {
  misc::MessageLogger log_("lab_cs_to_cm()");
  TInterpPtr angle_interp = make_lab2cm_angle_interp(ds_type);
  const double lab_theta_min =
    (lab_cs.GetBinArg(0) - lab_cs.GetBinWidth() / 2);
  const double lab_theta_max =
    (lab_cs.GetBinArg(lab_cs.GetBinsNumber() - 1)
     + lab_cs.GetBinWidth() / 2);
  // const double cm_theta_max = lab_theta_angle_to_cm(lab_theta_min);
  // const double cm_theta_min = lab_theta_angle_to_cm(lab_theta_max);
  const double cm_theta_max = angle_interp->Eval(lab_theta_min);
  const double cm_theta_min = angle_interp->Eval(lab_theta_max);

  log_.debug("CM angular range: [%.1f; %.1f]", cm_theta_min, cm_theta_max);
  ScatteringYield cm_cs(lab_cs.GetBinsNumber(), cm_theta_min, cm_theta_max);

  log_.debug("Converting cross-section...");
  for (int i = cm_cs.GetBinsNumber() - 1; i >= 0; --i) {
    int cm_cs_idx = cm_cs.GetBinsNumber() - 1 - i;
    // Double_t cm_angle = lab_theta_angle_to_cm(lab_cs.GetBinArg(i));
    Double_t cm_angle = angle_interp->Eval(lab_cs.GetBinArg(i));
    log_.debug("\t Computing CM MS for angle: %.3f; angle in scattering "
               "yield object: %.3f", cm_angle, cm_cs.GetBinArg(cm_cs_idx));

    Double_t solid_angle_factor_nonrel = 1 /
      (4 * TMath::Sin((cm_angle * gk_d2r) / 2));
    Double_t solid_angle_factor =
      lab2cm_solid_angle_factor(cm_angle * s13::gk_d2r, ds_type);
    log_.debug("SA factors (rel; non-rel): %.4f, %.4f",
               solid_angle_factor, solid_angle_factor_nonrel);
    Double_t cm_cs_value = lab_cs.GetBinValue(i) * solid_angle_factor;
    log_.debug("\t CS value in lab frame: %.3f; in CM frame: %.3f",
               lab_cs.GetBinValue(i), cm_cs_value);

    cm_cs.SetBinValue(cm_angle, cm_cs_value);
    cm_cs.SetBinError(cm_angle, lab_cs.GetBinError(i) * solid_angle_factor);
  }

  return cm_cs;
}

s13::ana::ScatteringYield
s13::ana::cs_from_yields(s13::ana::ScatteringYield& yields,
                         s13::ana::SolidAngle& solid_angle,
                         double num_beam_incident,
                         double tgt_areal_density) {
  ScatteringYield cs = yields;
  for (int i = 0; i < yields.GetBinsNumber(); i++) {
    double theta = yields.GetBinArg(i);
    double sa = solid_angle.at(theta);
    Double_t norm_factor =
      1 / (num_beam_incident * tgt_areal_density * sa);
    // from cm2 to barns
    norm_factor *= 1e24;
    // now into millibarns
    norm_factor *= 1e3;

    cs.ScaleBin(i, norm_factor);
  }
  return cs;
}
