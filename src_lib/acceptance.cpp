

#include "acceptance.hpp"


bool s13::ana::operator<(const PolarPoint& rhs, const PolarPoint& lhs) {
  if (rhs.theta < lhs.theta) {
    return true;
  } else if (rhs.theta > lhs.theta) {
    return false;
  }

  if (rhs.phi < lhs.phi) {
    return true;
  }

  return false;
}

bool s13::ana::operator==(const PolarPoint& rhs, const PolarPoint& lhs) {
  double epsilon = std::numeric_limits<double>::epsilon();
  if (std::abs(rhs.theta - lhs.theta) <= epsilon &&
      std::abs(rhs.phi - lhs.phi) <= epsilon) {
    return true;
  }
  return false;
}

bool s13::ana::operator<(const AbAngPoint& rhs, const AbAngPoint& lhs) {
  if (rhs.aang < lhs.aang) {
    return true;
  } else if (rhs.aang > lhs.aang) {
    return false;
  }

  if (rhs.bang < lhs.bang) {
    return true;
  }

  return false;
}

bool s13::ana::operator==(const AbAngPoint& rhs, const AbAngPoint& lhs) {
  double epsilon = std::numeric_limits<double>::epsilon();
  if (std::abs(rhs.aang - lhs.aang) <= epsilon &&
      std::abs(rhs.bang - lhs.bang) <= epsilon) {
    return true;
  }
  return false;
}

bool s13::ana::operator==(const AcceptanceElement& rhs,
                          const AcceptanceElement& lhs) {
  double epsilon = std::numeric_limits<double>::epsilon();
  if (std::abs(rhs.x1 - lhs.x1) <= epsilon &&
      std::abs(rhs.x2 - rhs.x2) <= epsilon) {
    return true;
  }
  return false;
}

bool s13::ana::operator<(const AcceptanceElement& rhs,
                         const AcceptanceElement& lhs) {
  if (rhs.x1 < lhs.x1) {
    return true;
  } else if (rhs.x1 > lhs.x1) {
    return false;
  }

  if (rhs.x2 < lhs.x2) {
    return true;
  }

  return false;
}

s13::ana::LabAcceptancePrefs
s13::ana::default_lab_acceptance_prefs() {
  return s13::ana::LabAcceptancePrefs();
}

s13::ana::AbAngAcceptancePrefs
s13::ana::default_abang_acceptance_prefs() {
  return s13::ana::AbAngAcceptancePrefs();
}

std::set<s13::ana::AbAngPoint>
s13::ana::lab_to_abang_acceptance(std::set<s13::ana::PolarPoint>& lab_acc) {
  using namespace s13::ana;
  s13::misc::MessageLogger logger("s13::ana::abang_acceptance()");
  std::set<AbAngPoint> abang_acc;
  std::for_each(lab_acc.begin(), lab_acc.end(),
                [&abang_acc, &lab_acc, &logger] (const PolarPoint& pt) {
                  TVector3 vec;
                  vec.SetMagThetaPhi(1,
                                     pt.theta * s13::gk_d2r,
                                     (pt.phi + 90) * s13::gk_d2r);
                  AbAngPoint ab_pt(TMath::ATan2(vec.X(), vec.Z()) / s13::gk_d2r,
                                   TMath::ATan2(vec.Y(), vec.Z()) / s13::gk_d2r);
                  // AbAngPoint ab_pt((vec.X() / vec.Z()) / s13::gk_d2r,
                  //                  (vec.Y() / vec.X()) / s13::gk_d2r);
                  logger.debug("Adding AbAngPoint with angles: %.1f, %.1f",
                               ab_pt.aang, ab_pt.bang);
                  abang_acc.insert(ab_pt);
                });
  return abang_acc;
}

/*
  Converts CM theta into lab theta of recoil particle (lab theta2);
*/
static double cm_theta_to_lab_theta2(double cm_theta) {
  return ((TMath::Pi() - (cm_theta * s13::gk_d2r)) / 2) / s13::gk_d2r;
}

/*
  Converts CM theta to lab theta of scattered fragment (lab theta1)
*/
static double cm_theta_to_lab_theta1(double cm_theta) {
  double cm_theta_rad = cm_theta * s13::gk_d2r;
  double tan_lab_he_angle = TMath::Sin(cm_theta_rad);
  tan_lab_he_angle /= TMath::Cos(cm_theta_rad) + 6./1.;
  double lab_he_angle_rad = TMath::ATan(tan_lab_he_angle);
  return lab_he_angle_rad / s13::gk_d2r;
}

static TVector3 projection_yz(TVector3 vec) {
  vec.SetX(0);
  return vec;
}

static TVector3 projection_xz(TVector3 vec) {
  vec.SetY(0);
  return vec;
}

static TVector3 get_direction_vector(s13::ana::AcceptanceTypes type,
                                     double x1, double x2) {
  s13::misc::MessageLogger logger("get_direction_vector()");
  TVector3 direction(0, 0, 1);
  if (type == s13::ana::AcceptanceTypes::lab_ab) {
    // logger.trace("Getting lab_ab direction vector for angles: %.1f, %.1f", x1, x2);
    direction.RotateX(x2 * s13::gk_d2r);
    direction.RotateY(x1 * s13::gk_d2r);
  } else if (type == s13::ana::AcceptanceTypes::lab_theta_phi) {
    direction.SetMagThetaPhi(1, x1 * s13::gk_d2r, (x2 + 90) * s13::gk_d2r);
  } else if (type == s13::ana::AcceptanceTypes::cm_theta_phi_recoil) {
    double lab_theta2 = cm_theta_to_lab_theta2(x1);
    // logger.trace("Converting CM theta to lab theta2 (CM -> lab): (%.1f -> %.1f)",
    //              x1, lab_theta2);
    return get_direction_vector(s13::ana::AcceptanceTypes::lab_theta_phi,
                                lab_theta2, x2);
  } else if (type == s13::ana::AcceptanceTypes::cm_theta_phi_fragment) {
    double lab_theta1 = cm_theta_to_lab_theta1(x1);
    return get_direction_vector(s13::ana::AcceptanceTypes::lab_theta_phi,
                                lab_theta1, x2);
  }
  return direction;
}

s13::ana::AcceptanceSet
s13::ana::compute_acceptance(TVector3 vertex_pos,
                             s13::ana::GeometricalAcceptance& detector_acceptance,
                             s13::ana::AcceptanceTypes type,
                             s13::ana::AcceptancePrefs prefs) {
  s13::misc::MessageLogger logger("s13::ana::compute_acceptance()");
  s13::ana::AcceptanceSet res;

  for (double x1 = prefs.min_x1;
       x1 <= prefs.max_x1; x1 += prefs.step_x1) {
    for (double x2 = prefs.min_x2;
         x2 <= prefs.max_x2; x2 += prefs.step_x2) {
      TVector3 direction = get_direction_vector(type, x1, x2);
      if (detector_acceptance.in_lab_acceptance(vertex_pos, direction)) {
        // logger.trace("Direction vector with angles (x1, x2): "
        //              "(%.1f, %.1f) in acceptance.", x1, x2);
          s13::ana::AcceptanceElement el(x1, x2);
        res.insert(el);
      }
    }
  }

  return res;
}

std::set<s13::ana::PolarPoint>
s13::ana::lab_acceptance(TVector3 vertex_pos,
                         s13::ana::GeometricalAcceptance& detector_acceptance,
                         LabAcceptancePrefs* prefs) {
  s13::misc::MessageLogger("s13::ana::lab_acceptance()");
  LabAcceptancePrefs acc_prefs = default_lab_acceptance_prefs();
  if (prefs != nullptr) {
    acc_prefs = *prefs;
  }

  std::set<s13::ana::PolarPoint> res;

  for (double theta = acc_prefs.min_theta;
       theta <= acc_prefs.max_theta; theta += acc_prefs.step_theta) {
    for (double x2 = acc_prefs.min_phi;
         x2 <= acc_prefs.max_phi; x2 += acc_prefs.step_phi) {
      TVector3 direction(0, 0, 0);
      direction.SetMagThetaPhi(1, theta * s13::gk_d2r, (x2 + 90) * s13::gk_d2r);
      if (detector_acceptance.in_lab_acceptance(vertex_pos, direction)) {
        s13::ana::PolarPoint pt(theta, x2);
        res.insert(pt);
      }
    }
  }

  return res;
}

std::set<s13::ana::AbAngPoint>
s13::ana::abang_acceptance(TVector3 vertex_pos,
                           s13::ana::GeometricalAcceptance& detector_acceptance,
                           AbAngAcceptancePrefs* prefs) {
  using namespace s13::ana;
  s13::misc::MessageLogger logger("s13::ana::abang_acceptance()");
  AbAngAcceptancePrefs acc_prefs = default_abang_acceptance_prefs();
  if (prefs != nullptr) {
    acc_prefs = *prefs;
  }

  std::set<AbAngPoint> res;

  for (double aang = acc_prefs.min_aang;
       aang <= acc_prefs.max_aang; aang += acc_prefs.step_aang) {
    for (double bang = acc_prefs.min_bang;
         bang <= acc_prefs.max_bang; bang += acc_prefs.step_bang) {
      TVector3 direction(0, 0, 1);
      // direction.SetXYZ(TMath::Tan(aang * gk_d2r), TMath::Tan(bang * gk_d2r), 1);
      direction.RotateX(bang * gk_d2r);
      direction.RotateY(aang * gk_d2r);
      if (detector_acceptance.in_lab_acceptance(vertex_pos, direction)) {
        AbAngPoint pt(aang, bang);
        res.insert(pt);
      }
    }
  }

  return res;
}


///
/// Functions to compute specific detector acceptances
///

TH2*
s13::ana::make_upstream_lab_acceptance_hist(const s13::ana::AcceptanceSet& acc,
                                            double bin_size,
                                            TString title,
                                            Espri espri_selector) {
  assert(espri_selector != Espri::both);
  double constexpr min_theta = 25;
  double constexpr max_theta = 85;
  double min_phi = 50;
  double max_phi = 130;
  if (espri_selector == Espri::right) {
    min_phi = -130;
    max_phi = -50;
  }

  int nbin_x = (max_theta - min_theta) / bin_size;
  int nbin_y = (max_phi - min_phi) / bin_size;
  TH2* hist_lab = s13::ana::make_th2(nbin_x, min_theta, max_theta,
                                     nbin_y, min_phi, max_phi,
                                     title,
                                     "#theta [lab. deg.]", "#x2 [lab.deg.]");

  for (auto& el : acc) {
    hist_lab->Fill(el.x1, el.x2);
  }

  return hist_lab;
}

TH2*
s13::ana::make_upstream_lababang_acceptance_hist(const s13::ana::AcceptanceSet& acc,
                                                 double bin_size,
                                                 TString title,
                                                 Espri espri_selector) {
  double min_aang = 25;
  double max_aang = 85;
  double min_bang = -40;
  double max_bang = 40;
  if (espri_selector == Espri::right) {
    min_aang = -85;
    max_aang = -25;
  }

  int nbin_x = (max_aang - min_aang) / bin_size;
  int nbin_y = (max_bang - min_bang) / bin_size;
  TH2* hist_lab =
    s13::ana::make_th2(nbin_x, min_aang, max_aang,
                       nbin_y, min_bang, max_bang,
                       title,
                       "aang [lab. deg.]", "#bang [lab.deg.]");

  for (auto& el : acc) {
    hist_lab->Fill(el.x1, el.x2);
  }

  return hist_lab;
}

TH2*
s13::ana::make_upstream_cm_acceptance_hist(const s13::ana::AcceptanceSet& acc,
                                           double bin_size,
                                           TString title,
                                           Espri espri_selector) {
  double constexpr min_theta = 10;
  double constexpr max_theta = 130;
  double min_phi = 50;
  double max_phi = 130;
  // the switch to select FDC2 acceptance area
  // which corresponds to ESPRI left/right
  // ESPRI right - left scattered He
  // ESPRI left - right scattered He
  if (espri_selector == Espri::right) {
    min_phi = -140;
    max_phi = -50;
  }

  int nbin_x = (max_theta - min_theta) / bin_size;
  int nbin_y = (max_phi - min_phi) / bin_size;
  TH2* hist_cm = s13::ana::make_th2(nbin_x, min_theta, max_theta,
                                     nbin_y, min_phi, max_phi,
                                     title,
                                     "#theta [lab. deg.]", "#x2 [lab.deg.]");

  for (auto& el : acc) {
    hist_cm->Fill(el.x1, el.x2);
  }

  return hist_cm;
}

TH2*
s13::ana::make_upstream_acceptance_hist(s13::ana::AcceptanceTypes type,
                                        const s13::ana::AcceptanceSet& acc_set,
                                        double bin_size,
                                        TString title,
                                        Espri espri_selector) {
  if (type == s13::ana::AcceptanceTypes::lab_ab) {
    title = s13::ana::tstrfmt("%s (lab a & b angles)", title.Data());
    return make_upstream_lababang_acceptance_hist(acc_set, bin_size,
                                                  title, espri_selector);
  } else if (type == s13::ana::AcceptanceTypes::lab_theta_phi) {
    title = s13::ana::tstrfmt("%s (lab #theta & #phi angles)", title.Data());
    return make_upstream_lab_acceptance_hist(acc_set, bin_size,
                                             title, espri_selector);
  } else if (type == s13::ana::AcceptanceTypes::cm_theta_phi_recoil) {
    title = s13::ana::tstrfmt("%s (CM #theta & #phi angles)", title.Data());
    return make_upstream_cm_acceptance_hist(acc_set, bin_size,
                                            title, espri_selector);
  } else if (type == s13::ana::AcceptanceTypes::cm_theta_phi_fragment) {
    title = s13::ana::tstrfmt("%s (CM #theta & #phi angles)", title.Data());
    return make_upstream_cm_acceptance_hist(acc_set, bin_size,
                                            title, espri_selector);
  } else {
    throw std::invalid_argument("Unsupported coordinate type for this function!");
  }
}

TH2*
s13::ana::make_downstream_lab_acceptance_hist(const s13::ana::AcceptanceSet& acc,
                                              double bin_size,
                                              TString title) {
  double constexpr min_theta = 0;
  double constexpr max_theta = 13;
  double constexpr min_phi = -180;
  double constexpr max_phi = 180;

  int nbin_x = (max_theta - min_theta) / bin_size;
  int nbin_y = (max_phi - min_phi) / bin_size;
  TH2* hist_lab = s13::ana::make_th2(nbin_x, min_theta, max_theta,
                                     nbin_y, min_phi, max_phi,
                                     title,
                                     "#theta [lab. deg.]", "#x2 [lab.deg.]");

  for (auto& el : acc) {
    hist_lab->Fill(el.x1, el.x2);
  }

  return hist_lab;
}

TH2*
s13::ana::make_downstream_lababang_acceptance_hist(const s13::ana::AcceptanceSet& acc,
                                                   double bin_size,
                                                   TString title) {
  double constexpr min_aang = -15;
  double constexpr max_aang = 15;
  double constexpr min_bang = -12;
  double constexpr max_bang = 12;

  int nbin_x = (max_aang - min_aang) / bin_size;
  int nbin_y = (max_bang - min_bang) / bin_size;
  TH2* hist_lab =
    s13::ana::make_th2(nbin_x, min_aang, max_aang,
                       nbin_y, min_bang, max_bang,
                       title,
                       "aang [lab. deg.]", "#bang [lab.deg.]");

  for (auto& el : acc) {
    hist_lab->Fill(el.x1, el.x2);
  }

  return hist_lab;
}

TH2*
s13::ana::make_downstream_cm_acceptance_hist(const s13::ana::AcceptanceSet& acc,
                                             double bin_size,
                                             TString title) {
  double constexpr min_theta = 10;
  double constexpr max_theta = 130;
  double constexpr min_phi = -180;
  double constexpr max_phi = 180;

  int nbin_x = (max_theta - min_theta) / bin_size;
  int nbin_y = (max_phi - min_phi) / bin_size;
  TH2* hist_cm = s13::ana::make_th2(nbin_x, min_theta, max_theta,
                                     nbin_y, min_phi, max_phi,
                                     title,
                                     "#theta [lab. deg.]", "#x2 [lab.deg.]");

  for (auto& el : acc) {
    hist_cm->Fill(el.x1, el.x2);
  }

  return hist_cm;
}

TH2*
s13::ana::make_downstream_acceptance_hist(s13::ana::AcceptanceTypes type,
                                          const s13::ana::AcceptanceSet& acc_set,
                                          double bin_size,
                                          TString title) {
  if (type == s13::ana::AcceptanceTypes::lab_ab) {
    title = s13::ana::tstrfmt("%s (lab a & b angles)", title.Data());
    return make_downstream_lababang_acceptance_hist(acc_set, bin_size, title);
  } else if (type == s13::ana::AcceptanceTypes::lab_theta_phi) {
    title = s13::ana::tstrfmt("%s (lab #theta & #phi angles)", title.Data());
    return make_downstream_lab_acceptance_hist(acc_set, bin_size, title);
  } else if (type == AcceptanceTypes::cm_theta_phi_fragment
             || type == AcceptanceTypes::cm_theta_phi_recoil) {
    title = s13::ana::tstrfmt("%s (CM #theta & #phi angles)", title.Data());
    return make_downstream_cm_acceptance_hist(acc_set, bin_size, title);
  } else {
    throw std::invalid_argument("Unsupported coordinate type for this function!");
  }
}

s13::ana::AcceptanceSet
s13::ana::compute_tgt_wnd_acceptance(s13::ana::AcceptanceTypes type,
                                     double vertexY,
                                     double step) {
  // window displacement from origin along Z axis in mm (from CAD drawings)
  constexpr double k_wnd_disp_z = 10;

  s13::misc::MessageLogger logger("s13::ana::compute_tgt_wnd_acc()");
  s13::ana::TargetWndAcceptance wnd_acc(k_wnd_disp_z);
  s13::ana::AcceptancePrefs acc_prefs;
  s13::ana::AbAngAcceptancePrefs acc_prefs_abang;

  acc_prefs.step_x1 = step;
  acc_prefs.step_x2 = step;

  if (type == s13::ana::AcceptanceTypes::lab_ab) {
    acc_prefs.min_x1 = 0;
    acc_prefs.max_x1 = 90;
  } else if (type == s13::ana::AcceptanceTypes::lab_theta_phi) {
    acc_prefs.min_x1 = 30;
    acc_prefs.max_x1 = 80;
  } else if (type == s13::ana::AcceptanceTypes::cm_theta_phi_recoil) {
    acc_prefs.min_x1 = 10;
    acc_prefs.max_x1 = 120;
  } else  {
    throw std::invalid_argument("Unsupported type of coordinates for this function!");
  }

  if (type == s13::ana::AcceptanceTypes::lab_ab) {
    acc_prefs.min_x2 = -40;
    acc_prefs.max_x2 = 0;
  } else if (type == s13::ana::AcceptanceTypes::lab_theta_phi) {
    acc_prefs.min_x2 = -90;
    acc_prefs.max_x2 = 0;
  } else if (type == s13::ana::AcceptanceTypes::cm_theta_phi_recoil) {
    acc_prefs.min_x2 = -90;
    acc_prefs.max_x2 = 0;
  }
  // Y-max vertex position
  auto y_max_edge_acc = compute_acceptance(TVector3(0, vertexY, 0),
                                           wnd_acc, type, acc_prefs);

  if (type == s13::ana::AcceptanceTypes::lab_ab) {
    acc_prefs.min_x2 = 0;
    acc_prefs.max_x2 = 40;
  } else if (type == s13::ana::AcceptanceTypes::lab_theta_phi) {
    acc_prefs.min_x2 = -180;
    acc_prefs.max_x2 = -90;
  } else if (type == s13::ana::AcceptanceTypes::cm_theta_phi_recoil) {
    acc_prefs.min_x2 = -180;
    acc_prefs.max_x2 = -90;
  }
  // Y-min vertex position
  auto y_min_edge_acc = compute_acceptance(TVector3(0, -vertexY, 0),
                                           wnd_acc, type, acc_prefs);

  AcceptanceSet tgt_y_edge_acc;
  tgt_y_edge_acc.insert(y_max_edge_acc.begin(), y_max_edge_acc.end());
  tgt_y_edge_acc.insert(y_min_edge_acc.begin(), y_min_edge_acc.end());
  logger.debug("Target window acceptance computed, size of acceptance set: %d",
               tgt_y_edge_acc.size());
  return tgt_y_edge_acc;
}

s13::ana::AcceptanceSet
s13::ana::compute_rdc_acceptance(s13::ana::AcceptanceTypes type,
                                 double step) {
  s13::misc::MessageLogger logger("s13::ana::compute_rdc_acc()");
  s13::ana::RdcAcceptance rdc_acc(s13::gk_dist_rdc_target,
                                  s13::gk_rdc_center_angle * s13::gk_d2r);
  s13::ana::AcceptancePrefs acc_prefs;
  s13::ana::AbAngAcceptancePrefs acc_prefs_abang;

  acc_prefs.step_x1 = step;
  acc_prefs.step_x2 = step;

  if (type == s13::ana::AcceptanceTypes::lab_ab) {
    acc_prefs.min_x1 = 0;
    acc_prefs.max_x1 = 90;
    acc_prefs.min_x2 = -40;
    acc_prefs.max_x2 = 40;
  } else if (type == s13::ana::AcceptanceTypes::lab_theta_phi) {
    acc_prefs.min_x1 = 30;
    acc_prefs.max_x1 = 80;
    acc_prefs.min_x2 = -130;
    acc_prefs.max_x2 = -50;
  } else if (type == s13::ana::AcceptanceTypes::cm_theta_phi_recoil) {
    acc_prefs.min_x1 = 10;
    acc_prefs.max_x1 = 120;
    acc_prefs.min_x2 = -130;
    acc_prefs.max_x2 = -50;
  } else  {
    throw std::invalid_argument("Unsupported type of coordinates for this function!");
  }
  auto acc = compute_acceptance(TVector3(0, 0, 0),
                                rdc_acc, type, acc_prefs);
  return acc;
}

s13::ana::AcceptanceSet
s13::ana::compute_nai_acceptance(s13::ana::AcceptanceTypes type,
                                 double step) {
  s13::misc::MessageLogger logger("s13::ana::compute_nai_acc()");
  s13::ana::NaiAcceptance nai_acc(s13::gk_dist_rdc_target,
                                  s13::gk_rdc_center_angle * s13::gk_d2r);
  s13::ana::AcceptancePrefs acc_prefs;

  acc_prefs.step_x1 = step;
  acc_prefs.step_x2 = step;

  if (type == s13::ana::AcceptanceTypes::lab_ab) {
    acc_prefs.min_x1 = 0;
    acc_prefs.max_x1 = 90;
    acc_prefs.min_x2 = -40;
    acc_prefs.max_x2 = 40;
  } else if (type == s13::ana::AcceptanceTypes::lab_theta_phi) {
    acc_prefs.min_x1 = 30;
    acc_prefs.max_x1 = 80;
    acc_prefs.min_x2 = -130;
    acc_prefs.max_x2 = -50;
  } else if (type == s13::ana::AcceptanceTypes::cm_theta_phi_recoil) {
    acc_prefs.min_x1 = 10;
    acc_prefs.max_x1 = 120;
    acc_prefs.min_x2 = -130;
    acc_prefs.max_x2 = -50;
  } else  {
    throw std::invalid_argument("Unsupported type of coordinates for this function!");
  }
  auto acc = compute_acceptance(TVector3(0, 0, 0),
                                nai_acc, type, acc_prefs);
  return acc;
}

s13::ana::AcceptanceSet
s13::ana::compute_fdc2_acceptance(s13::ana::AcceptanceTypes type,
                                  double step) {
  s13::misc::MessageLogger logger("s13::ana::compute_fdc2_acceptance()");
  s13::ana::AcceptancePrefs acc_prefs;

  if (type == s13::ana::AcceptanceTypes::lab_ab) {
    acc_prefs.min_x1 = -9.99;
    acc_prefs.max_x1 = 9.99;
    acc_prefs.min_x2 = -8;
    acc_prefs.max_x2 = 8;
  } else if (type == s13::ana::AcceptanceTypes::lab_theta_phi) {
    acc_prefs.min_x1 = 0;
    acc_prefs.max_x1 = 10;
    acc_prefs.min_x2 = -180;
    acc_prefs.max_x2 = 180;
  } else if (type == s13::ana::AcceptanceTypes::cm_theta_phi_fragment) {
    acc_prefs.min_x1 = 10;
    acc_prefs.max_x1 = 120;
    acc_prefs.min_x2 = -180;
    acc_prefs.max_x2 = 180;
  } else  {
    throw std::invalid_argument("Unsupported type of coordinates for this function!");
  }

  s13::ana::AcceptanceSet res;
  s13::ana::HeAbangFdc2YConverter conv("data/heabang_fdc2y_conv_factors.csv");
  for (double x1 = acc_prefs.min_x1; x1 < acc_prefs.max_x1; x1 += acc_prefs.step_x1) {
    for (double x2 = acc_prefs.min_x2; x2 < acc_prefs.max_x2; x2 += acc_prefs.step_x2) {
      TVector3 direction = get_direction_vector(type, x1, x2);
      TVector3 z_dir(0, 0, 1);
      double aang = projection_xz(direction).Angle(z_dir) / s13::gk_d2r;
      double bang = projection_yz(direction).Angle(z_dir) / s13::gk_d2r;
      double phi = direction.Phi() / s13::gk_d2r;

      if (phi < 0) { bang = -bang; }
      if ((phi > -180 && phi < -90) || (phi > 90 && phi < 180)) {
        aang = -aang;
      }

      double ypos = conv.fdc2_ypos(aang, bang);
      if (std::abs(ypos) < s13::gk_fdc2_acceptance_height/2) {
        double s13_phi = phi - 90;
        if (s13_phi < -180) {
          s13_phi += 360;
        }
        res.insert(AcceptanceElement(x1, x2));
      }

    }
  }

  return res;
}

TH2*
s13::ana::compute_tgt_wnd_acceptance_hist(s13::ana::AcceptanceTypes type,
                                          double vertexY,
                                          double step) {
  s13::ana::AcceptanceSet acc_set =
    compute_tgt_wnd_acceptance(type, vertexY, step);
  return make_downstream_acceptance_hist(type, acc_set, step, "Target window acceptance");
}

TH2*
s13::ana::compute_rdc_acceptance_hist(s13::ana::AcceptanceTypes type,
                                      double step) {
  s13::ana::AcceptanceSet acc_set =
    compute_rdc_acceptance(type, step);
  return make_downstream_acceptance_hist(type, acc_set, step, "RDC acceptance");
}

TH2*
s13::ana::compute_nai_acceptance_hist(s13::ana::AcceptanceTypes type,
                                      double step) {
  s13::ana::AcceptanceSet acc_set =
    compute_nai_acceptance(type, step);
  return make_downstream_acceptance_hist(type, acc_set, step, "NaIs acceptance");
}

TH2*
s13::ana::compute_fdc2_acceptance_hist(AcceptanceTypes type,
                                       double step) {
  AcceptanceSet acc_set = compute_fdc2_acceptance(type, step);
  return make_downstream_acceptance_hist(type, acc_set, step, "FDC2 acceptance");
}
