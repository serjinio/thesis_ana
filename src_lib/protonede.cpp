/**
 *   \file protonede.cpp
 *   \brief Implementation of proton dE-E computation.
 *
 */


#include "protonede.hpp"


std::unique_ptr<ROOT::Math::Interpolator>
s13::ana::load_he_kin_e(std::string filename) {
  std::ifstream ifile(filename);
  std::vector<Double_t> p_angles, p_energies;

  for (CSVIterator iter(ifile, 2); iter != CSVIterator(); iter++) {
    auto csv_row = *iter;
    assert(csv_row.size() == 4 && "Invalid input file format");

    Double_t p_a = std::stod(csv_row[2]) * D2R;
    Double_t p_e = std::stod(csv_row[3]);

    p_angles.push_back(p_a);
    p_energies.push_back(p_e);
  }

  auto interp = new ROOT::Math::Interpolator {
    static_cast<unsigned int>(p_angles.size()),
    ROOT::Math::Interpolation::kLINEAR};
  interp->SetData(p_angles, p_energies);

  using interp_ptr_t = std::unique_ptr<ROOT::Math::Interpolator>;
  return interp_ptr_t(interp);
}

std::unique_ptr<ROOT::Math::Interpolator>
s13::ana::load_he_kin_e2(std::string filename) {
  std::ifstream ifile(filename);
  std::vector<Double_t> p_angles, p_energies;

  for (CSVIterator iter(ifile, 2); iter != CSVIterator(); iter++) {
    auto csv_row = *iter;
    assert(csv_row.size() == 2 && "Invalid input file format");

    Double_t p_a = std::stod(csv_row[0]) * D2R;
    Double_t p_e = std::stod(csv_row[1]);

    p_angles.push_back(p_a);
    p_energies.push_back(p_e);
  }

  auto interp = new ROOT::Math::Interpolator {
    static_cast<unsigned int>(p_angles.size()),
    ROOT::Math::Interpolation::kLINEAR};
  interp->SetData(p_angles, p_energies);

  using interp_ptr_t = std::unique_ptr<ROOT::Math::Interpolator>;
  return interp_ptr_t(interp);
}

std::unique_ptr<ROOT::Math::Interpolator>
s13::ana::load_stopping_powers(std::string filename) {
  std::ifstream ifile(filename);
  std::vector<Double_t> p_energies, stopping_powers;

  for (CSVIterator iter(ifile, 2); iter != CSVIterator(); iter++) {
    auto csv_row = *iter;
    assert(csv_row.size() == 4 && "Invalid input file format");

    Double_t p_e = std::stod(csv_row[0]);
    Double_t stp = std::stod(csv_row[3]);

    p_energies.push_back(p_e);
    stopping_powers.push_back(stp);
  }

  auto interp = new ROOT::Math::Interpolator{
    static_cast<unsigned int>(p_energies.size()),
    ROOT::Math::Interpolation::kLINEAR};
  interp->SetData(p_energies, stopping_powers);

  using interp_ptr_t = std::unique_ptr<ROOT::Math::Interpolator>;
  return interp_ptr_t(interp);
}

Double_t
s13::ana::compute_proton_pla_de(Double_t proton_e,
                                std::unique_ptr<ROOT::Math::Interpolator>& stp_powers) {
  misc::MessageLogger logger("s13::ana::compute_proton_pla_de()");
  static const Double_t dx_step = 0.01;
  static const Double_t pla_thickness = 4;

  Double_t p_de = 0;
  Double_t p_e_cur = proton_e;
  for (Double_t z = 0; z < pla_thickness; z += dx_step) {
    if (p_e_cur < 1.01) {
      break;
    }
    Double_t de = dx_step * stp_powers->Eval(p_e_cur);
    p_de += de;
    p_e_cur -= de;
  }

  return p_de;
}

/*
  Returns two interpolator objects to obtain dE & E as a function of proton angle.
*/
std::pair<
  std::shared_ptr<ROOT::Math::Interpolator>,
  std::shared_ptr<ROOT::Math::Interpolator> >
s13::ana::compute_he6_proton_ref_de_e() {
  misc::MessageLogger logger("s13::ana::compute_he6_proton_ref_de_e()");
  std::string p_in_nai_stp_powers_file = "data/p_in_BC400.csv";
  std::string he6_kin_file = "data/he6_kin_E.csv";
  auto p_stp_powers = load_stopping_powers(p_in_nai_stp_powers_file);
  auto p_energies = load_he_kin_e(he6_kin_file);

  std::vector<Double_t> angles, des, es;
  for (Double_t p_angle = 50; p_angle < 80; p_angle += 0.5) {
    Double_t angle_rad = p_angle * D2R;
    Double_t p_e_kin = p_energies->Eval(angle_rad);
    Double_t p_de = compute_proton_pla_de(p_e_kin, p_stp_powers);
    angles.push_back(p_angle);
    des.push_back(p_de);
    es.push_back(p_e_kin);
    logger.trace("Evaluted proton E & dE @ %.1f deg.: %.1f; %.1f MeV",
                 p_angle, p_e_kin, p_de);
  }

  using interp_ptr_t = std::shared_ptr<ROOT::Math::Interpolator>;

  auto interp_de = new ROOT::Math::Interpolator{
    static_cast<unsigned int>(angles.size()),
    ROOT::Math::Interpolation::kLINEAR};
  interp_de->SetData(angles, des);
  auto p_interp_de = interp_ptr_t(interp_de);

  auto interp_e = new ROOT::Math::Interpolator{
    static_cast<unsigned int>(angles.size()),
    ROOT::Math::Interpolation::kLINEAR};
  interp_e->SetData(angles, es);
  auto p_interp_e = interp_ptr_t(interp_e);

  return make_pair(p_interp_de, p_interp_e);
  }

TGraph*
s13::ana::make_de_e_graph(std::shared_ptr<ROOT::Math::Interpolator>& angle_de,
                          std::shared_ptr<ROOT::Math::Interpolator>& angle_e) {
  std::vector<Double_t> des, es;
  Double_t *p_des, *p_es;

  for (Double_t angle = 50; angle <= 75; ++angle) {
    Double_t p_de = angle_de->Eval(angle);
    Double_t p_e = angle_e->Eval(angle);
    des.push_back(p_de);
    es.push_back(p_e);
  }

  p_des = new Double_t[des.size()];
  p_es = new Double_t[des.size()];
  std::copy(des.begin(), des.end(), p_des);
  std::copy(es.begin(), es.end(), p_es);

  TGraph* graph = new TGraph(des.size(), p_des, p_es);
  graph->SetLineColor(kRed);
  graph->SetLineWidth(4);
  graph->SetMarkerColor(kRed);
  graph->SetMarkerStyle(21);
  graph->SetTitle("dE-E");
  graph->GetXaxis()->SetTitle("dE [MeV]");
  graph->GetYaxis()->SetTitle("E [MeV]");
  return graph;
}

TGraph*
s13::ana::make_he6_proton_de_e_graph() {
  auto ref_angle_de_e_pair = compute_he6_proton_ref_de_e();
  auto ref_angle_de = std::get<0>(ref_angle_de_e_pair);
  auto ref_angle_e = std::get<1>(ref_angle_de_e_pair);
  auto graph = make_de_e_graph(ref_angle_de, ref_angle_e);
  return graph;
}

TGraph*
s13::ana::make_proton_angle_energy_graph(DatasetType ds_type,
                                         Espri espri) {
  const std::string he4_esr_fname = "he4_esr_kin_E_with_mtrls_and_E-deg_loss.csv";
  const std::string he4_esl_fname = "he4_esl_kin_E_with_mtrls_and_E-deg_loss.csv";
  const std::string he6_esr_fname = "he6_esr_kin_E_with_mtrls_and_E-deg_loss.csv";
  const std::string he6_esl_fname = "he6_esl_kin_E_with_mtrls_and_E-deg_loss.csv";
  std::string kin_fname = "data/";
  if (ds_type == DatasetType::he6) {
    kin_fname += espri == Espri::left ? he6_esl_fname : he6_esr_fname;
  } else if (ds_type == DatasetType::he4) {
    kin_fname += espri == Espri::left ? he4_esl_fname : he4_esr_fname;
  } else {
    throw std::invalid_argument("Invalid dataset type for this graph!");
  }

  auto angle_energy_interp = load_he_kin_e(kin_fname);
  std::vector<Double_t> angles, energies;
  Double_t *p_angles, *p_energies;

  for (Double_t angle = 53.50; angle < 73.00; ++angle) {
    Double_t p_e = angle_energy_interp->Eval(angle * D2R);
    angles.push_back(angle);
    energies.push_back(p_e);
  }

  p_angles = new Double_t[angles.size()];
  p_energies = new Double_t[angles.size()];
  std::copy(angles.begin(), angles.end(), p_angles);
  std::copy(energies.begin(), energies.end(), p_energies);

  TGraph* graph = new TGraph(angles.size(), p_angles, p_energies);
  graph->SetLineColor(kRed);
  graph->SetLineWidth(4);
  graph->SetMarkerColor(kRed);
  graph->SetMarkerStyle(21);
  graph->SetTitle("p-He4 kin. E");
  graph->GetXaxis()->SetTitle("Lab angle [deg]");
  graph->GetYaxis()->SetTitle("E [MeV]");
  return graph;
}
