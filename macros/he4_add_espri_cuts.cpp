
#include <fstream>

#include "TChain.h"
#include "TGraph.h"
#include "Math/Interpolator.h"

#include "init_4he_ds.hpp"
#include "macroincludes.hpp"


using namespace s13::ana;


std::unique_ptr<ROOT::Math::Interpolator>
load_he6_kin_e(std::string filename) {
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
load_stopping_powers(std::string filename) {
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
compute_proton_pla_de(Double_t proton_e,
                      std::unique_ptr<ROOT::Math::Interpolator>& stp_powers) {
  static const Double_t dx_step = 0.01;
  static const Double_t pla_thickness = 4;

  Double_t p_de = 0;
  Double_t p_e_cur = proton_e;
  for (Double_t z = 0; z < pla_thickness; z += dx_step) {
    Double_t de = dx_step * stp_powers->Eval(p_e_cur);
    p_de += de;
    p_e_cur -= de;
    if (p_e_cur < 1.0) {
      break;
    }
  }

  return p_de;
}

/*
  Returns two interpolator objects to obtain dE & E as a function of proton angle.
*/
std::pair<
  std::shared_ptr<ROOT::Math::Interpolator>,
  std::shared_ptr<ROOT::Math::Interpolator> >
compute_ref_de_e() {
  std::string p_in_nai_stp_powers_file = "data/p_in_BC400.csv";
  std::string he6_kin_file = "data/he4_kin_E.csv";
  auto p_stp_powers = load_stopping_powers(p_in_nai_stp_powers_file);
  auto p_energies = load_he6_kin_e(he6_kin_file);

  std::vector<Double_t> angles, des, es;
  for (Double_t p_angle = 52; p_angle < 75; p_angle += 0.5) {
    Double_t angle_rad = p_angle * D2R;
    Double_t p_e_kin = p_energies->Eval(angle_rad);
    Double_t p_de = compute_proton_pla_de(p_e_kin, p_stp_powers);
    Double_t p_e = p_e_kin - p_de;
    angles.push_back(p_angle);
    des.push_back(p_de);
    es.push_back(p_e);
//    std::cout << "angle: " << p_angle << " dE: " << p_de <<
//      " E: " << p_e << std::endl;
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
make_de_e_graph(std::shared_ptr<ROOT::Math::Interpolator>& angle_de,
                std::shared_ptr<ROOT::Math::Interpolator>& angle_e) {
  std::vector<Double_t> des, es;
  Double_t *p_des, *p_es;

  for (Double_t angle = 55; angle <= 70; ++angle) {
    Double_t p_de = angle_de->Eval(angle);
    Double_t p_e = angle_e->Eval(angle);
    des.push_back(p_de);
    es.push_back(p_e);
//    cout << "angle: " << angle << "; " << "de: " << p_de <<
//      "; " << "e: " << p_e << std::endl;
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
  // graph->SetMinimum(-1.2);
  // graph->SetMaximum(1.2);

  return graph;
}


void plot_ref_e_de() {
  auto ref_angle_de_e_pair = compute_ref_de_e();
  auto ref_angle_de = std::get<0>(ref_angle_de_e_pair);
  auto ref_angle_e = std::get<1>(ref_angle_de_e_pair);
  auto graph = make_de_e_graph(ref_angle_de, ref_angle_e);

  graph->Draw("C SAME");
}

void plot_espri_e_de(RootScript& script, ElasticScatteringCuts& cuts) {
  static int hist_num = 1;

  TString prefix = cuts.espri_selector == Espri::left ? "esl" : "esr";
  TString title = TString::Format("dE-E (%s)",
                                  cuts.GetTotalCutTitle().Data());
  TString hist_name = TString::Format("de_e%i", hist_num);
  TString draw_cmd = TString::Format("p_e_eff:p_de");
  TCut user_cut = (TString::Format("p_de > 1 && p_e_eff > 5")).Data();

  //std::cout << "drawing de-e with command: " << draw_cmd << std::endl;
  hdraw(g_chain_total, hist_name, draw_cmd, "(100,0,12,100,0,180)",
        cuts.GetTotalCut() && user_cut, title, "dE (MeV)", "E (MeV)", "colz");
  gPad->SetLogz();

  ++hist_num;
}

void plot_espri_e_de_sim(RootScript& script, ElasticScatteringCuts& cuts) {
  assert(cuts.espri_selector != Espri::both && "Cannot plot for both espri's");

  static int hist_num = 1;
  TString prefix = cuts.espri_selector == Espri::left ? "esl" : "esr";
  TString title = TString::Format("dE-E (%s)",
                                  cuts.GetTotalCutTitle().Data());
  TString hist_name = TString::Format("%s_de_e%i", prefix.Data(), hist_num);
  TString draw_cmd = TString::Format("p_e_sim:p_de_sim");

  hdraw(g_chain_total, hist_name, draw_cmd, "(100,0,12,100,0,180)",
        cuts.GetTotalCut(), title, "dE (MeV)", "E (MeV)", "colz");
  gPad->SetLogz();

  ++hist_num;
}

void plot_espri_e_de_1d(RootScript& script, ElasticScatteringCuts& cuts,
                        Double_t espri_e_cut_width_marker = -1) {
  static int hist_num = 1;
  TString title = TString::Format("dE-E (%s)",
                                  cuts.GetTotalCutTitle().Data());
  TString hist_name = TString::Format("de_e_1d%i", hist_num);
  TString draw_cmd = TString::Format("(p_e_eff - p_e_sim) / p_e_sim");
  TCut user_cut = (TString::Format("p_de > 1 && p_e_eff > 5")).Data();

  TH1* hist = hdraw(g_chain_total, hist_name, draw_cmd, "(100,-1.5,1.5)",
                    cuts.GetTotalCut() && user_cut, title, "(p_e_eff - p_e_sim)/p_e_sim (MeV)",
                    "Counts", "");

  if (espri_e_cut_width_marker > 0) {
    Double_t cut_line_x = espri_e_cut_width_marker / 2;
    gPad->Update();
    draw_marker_line(-cut_line_x, 0, -cut_line_x, gPad->GetUymax());
    draw_marker_line(cut_line_x, 0, cut_line_x, gPad->GetUymax());
  }

  ++hist_num;
}

void he4_add_espri_cuts() {
  init_dataset();
  RootScript script("he4_add_espri_cuts");
  ElasticScatteringCuts cuts;
  //cuts.phi_width = 4;

  cuts.apply_vertexXY_cut = false;
  cuts.apply_espri_aang_cut = false;
  cuts.apply_phi_cut = false;
  cuts.apply_theta_cut = false;

  script.NewPage(2,2);

  script.cd(1);
  cuts.espri_selector = Espri::left;
  plot_espri_e_de(script, cuts);
  plot_ref_e_de();

  script.cd(2);
  cuts.espri_selector = Espri::right;
  plot_espri_e_de(script, cuts);
  plot_ref_e_de();

  cuts.apply_vertexXY_cut = true;
  cuts.apply_espri_aang_cut = true;
  cuts.apply_phi_cut = true;
  cuts.apply_theta_cut = false;

  script.cd(3);
  cuts.espri_selector = Espri::left;
  plot_espri_e_de(script, cuts);
  plot_ref_e_de();

  script.cd(4);
  cuts.espri_selector = Espri::right;
  plot_espri_e_de(script, cuts);
  plot_ref_e_de();

  script.NewPage(2,2);

  script.cd(1);
  plot_espri_e_de_sim(script, cuts);
  plot_ref_e_de();

  script.cd(2);
  cuts.espri_selector = Espri::both;
  plot_espri_e_de(script, cuts);
  plot_ref_e_de();

  script.cd(3);
  cuts.espri_selector = Espri::left;
  plot_espri_e_de_1d(script, cuts, 20);

  script.cd(4);
  cuts.espri_selector = Espri::right;
  plot_espri_e_de_1d(script, cuts, 20);

  script.NewPage(2,2);
  cuts.espri_selector = Espri::both;
  TCut e_cut = "abs(p_e_eff - p_e_sim) < 22";

  script.cd(1);
  draw_theta_corr_2d(g_chain_total, cuts.GetTotalCut(), "No E cut");

  script.cd(2);
  draw_theta_corr_2d(g_chain_total, cuts.GetTotalCut() && e_cut,
                     "|E - E_sim| < 1 sigma");

  script.cd(3);
  draw_theta_corr_1d(g_chain_total, cuts.GetTotalCut(), "No E cut");

  script.cd(4);
  draw_theta_corr_1d(g_chain_total, cuts.GetTotalCut() && e_cut,
                     "|E - E_sim| < 1 sigma");

  script.NewPage(4,3);
  cuts.espri_selector = Espri::both;

  script.cd(1);
  plot_espri_e_de_1d(script, cuts, 0.1);
  script.cd(5);
  e_cut = "abs(p_e_eff - p_e_sim) / p_e_sim < 0.1";
  draw_theta_corr_1d(g_chain_total, cuts.GetTotalCut() && e_cut);
  script.cd(9);
  draw_theta_corr_2d(g_chain_total, cuts.GetTotalCut() && e_cut);

  script.cd(2);
  plot_espri_e_de_1d(script, cuts, 0.2);
  script.cd(6);
  e_cut = "abs(p_e_eff - p_e_sim) / p_e_sim < 0.2";
  draw_theta_corr_1d(g_chain_total, cuts.GetTotalCut() && e_cut);
  script.cd(10);
  draw_theta_corr_2d(g_chain_total, cuts.GetTotalCut() && e_cut);

  script.cd(3);
  plot_espri_e_de_1d(script, cuts, 0.4);
  script.cd(7);
  e_cut = "abs(p_e_eff - p_e_sim) / p_e_sim < 0.4";
  draw_theta_corr_1d(g_chain_total, cuts.GetTotalCut() && e_cut);
  script.cd(11);
  draw_theta_corr_2d(g_chain_total, cuts.GetTotalCut() && e_cut);

  script.cd(4);
  plot_espri_e_de_1d(script, cuts, 0.8);
  script.cd(8);
  e_cut = "abs(p_e_eff - p_e_sim) / p_e_sim < 0.8";
  draw_theta_corr_1d(g_chain_total, cuts.GetTotalCut() && e_cut);
  script.cd(12);
  draw_theta_corr_2d(g_chain_total, cuts.GetTotalCut() && e_cut);
}

