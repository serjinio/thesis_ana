
#include <fstream>
#include <future>

#include "TChain.h"
#include "TGraph.h"
#include "Math/Interpolator.h"

#include "init_ds.hpp"
#include "macroincludes.hpp"


using namespace s13::ana;


TChain *g_local_chain = init_6he_partial_dataset(50);

// center of the phi_corr: p_theta_eff - He_theta [deg]
// ideally should be 180 [deg]
const Double_t gk_phi_corr_center = 177.5;

// proton center angle, angle width & energy from kinematics
const Double_t gk_center_angle = 67.5 * D2R;
const Double_t gk_pde_center_angle_width = 1.0 * D2R;
const Double_t gk_nai_center_angle_width = 2.0 * D2R;

// proton center angles & energy for pdE calibration
const Double_t gk_pde_calib_num_points = 4;
const Double_t gk_pde_theta_center_angles[] = {60 * D2R, 62.58 * D2R,
                                               65.54 * D2R, 67 * D2R};
const Double_t gk_pde_theta_center_angle_energies[] = {123.36, 103.84,
                                                       83.01, 73.34};
const Double_t gk_pde_theta_center_angle_des[] = {3.12, 3.55,
                                                  4.24, 4.69};

// e of proton from kinematics
const Double_t gk_center_angle_energy = 70.23;
// e loss in plastic scintillator
const Double_t gk_center_de = 3.99;
const Double_t gk_center_e_at_nai = gk_center_angle_energy - gk_center_de;


struct pde_pedestals {
  Double_t right;
  Double_t left;
};

struct nai_pedestals {
  Double_t esl_nai[7];
  Double_t esr_nai[7];
};


TH1* plot_nai_e(RootScript& script, ElasticScatteringCuts& cuts, int nai_idx,
                Double_t center_angle = gk_center_angle,
                Double_t center_angle_width = gk_nai_center_angle_width) {
  assert(cuts.espri_selector != Espri::both && "Cannot plot for both espri's");

  static int hist_num = 1;
  TString prefix = cuts.espri_selector == Espri::left ? "esl" : "esr";
  TString title = TString::Format("NaI %i E (%s)",
                                  nai_idx,
                                  cuts.GetTotalCutTitle().Data());
  TString hist_name = TString::Format("%s_de_e%i", prefix.Data(), hist_num);
  TString draw_cmd = TString::Format("%s_naie[%i]", prefix.Data(), nai_idx);
  Double_t min_angle = center_angle - center_angle_width / 2;
  Double_t max_angle = center_angle + center_angle_width / 2;
  TCut user_cut = (TString::Format("p_theta_eff > %f && p_theta_eff < %f && %s_naie[%i] > 100",
                                   min_angle, max_angle,
                                   prefix.Data(), nai_idx)).Data();
  std::cout << "plot_nai_e(): proton angle cut: " << user_cut << std::endl;

  TH1* hist = hdraw(*g_local_chain, hist_name, draw_cmd, "(100, 0, 4000)",
        cuts.GetTotalCut() && user_cut, title, "E (a.u.)", "Counts", "colz");
  gPad->SetLogz();

  ++hist_num;
  return hist;
}

TH1* plot_pde(TCut cuts, Espri espri_sel,
              Double_t center_angle = gk_center_angle,
              Double_t center_angle_width = gk_pde_center_angle_width,
              TString binning = "(100, 0, 4000)") {
 static int hist_num = 1;

 TString prefix = espri_sel == Espri::left ? "esl" : "esr";
 TCut espri_sel_cut = espri_sel == Espri::left ? "p_phi < 0" : "p_phi > 0";
 TString title = "";
 TCut total_cut = cuts && espri_sel_cut;
 TString hist_name = TString::Format("%s_pde_h%i", prefix.Data(), hist_num);
 TString draw_cmd = TString::Format("%s_de", prefix.Data());

 if (center_angle > 0) {
   Double_t min_angle = center_angle - center_angle_width / 2;
   Double_t max_angle = center_angle + center_angle_width / 2;
   TCut user_cut = TString::
     Format("p_theta_eff*57.3 > %f && p_theta_eff*57.3 < %f",
            min_angle, max_angle).Data();

   std::cout << "plot_pde(): proton angle cut: " << user_cut << std::endl;
   total_cut = total_cut && user_cut;

   title = TString:: Format("pdE #theta_p:(%.1f,%.1f) (%s)",
                            min_angle, max_angle, prefix.Data());
 } else {
   title = TString:: Format("pdE #theta_p:(%.1f,%.1f) (%s)",
                            55., 70., prefix.Data());
 }

 TH1* hist = hdraw(*g_local_chain, hist_name, draw_cmd, binning,
                   total_cut, title, "dE (a.u.)", "Counts", "colz");
 //gPad->SetLogy();
 //hist->GetYaxis()->SetRangeUser(0, 200);
 //hist->Draw();
 //gPad->Update();

 ++hist_num;
 return hist;
}

TH1* plot_pde_ylim(TCut cuts, Espri espri_sel,
                   Double_t ylim,
                   Double_t center_angle = gk_center_angle,
                   Double_t center_angle_width = gk_pde_center_angle_width,
                   TString binning = "(100, 0, 4000)") {
  TH1* hist = plot_pde(cuts, espri_sel, center_angle, center_angle_width, binning);
  hist->GetYaxis()->SetRangeUser(0, ylim);
  // gPad->SetLogy();
  hist->Draw();
  gPad->Update();
  return hist;
}

std::ostream& operator<<(std::ostream& os, const nai_pedestals& obj) {
  os << "ESL NaI:";
  for (int i = 0; i < 7; i++) {
    os << " " << obj.esl_nai[i];
  }
  os << std::endl;

  os << "ESR NaI:";
  for (int i = 0; i < 7; i++) {
    os << " " << obj.esr_nai[i];
  }
  os << std::endl;

  return os;
}

Double_t find_nai_pedestal(Espri espri_sel, int nai_idx) {
  assert(espri_sel != Espri::both);
  static int hist_num = 1;

  auto empty_tgt_run = init_6he_partial_dataset(1, 351);
  TString prefix = espri_sel == Espri::left ? "esl" : "esr";
  TString hist_name = TString::Format("%s_nai_%i", prefix.Data(), hist_num);
  TString title = TString::Format("%s NaI %i pedestal distribution",
                                  espri_sel == Espri::left ? "ESL" : "ESR",
                                  nai_idx);
  TString draw_cmd = TString::Format("%s_naie[%i]", prefix.Data(), nai_idx);
  // just beam trigger and make sure pdE is NOT fired
  TCut cuts = "triggers[1] == 1 && triggers[5] == 0";

  TH1 *hist = hdraw(*empty_tgt_run, hist_name, draw_cmd, "(500, -100, 200)",
                    cuts, title,
                    "E [ch]", "Counts");
  gPad->SetLogz();

  auto fit_params = fit_hist(hist);
  Double_t mean_value = fit_params.params[1];

  hist_num += 1;
  return mean_value;
}

nai_pedestals
find_nai_pedestals(RootScript& script) {
  nai_pedestals peds;

  script.NewPage(4,2).cd(PadSeq::row);
  for (int i = 0; i < 7; i++) {
    script.cd();
    peds.esl_nai[i] = find_nai_pedestal(Espri::left, i);
  }

  script.NewPage(4,2).cd(PadSeq::row);
  for (int i = 0; i < 7; i++) {
    script.cd();
    peds.esr_nai[i] = find_nai_pedestal(Espri::right, i);
  }

  return peds;
}

void compute_nai_cal_consts(RootScript& script, ElasticScatteringCuts& cuts,
                            Double_t *nai_peds, std::ofstream& ofstr) {
  assert(cuts.espri_selector != Espri::both);

  std::vector<Double_t> cal_consts;
  script.NewPage(4,2).cd(PadSeq::row);
  for (int i = 0; i < 7; i++) {
    std::cout << "Computing calibration constant for NaI " << i << "..." <<
      std::endl;
    script.cd();
    ofstr << i << "," << nai_peds[i] << ",";
    TH1 *hist = plot_nai_e(script, cuts, i);
    FitParams fit_params = fit_hist(hist);
    std::cout << "Fit params: " << fit_params << std::endl;
    Double_t cal_const = gk_center_e_at_nai / (fit_params.params[1] - nai_peds[i]);
    cal_consts.push_back(cal_const);
    ofstr << cal_const << std::endl;
    std::cout << "Calibration const: " << cal_const << std::endl;
  }

  int nai_idx = 0;
  std::for_each(cal_consts.begin(), cal_consts.end(),
                [&nai_idx](Double_t &el) {
                  std::cout << "NaI " << nai_idx << ": " << el << std::endl;
                  ++nai_idx;
                });
}

Double_t compute_pde_mean_value(TH1* hist) {
  FitParams fit_params = fit_hist(hist);
  std::cout << "Fit params: " << fit_params << std::endl;
  Double_t mean_value = fit_params.params[1];
  std::cout << "Mean value: " << mean_value << std::endl;
  return mean_value;
}

TGraph* configure_pde_fit_graph(int num_points, Double_t *xes, Double_t *ys,
                                TString title = "PdE fit") {
  TGraph* graph = new TGraph(num_points, xes, ys);
  graph->SetLineColor(kGreen);
  graph->SetLineWidth(1);
  graph->SetMarkerColor(kBlue);
  graph->SetMarkerStyle(21);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitle("Channel [a.u.]");
  graph->GetYaxis()->SetTitle("Proton dE [MeV]");
  graph->SetMinimum(0.0);
  graph->SetMaximum(5.0);

  return graph;
}

FitParams fit_pde_cal_const(TGraph* gr) {
  FitParams fit_params(2);
  gr->Fit("pol1", "W");
  assert(gr->GetFunction("pol1") != nullptr && "fit function cannot be NULL");
  gr->GetFunction("pol1")->Print();
  std::cout << "Fit function param0: " << gr->GetFunction("pol1")->GetParameter(0)
            << std::endl;
  std::cout << "Fit function param1: " << gr->GetFunction("pol1")->GetParameter(1)
            << std::endl;
  Double_t params[] = {gr->GetFunction("pol1")->GetParameter(0),
                       gr->GetFunction("pol1")->GetParameter(1)};
  Double_t errors[] = {gr->GetFunction("pol1")->GetParError(0),
                       gr->GetFunction("pol1")->GetParError(1)};
  fit_params.SetParams(params);
  fit_params.SetErrors(errors);
  return fit_params;
}

FitParams compute_pde_cal_const_multi_point(
    RootScript &script, ElasticScatteringCuts &cuts,
    Double_t pedestal) {
  TString prefix = cuts.espri_selector == Espri::left ? "esl" : "esr";
  TString hist_name = TString::Format("%s_pde_fit", prefix.Data());
  TString hist_title = TString::Format("%s pdE fit", prefix.Data());

  Double_t xes[] = {0,0,0,0};
  Double_t ys[] = {0,0,0,0};
  for (int i = 0; i < gk_pde_calib_num_points; ++i) {
    std::cout << "\tComputing for angle: " <<
              gk_pde_theta_center_angles[i] * 1 / D2R << " deg." << std::endl;
    script.cd(i + 1);
    TH1 *hist = plot_pde(cuts.GetTotalCut(), cuts.espri_selector,
                         gk_pde_theta_center_angles[i] / D2R);
    xes[i] = compute_pde_mean_value(hist);
    ys[i] = gk_pde_theta_center_angle_des[i];
    std::cout << "Setting mean value on fit graph: " << xes[i] <<
              " with energy: " << ys[i] << " MeV" << std::endl;
  }

  TGraph* graph_means = configure_pde_fit_graph(gk_pde_calib_num_points, xes, ys);
  FitParams cal_params = fit_pde_cal_const(graph_means);
  script.cd(gk_pde_calib_num_points + 1);
  graph_means->Draw();
  graph_means->GetXaxis()->SetRangeUser(0.0,1000);
  graph_means->Draw();

  return cal_params;
}

Double_t compute_pde_cal_const(
    RootScript &script, ElasticScatteringCuts &cuts, Double_t pedestal) {
  TString prefix = cuts.espri_selector == Espri::left ? "esl" : "esr";
  TString hist_name = TString::Format("%s_pde_fit", prefix.Data());
  TString hist_title = TString::Format("%s pdE fit", prefix.Data());

  std::cout << "\tComputing pdE cal. const at angle: " <<
            gk_center_angle * 1 / D2R << " deg." << std::endl;
  script.cd();
  TH1 *hist = plot_pde(cuts.GetTotalCut(), cuts.espri_selector,
                       gk_center_angle/D2R);
  Double_t mean_value = compute_pde_mean_value(hist);
  std::cout << "Obtained mean value: " << mean_value <<
            " with dE: " << gk_center_de << " MeV" << std::endl;
  Double_t cal_fact = gk_center_de / (mean_value - pedestal);

  return cal_fact;
}

Double_t find_pde_pedestal(Espri espri_sel) {
  assert(espri_sel != Espri::both);
  static int hist_num = 1;

  auto empty_tgt_run = init_6he_partial_dataset(1, 351);
  TString prefix = espri_sel == Espri::left ? "esl" : "esr";
  TString hist_name = TString::Format("%s_de_%i", prefix.Data(), hist_num);
  TString title = TString::Format("%s pdE pedestal distribution",
                                  espri_sel == Espri::left ? "ESL" : "ESR");
  TString draw_cmd = TString::Format("%s_de", prefix.Data());
  // just beam trigger and make sure pdE is NOT fired
  TCut cuts = "triggers[1] == 1 && triggers[5] == 0";

  TH1 *hist = hdraw(*empty_tgt_run, hist_name, draw_cmd, "(500, 0, 500)",
                    cuts, title,
                    "dE [ch]", "Counts");
  gPad->SetLogz();

  auto fit_params = fit_hist(hist);
  Double_t mean_value = fit_params.params[1];

  hist_num += 1;
  return mean_value;
}

pde_pedestals
find_pde_pedestals(RootScript &script) {
  pde_pedestals peds;

  script.cd();
  peds.left = find_pde_pedestal(Espri::left);
  script.cd();
  peds.right = find_pde_pedestal(Espri::right);

  return peds;
}

void compute_pde_cal_consts(RootScript& script, ElasticScatteringCuts& cuts,
                            std::ofstream& ofstr) {
  script.NewPage(2,1).cd(PadSeq::column);
  auto pde_pedestals = find_pde_pedestals(script);

  script.NewPage(2,1).cd(PadSeq::column);
  std::cout << "Computing calibration constant for pdE left..." << std::endl;
  cuts.espri_selector = Espri::left;
  Double_t left_pde_cal_fact = compute_pde_cal_const(script, cuts, pde_pedestals.left);
  ofstr << "L, " << pde_pedestals.left << "," << left_pde_cal_fact << std::endl;

  std::cout << "Computing calibration constant for pdE right..." << std::endl;
  cuts.espri_selector = Espri::right;
  Double_t right_pde_cal_fact = compute_pde_cal_const(script, cuts, pde_pedestals.right);
  ofstr << "R, " << pde_pedestals.right << "," << right_pde_cal_fact << std::endl;
}

/*
  Draws ESPRIs pdEs dists. with two different triggers given:
  for reaction and beam triggers. To compare for threshold value.
*/
void draw_pdes_btrt(TCut trig1, TCut trig2,
                    Espri espri_sel,
                    double theta_center = 62.5, double theta_width = 15) {
  assert(espri_sel != Espri::both);
  TString prefix = espri_sel == Espri::left ? "ESL" : "ESR";

  TH1* hist_rt = plot_pde(trig1, espri_sel,
                          theta_center, theta_width, "(100, 150, 900)");
  hist_rt->SetTitle(TString::Format("RT - %s pdE dist @%.1f deg.",
                                    prefix.Data(), theta_center));
  hist_rt->SetLineColor(kGreen);

  TH1* hist_bt = plot_pde(trig2, espri_sel,
                          theta_center, theta_width, "(100, 150, 900)");
  hist_bt->SetTitle("ELA");
  hist_bt->SetLineColor(kRed);

  std::cout << "Scaling RT hist by: " << hist_bt->Integral() / hist_rt->Integral()
            << std::endl;
  hist_rt->Scale(hist_bt->Integral() / hist_rt->Integral());

  hist_rt->Draw();
  hist_bt->Draw("SAME");

  gPad->SetTitle();
  gPad->BuildLegend(0.8, 0.4, 0.95, 0.58);
  // gPad->SetLogy();
}

void espri_e_calib() {
  //init_dataset();
  RootScript script("espri_e_calib_pde_dists_rtela_zoom_2_lin");
  ElasticScatteringCuts cuts;
  cuts.vertexXY_radius = 10;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_phi_cut = true;
  cuts.apply_theta_cut = true;
  cuts.apply_hodf_he6_cut = true;
  cuts.apply_espri_aang_cut = true;

  std::ofstream left_espri_consts, right_espri_consts, pde_consts;
  // left_espri_consts.open("data/esl_nai_calib.csv.test");
  // right_espri_consts.open("data/esr_nai_calib.csv.test");
  // pde_consts.open("data/es_pde_calib.csv.test");

  // script.NewPage(2, 2).cd(PadSeq::column);
  // cuts.apply_theta_cut = false;
  // cuts.espri_selector = Espri::left;
  // script.cd();
  // draw_theta_corr_2d(*g_local_chain, cuts.GetTotalCut(), cuts.GetTotalCutTitle());
  // script.cd();
  // draw_theta_corr_1d(*g_local_chain, cuts.GetTotalCut(), cuts.GetTotalCutTitle());
  // cuts.espri_selector = Espri::right;
  // script.cd();
  // draw_theta_corr_2d(*g_local_chain, cuts.GetTotalCut(), cuts.GetTotalCutTitle());
  // script.cd();
  // draw_theta_corr_1d(*g_local_chain, cuts.GetTotalCut(), cuts.GetTotalCutTitle());

  cuts.apply_theta_cut = true;
  cuts.espri_selector = Espri::both;
  script.NewPage(4,2).cd(PadSeq::column);

  TCut bt_cut = ("triggers[1]==1");
  TCut rt_cut = ("triggers[5]==1");
  TCut ela_cut = cuts.GetTotalCut();

  // script.cd();
  // draw_pdes_btrt(rt_cut, ela_cut, Espri::left);
  // script.cd();
  // draw_pdes_btrt(rt_cut, ela_cut, Espri::right);

  script.cd();
  draw_pdes_btrt(rt_cut, ela_cut, Espri::left, 56, 3.75);
  script.cd();
  draw_pdes_btrt(rt_cut, ela_cut, Espri::right, 56, 3.75);

  script.cd();
  draw_pdes_btrt(rt_cut, ela_cut, Espri::left, 60, 3.75);
  script.cd();
  draw_pdes_btrt(rt_cut, ela_cut, Espri::right, 60, 3.75);

  script.cd();
  draw_pdes_btrt(rt_cut, ela_cut, Espri::left, 64, 3.75);
  script.cd();
  draw_pdes_btrt(rt_cut, ela_cut, Espri::right, 64, 3.75);

  script.cd();
  draw_pdes_btrt(rt_cut, ela_cut, Espri::left, 68, 3.75);
  script.cd();
  draw_pdes_btrt(rt_cut, ela_cut, Espri::right, 68, 3.75);

  //compute_pde_cal_consts(script, cuts, pde_consts);

  // auto nai_pedestals = find_nai_pedestals(script);
  // cuts.espri_selector = Espri::left;
  // compute_nai_cal_consts(script, cuts, nai_pedestals.esl_nai, left_espri_consts);
  // cuts.espri_selector = Espri::right;
  // compute_nai_cal_consts(script, cuts, nai_pedestals.esr_nai, right_espri_consts);
}

int main() {
  espri_e_calib();
}
