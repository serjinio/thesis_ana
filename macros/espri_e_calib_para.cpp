//
// Created by serj on 16/10/25.
//


#include <fstream>
#include <future>

#include "TChain.h"
#include "TGraph.h"
#include "Math/Interpolator.h"

#include "init_6he_ds.hpp"
#include "macroincludes.hpp"


using namespace s13::ana;


// center of the phi_corr: p_theta_eff - He_theta [deg]
// ideally should be 180 [deg]
const Double_t gk_phi_corr_center = 177.5;

// proton center angle, angle width & energy from kinematics
const Double_t gk_center_angle = 67.5 * D2R;
const Double_t gk_pde_center_angle_width = 2.0 * D2R;
const Double_t gk_pde_min_angle = (gk_center_angle -
                                 (gk_pde_center_angle_width / 2));
const Double_t gk_pde_max_angle = (gk_center_angle +
                                 (gk_pde_center_angle_width / 2));

// proton center angles & energy for pdE calibration
const int gk_pde_calib_num_points = 4;
const Double_t gk_pde_theta_center_angles[] = {60 * D2R, 62.58 * D2R,
                                               65.54 * D2R, 67 * D2R};
const Double_t gk_pde_theta_center_angle_energies[] = {123.36, 103.84,
                                                       83.01, 73.34};
const Double_t gk_pde_theta_center_angle_des[] = {3.12, 3.55,
                                                  4.24, 4.69};

const Double_t gk_center_angle_energy = 70.0;
// e loss in plastic scintillator
const Double_t gk_center_de = 4.0;
const Double_t gk_center_e_at_nai = gk_center_angle_energy - gk_center_de;


// hist set
TH1 *esl_pde_calib_hists[gk_pde_calib_num_points];
TH1 *esr_pde_calib_hists[gk_pde_calib_num_points];
TGraph *esl_pde_fit_hist;
TGraph *esr_pde_fit_hist;

void init_hists() {
  for (int i = 0; i < gk_pde_calib_num_points; ++i) {
    esl_pde_calib_hists[i] = new TH1F(TString::Format("esl_pde_%f", gk_pde_theta_center_angles[i]).Data(),
                                     TString::Format("ESL pdE %f deg", gk_pde_theta_center_angles[i] * 1/D2R),
                                     100,0,4000);
    esr_pde_calib_hists[i] = new TH1F(TString::Format("esr_pde_%f", gk_pde_theta_center_angles[i]).Data(),
                                      TString::Format("ESR pdE %f deg", gk_pde_theta_center_angles[i] * 1/D2R),
                                      100,0,4000);
  }
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

FitParams compute_pde_cal_const(ElasticScatteringCuts& cuts) {
  TH1 **hists_set = cuts.espri_selector == Espri::left ? esl_pde_calib_hists : esr_pde_calib_hists;
  TGraph **hist_fit = cuts.espri_selector == Espri::left ? &esl_pde_fit_hist : &esr_pde_fit_hist;
  TString prefix = cuts.espri_selector == Espri::left ? "esl" : "esr";
  TString hist_name = TString::Format("%s_pde_fit", prefix.Data());
  TString hist_title = TString::Format("%s pdE fit", prefix.Data());

  Double_t xes[] = {0,0,0,0};
  Double_t ys[] = {0,0,0,0};
  for (int i = 0; i < gk_pde_calib_num_points; ++i) {
    std::cout << "\tComputing for angle: " <<
              gk_pde_theta_center_angles[i] * 1 / D2R << " deg." << std::endl;
    TH1 *hist = hists_set[i];
    xes[i] = compute_pde_mean_value(hist);
    ys[i] = gk_pde_theta_center_angle_des[i];
    std::cout << "Setting mean value on fit graph: " << xes[i] <<
              " with energy: " << ys[i] << " MeV" << std::endl;
  }

  *hist_fit = configure_pde_fit_graph(gk_pde_calib_num_points, xes, ys);
  FitParams cal_params = fit_pde_cal_const(*hist_fit);

  return cal_params;
}

void compute_pde_cal_consts(ElasticScatteringCuts& cuts,
                            std::ofstream& ofstr) {
  FitParams cal_params;

  std::cout << "Computing calibration constant for pdE left..." << std::endl;
  cuts.espri_selector = Espri::left;
  cal_params = compute_pde_cal_const(cuts);
  assert(cal_params.params.size() == 2);
  ofstr << "L, " << cal_params.params.at(0) << "," << cal_params.params.at(1) << ","
        << cal_params.errors.at(0) << "," << cal_params.errors.at(1) << std::endl;

  std::cout << "Computing calibration constant for pdE right..." << std::endl;
  cuts.espri_selector = Espri::right;
  cal_params = compute_pde_cal_const(cuts);
  ofstr << "R, " << cal_params.params.at(0) << "," << cal_params.params.at(1) << ","
        << cal_params.errors.at(0) << "," << cal_params.errors.at(1) << std::endl;
}

void fill_hists(ElasticScatteringCuts& cuts, ScatteringEvent& evt) {
  for (int i = 0; i < gk_pde_calib_num_points; ++i) {
    Double_t center_angle_deg = gk_pde_theta_center_angles[i] * 1/D2R;
    if (std::abs(evt.p_theta_eff - center_angle_deg) < gk_pde_center_angle_width/2 * 1/D2R) {
      if (cuts.IsLeftEspriEvent(evt)) {
        esl_pde_calib_hists[i]->Fill(evt.p_de_raw);
      } else {
        esr_pde_calib_hists[i]->Fill(evt.p_de_raw);
      }
    }
  }
}

void compute_results(ElasticScatteringCuts& cuts) {
  std::ofstream left_espri_consts, right_espri_consts, pde_consts;
  left_espri_consts.open("data/esl_nai_calib.csv.test");
  right_espri_consts.open("data/esr_nai_calib.csv.test");
  pde_consts.open("data/es_pde_calib.csv.test");
  compute_pde_cal_consts(cuts, pde_consts);
}

void draw_results(RootScript& script) {
  script.NewPage(3,2);
  for (int i = 0; i < gk_pde_calib_num_points; ++i) {
    script.cd(i+1);
    esl_pde_calib_hists[i]->Draw();
  }
  script.cd(gk_pde_calib_num_points+1);
  esl_pde_fit_hist->Draw();
  script.NewPage(3,2);
  for (int i = 0; i < gk_pde_calib_num_points; ++i) {
    script.cd(i+1);
    esr_pde_calib_hists[i]->Draw();
  }
  script.cd(gk_pde_calib_num_points+1);
  esr_pde_fit_hist->Draw();
}

void walk_tree(ElasticScatteringCuts& cuts) {
  init_hists();

  std::cout << "Filling hists with data..." << std::endl;
  ScatteringTreeWalker tree_walker;
  tree_walker.Walk(g_chain_total, [&] (ScatteringEvent& evt) {
    if (!cuts.IsElasticEvent(evt)) {
      return;
    }
    fill_hists(cuts, evt);
  });

  std::cout << "Computing result..." << std::endl;
  compute_results(cuts);
}

void espri_e_calib_para() {
  init_dataset();
  RootScript script("espri_e_calib_para");
  ElasticScatteringCuts cuts;
  cuts.phi_width = 4; // gk_phi_cut_width;

  cuts.apply_vertexXY_cut = true;
  cuts.apply_phi_cut = true;
  cuts.apply_theta_cut = true;
  cuts.apply_hodf_he6_cut = true;

  walk_tree(cuts);
  draw_results(script);
}

int main() {
  espri_e_calib_para();
}