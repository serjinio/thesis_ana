
#include "TChain.h"
#include "TGraph.h"

#include "init_ds.hpp"
#include "macroincludes.hpp"
#include "drawcuts.h"

using namespace s13::ana;
using namespace s13::ana::algs;


constexpr Double_t gk_he4_beam_dsc_factor = 2000;
// R < 10 mm
//constexpr Double_t gk_he4_total_incident = 5.248e6 * gk_he4_beam_dsc_factor;
// R < 6 mm
constexpr Double_t gk_he4_total_incident = 2.855e6 * gk_he4_beam_dsc_factor;


auto g_moss_cs = LoadHe4MossLabCs("data/he4_200mev_moss_cs_ay.csv");
//auto g_moss_cms_cs = LoadHe4MossCmsCs("data/he4_200mev_moss_cs.csv");
//TChain *g_local_chain = init_4he_partial_dataset(1);
TChain* g_local_chain = &g_chain_total;

Double_t cs_fit_fn(Double_t *v, Double_t *par) {
  // std::cout << "Evaluating Moss CS at: " << v[0] <<
  //   " with param value: " << par[0] << " fn value: " <<
  //   g_moss_cs->Eval(v[0] * D2R) << std::endl;
  return g_moss_cs->Eval(v[0] * D2R) * par[0];
}

TH1* compute_he4_cs(ElasticScatteringCuts& cuts, Double_t scaling_factor = 1.0) {
  static int hist_num = 1;
  TString title = TString::Format("He4 CS (%s)",
                                  cuts.GetTotalCutTitle().Data());
  TString hist_name = TString::Format("cs_%i", hist_num);
  TString draw_cmd = TString::Format("p_theta_eff*57.3");
  TString user_cut = "p_theta_eff*57.3 >= 55 && p_theta_eff*57.3 <= 70";
  TH1* hist = hdraw(g_chain_total, hist_name, draw_cmd, "(32,55,70)",
                    cuts.GetTotalCut() && user_cut, title,
                    "Lab angle (deg.)", "Yield", "E1");
  //TH1* hist = hdraw();
  //hist->Scale(scaling_factor);
  //hist->Draw("E2");

  std::vector<int> v1;
  ++hist_num;
  return hist;
}

TH1* make_moss_cs_hist(int bin_num, int start_angle, int stop_angle, Double_t moss_cs_scaling_factor = 1) {
  static int hist_num = 1;
  TH1* hist = new TH1F(TString::Format("moss_cs_%i", hist_num), "",
                       bin_num, start_angle, stop_angle);

  for (int i = 1; i <= bin_num; ++i) {
    Double_t bin_center = hist->GetBinCenter(i);
    Double_t moss_cs_value = g_moss_cs->Eval(bin_center * D2R);
    //std::cout << "bin center: " << bin_center << "; weight: " << moss_cs_value << std::endl;
    for (int j = 0; j < moss_cs_scaling_factor * moss_cs_value; ++j) {
      hist->Fill(bin_center);
    }
  }

  hist_num += 1;
  return hist;
}

TH1* compute_he4_cs_error(TH1* exp_cs_hist, Double_t moss_data_fit_param) {
  static int hist_num = 1;
  TString title = TString::Format("He4 CS error vs. Moss");
  TString hist_name = TString::Format("cse_%i", hist_num);
  auto hist = (TH1*)exp_cs_hist->Clone();
  hist->SetTitle(title);
  hist->SetName(hist_name);
  TH1* moss_cs_hist = make_moss_cs_hist(32, 55, 70, moss_data_fit_param);
  moss_cs_hist->Print();
  hist->Add(moss_cs_hist, -1);
  // now draw hist with difference exp_cs - moss_cs
  hist->Draw("E1");

  ++hist_num;
  return hist;
}

Double_t fit_cs_hist(TH1* cs_hist) {
  //create a function with 1 parameter in the range [0, 1e6]
  TF1 *func = new TF1("cs_fit", cs_fit_fn, 0, 1e4, 1);
  func->SetParNames("Scaling_factor");
  func->SetParameter("Scaling_factor", 100);
  cs_hist->Fit("cs_fit");
  const Double_t *fit_params_raw = cs_hist->GetFunction("cs_fit")->GetParameters();
  const Double_t *fit_errors_raw = cs_hist->GetFunction("cs_fit")->GetParErrors();
  std::cout << "CS fit parameter: " << fit_params_raw[0] << std::endl;
  return fit_params_raw[0];
}

void draw_sn(RootScript& script, ElasticScatteringCuts cuts) {
  cuts.apply_theta_cut = false;
  script.NewPage(2, 1).cd(PadSeq::row);
  gStyle->SetOptStat(1111100);
  script.cd();
  draw_theta_corr_2d(*g_local_chain, cuts.GetTotalCut(), cuts.GetTotalCutTitle());
  script.cd();
  std::cout << "drawing 1d S/N ratio..." << std::endl << std::flush;
  draw_theta_corr_1d(*g_local_chain, cuts.GetTotalCut(), cuts.GetTotalCutTitle());
}

void draw_cs(RootScript& script, ElasticScatteringCuts cuts) {
  cuts.apply_theta_cut = false;
  script.NewPage(2, 1).cd(PadSeq::row);
  gStyle->SetOptStat(0000000);
  script.cd(1);
  auto cs_hist = compute_he4_cs(cuts);
  auto cs_hist_copy = (TH1*)cs_hist->Clone();
  Double_t cs_fit_param = fit_cs_hist(cs_hist);
  script.cd(2);
  compute_he4_cs_error(cs_hist_copy, cs_fit_param);
}

TH1* draw_rdc_hit_pattern(RootScript& script, ElasticScatteringCuts cuts) {
  static int hist_num = 1;
  TString title = TString::Format("RDC hit pattern (%s)",
                                  cuts.GetTotalCutTitle().Data());
  TString hist_name = TString::Format("cs_%i", hist_num);
  TString prefix = cuts.espri_selector == Espri::left ? "esl" : "esr";
  TString draw_cmd = TString::Format("%s_ypos", prefix.Data());
  TH1* hist = hdraw(*g_local_chain, hist_name, draw_cmd, "(200, -250, 250)",
                    cuts.GetTotalCut(), title,
                    "Y pos (mm)", "Counts");

  ++hist_num;
  return hist;
}

void draw_cs_bins_bg_conditions(RootScript& script, CrossSectionRetVal& retval) {
  gStyle->SetOptStat(1111111);
  for (int i = 0; i < retval.lab_cs.bins_num; ++i) {
    if (i % 4 == 0) {
      script.NewPage(4,2).cd(PadSeq::column);
    }

    auto hist_bg_scaled = static_cast<TH1*>(retval.hist_bgs.at(i)->Clone());
    auto bg_norm_fact_interp = MakeBgNormFactorInterpolator();
    Double_t norm_fact = bg_norm_fact_interp->Eval(retval.lab_cs.GetBinArg(i));
    hist_bg_scaled->Scale(norm_fact - 1);

    script.cd();
    Double_t cm_angle = (TMath::Pi() - (2 * retval.lab_cs.GetBinArg(i) * D2R) ) / D2R;
//    retval.hist_signals.at(i)->SetTitle(
//        TString::Format("He4 IN #phi @%.2f deg.", cm_angle));
    retval.hist_signals.at(i)->Draw();
    hist_bg_scaled->SetLineColor(kRed);
    hist_bg_scaled->Draw("SAME");
    script.cd();
    auto hist_signal_wo_bg = static_cast<TH1*>(retval.hist_signals.at(i)->Clone());
    hist_signal_wo_bg->SetName(TString::Format("%s-nobg", hist_signal_wo_bg->GetName()).Data());
//    hist_signal_wo_bg->SetTitle(
//        TString::Format("He4 w/o BG @%.2f deg.", cm_angle));
    hist_signal_wo_bg->Add(hist_bg_scaled, -1);
    hist_signal_wo_bg->SetMinimum(0);
    hist_signal_wo_bg->SetLineColor(kGreen);
    hist_signal_wo_bg->Draw();
  }
}

TGraph* draw_diff_with_moss_cs_cm(ScatteringYield cs_cm) {
  auto moss_cs_cm_interp = LoadHe4MossCmsCs("data/he4_200mev_moss_cs.csv");
  ScatteringYield moss_cs_cm{cs_cm.bins_num, cs_cm.range_start, cs_cm.range_end};
  for (int i = 0; i < moss_cs_cm.bins_num; ++i) {
    moss_cs_cm.SetBinValue(
        cs_cm.GetBinArg(i),
        moss_cs_cm_interp->Eval(cs_cm.GetBinArg(i) * D2R));
    moss_cs_cm.SetBinError(
        cs_cm.GetBinArg(i),
        0);
  }
  std::cout << "[DEBUG:draw_diff_with_moss_cs_cm] moss_cs_cm: "
            << moss_cs_cm << std::endl;

  auto diff_cs = (cs_cm - moss_cs_cm) / moss_cs_cm;
  auto diff_cs_graph = BuildTGraphErrors(diff_cs, "CS (S13 - Moss)",
                                         "CM angle [deg]", "Rel. error x100 [%]");
  diff_cs_graph->Draw("ACP");
  return diff_cs_graph;
}

TGraph* draw_cses_ratio(CrossSectionRetVal cs_retval1, CrossSectionRetVal cs_retval2) {
  auto ratio = cs_retval1.cm_cs / cs_retval2.cm_cs;
  auto name = TString::Format("CS (%s / %s)",
                              cs_retval1.name.Data(),
                              cs_retval2.name.Data());
  auto graph = BuildTGraphErrors(ratio, name,
                                 "CM angle [deg]",
                                 "Rel. error x100 [%]");
  graph->Draw("ACP");
  return graph;
}

CrossSectionRetVal compute_cs(ElasticScatteringCuts cuts) {
  CrossSectionUtils cs_utils(gk_he4_total_incident, cuts, 15);
  auto retval = CrossSection(*g_local_chain, cs_utils);
  std::cout << "CS log: " << retval.lab_cs << std::endl;
  std::cout << "CM CS log: " << retval.cm_cs << std::endl;

  return retval;
}

CrossSectionRetVal compute_cm_cs(ElasticScatteringCuts cuts) {
  CrossSectionUtils cs_utils(gk_he4_total_incident, cuts, 15);
  auto retval = CmCrossSection(*g_local_chain, cs_utils);
  //std::cout << "CS log: " << retval.lab_cs << std::endl;
  std::cout << "CM CS log: " << retval.cm_cs << std::endl;

  return retval;
}

TMultiGraph* draw_abs_cs(RootScript &script, CrossSectionRetVal cs_retval) {
  auto graph = BuildTGraphErrors(cs_retval.cm_cs, "CS p-4He@200MeV (S13)",
                                 "CM angle [deg]", "CS [mb/sr]");
  auto moss_cs_graph = BuildMossCmsCsGraph("data/he4_200mev_moss_cs.csv");

  auto cs_w_bg = cs_retval.yield_signal_total;
  cs_retval.utils.normalize_yields(cs_w_bg);
  auto cm_cs_w_bg = cs_retval.utils.lab_cs_to_cm(cs_w_bg);
  auto graph_w_bg = BuildTGraphErrors(cm_cs_w_bg, "CS p-4He@200MeV (S13 w. BG)",
                                      "CM angle [deg]", "CS [mb/sr]");
  graph_w_bg->SetLineColor(kRed);
  graph_w_bg->SetMarkerColor(kBlueYellow);

  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle(TString::Format("%s;%s;%s", "He4 CS",
                               "CM angle [deg.]", "CS [mb/sr]"));
  mg->Add(graph); mg->Add(moss_cs_graph);

  script.NewPage(1,1).cd(1);
  mg->Draw("ACP");
  gPad->BuildLegend();

  return mg;
}

TMultiGraph* draw_abs_cses(RootScript &script,
                           std::vector<CrossSectionRetVal> cs_retvals,
                           bool draw_moss = false) {
  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle(TString::Format("%s;%s;%s", "He4 CS",
                               "CM angle [deg.]", "CS [mb/sr]"));
  if (draw_moss) {
    auto moss_cs_graph = BuildMossCmsCsGraph("data/he4_200mev_moss_cs.csv");
    mg->Add(moss_cs_graph);
  }

  int idx = 1;
  for (auto& cs_retval : cs_retvals) {
    auto graph = BuildTGraphErrors(cs_retval.cm_cs, cs_retval.name,
                                   "CM angle [deg]", "CS [mb/sr]");
    graph->SetLineColor(idx + 4);
    graph->SetMarkerColor(idx + 6);
    mg->Add(graph);
    ++idx;
  }

  mg->Draw("ACP");
  gPad->BuildLegend();
  return mg;
}

void draw_cut_conditions(RootScript& script, ElasticScatteringCuts cuts) {
  script.NewPage(4,2).cd(PadSeq::column);
  s13::cuts::draw_phi_cut_range(script, cuts, 4);
  // script.NewPage(4,2).cd(PadSeq::column);
  // s13::cuts::draw_rdc_aang_cut_range(script, cuts, 4);
//  script.NewPage(4,2).cd(PadSeq::column);
//  s13::cuts::draw_espri_e_cut_range(script, cuts, 4);
}

void draw_bg_conditions(RootScript& script, ElasticScatteringCuts cuts) {
  script.NewPage(4, 2).cd(PadSeq::column);
  s13::cuts::draw_he4_bg_range(script, cuts, 4, 55, 59);
  script.NewPage(4, 2).cd(PadSeq::column);
  s13::cuts::draw_he4_bg_range(script, cuts, 4, 59, 63);
  script.NewPage(4, 2).cd(PadSeq::column);
  s13::cuts::draw_he4_bg_range(script, cuts, 4, 63, 67);
  script.NewPage(4, 2).cd(PadSeq::column);
  s13::cuts::draw_he4_bg_range(script, cuts, 4, 67, 71);
}

void he4_cs_au() {
  init_4he_dataset();
  RootScript script("he4_cs_phi38_cs_cm", g_local_chain);
  gStyle->SetOptFit();
  ElasticScatteringCuts cuts;
  cuts.espri_selector = Espri::right;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_phi_cut = true;
  cuts.apply_theta_cut = true;
  cuts.apply_espri_aang_cut = false;
  cuts.apply_espri_min_e_de_cut = false;
  cuts.apply_espri_e_cut = false;

  cuts.vertexXY_radius = 6;
  cuts.theta_width = 3.2;
  cuts.phi_width = 37.8;
  cuts.espri_aang_width = 180;

  cuts.apply_user_cut = true;
  cuts.user_cut = ("(p_phi < 0 && abs(esl_ypos) < 132) || "
                   "(p_phi > 0 && abs(esr_ypos) < 132)");

  draw_sn(script, cuts);
  // draw_cut_conditions(script, cuts);
  // draw_bg_conditions(script, cuts);

  std::cout << "Computing CM CS for ESR..." << std::endl;
  cuts.espri_selector = Espri::right;
  auto cm_cs_esr_retval = compute_cm_cs(cuts);
  cm_cs_esr_retval.name = "CM CS p-4He (S13 ESR)";
  std::cout << "Computing CM CS for ESL..." << std::endl;
  cuts.espri_selector = Espri::left;
  auto cm_cs_esl_retval = compute_cm_cs(cuts);
  cm_cs_esl_retval.name = "CM CS p-4He (S13 ESL)";

  std::cout << "Computing CS for ESR..." << std::endl;
  cuts.espri_selector = Espri::right;
  auto cs_esr_retval = compute_cs(cuts);
  cs_esr_retval.name = "CS p-4He (S13 ESR)";
  std::cout << "Computing CS for ESL..." << std::endl;
  cuts.espri_selector = Espri::left;
  auto cs_esl_retval = compute_cs(cuts);
  cs_esl_retval.name = "CS p-4He (S13 ESL)";

  // KLUDGE: temporary fix for RDC efficiency
  // cs_esl_retval.cm_cs = cs_esl_retval.cm_cs.multiply(2.083);
  // cs_esr_retval.cm_cs = cs_esr_retval.cm_cs.multiply(1.370);
  // if (cs_retval.utils.cuts().espri_selector == Espri::left) {
  //   cs_retval.cm_cs = cs_retval.cm_cs.multiply(2.083);
  // } else if (cs_retval.utils.cuts().espri_selector == Espri::right) {
  //   cs_retval.cm_cs = cs_retval.cm_cs.multiply(1.370);
  // }

  // draw_abs_cs(script, cs_retval);
  // script.NewPage(1).cd();
  // draw_diff_with_moss_cs_cm(cs_retval.cm_cs);
  // draw_cs_bins_bg_conditions(script, cs_retval);

  script.NewPage(1,1).cd();
  draw_abs_cses(script, {cs_esl_retval, cs_esr_retval}, true);
  script.NewPage(1,1).cd();
  draw_cses_ratio(cs_esl_retval, cs_esr_retval);
  draw_cs_bins_bg_conditions(script, cs_esr_retval);

  script.NewPage(1,1).cd();
  draw_abs_cses(script, {cm_cs_esl_retval, cm_cs_esr_retval}, true);
  script.NewPage(1,1).cd();
  draw_cses_ratio(cm_cs_esl_retval, cm_cs_esr_retval);
  draw_cs_bins_bg_conditions(script, cm_cs_esr_retval);

  script.NewPage(1,1).cd();
  draw_abs_cses(script,
                {cs_esl_retval, cm_cs_esl_retval,
                    cs_esr_retval, cm_cs_esr_retval}, true);

  //  draw_sn(script, cuts);
  //  script.NewPage(2,1).cd(PadSeq::row);
  //  script.cd();
  //  cuts.espri_selector = Espri::left;
  //  draw_rdc_hit_pattern(script, cuts);
  //  script.cd();
  //  cuts.espri_selector = Espri::right;
  //  draw_rdc_hit_pattern(script, cuts);
}
