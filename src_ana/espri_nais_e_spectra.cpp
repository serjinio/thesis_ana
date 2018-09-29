
#include "TROOT.h"
#include "TPDF.h"
#include "TPaveStats.h"

#include "common.hpp"
#include "cli.hpp"
#include "init_ds.hpp"
#include "consts.hpp"
#include "csalg.hpp"
#include "drawing.hpp"
#include "scattyield.hpp"
#include "protonede.hpp"
#include "treewalk.hpp"
#include "cuts_conf.hpp"
#include "rootscript.hpp"
#include "vertzalg.hpp"


static s13::misc::CliOptions cli_opts;


using Script = s13::ana::RootScript;
using Yield = s13::ana::ScatteringYield;
using DsType = s13::ana::DatasetType;
using PolDir = s13::ana::PolarizationDirection;
using Dsdcs = s13::ana::Dsdcs;
using ElaCuts = s13::ana::ElasticScatteringCuts;
using Evt = s13::ana::ScatteringEvent;
using FitParams = s13::ana::FitParams;
using Espri = s13::ana::Espri;
using LinearInterp = s13::ana::LinearInterp;


class NaisESpectra : public s13::ana::SimpleTTreeAlgorithmBase {

private:

  using AngDepHistsVec = std::vector< std::vector<TH1*> >;

  size_t bin_num_ = 7;

  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;

  std::shared_ptr<s13::ana::ERelMeanInterp> he6_e_rel_interp_ =
    s13::ana::make_he6_nai_erel_mean_interp();

  std::shared_ptr<s13::ana::ERelMeanInterp> he4_e_rel_interp_ =
    s13::ana::make_he4_nai_erel_mean_interp();

public:

  static constexpr double k_min_theta = s13::gk_cs_range_start;
  static constexpr double k_max_theta = s13::gk_cs_range_end;

  std::vector<TH2*> nai_adc_hists;
  std::vector<TH2*> nai_e_hists;
  std::vector<TH1*> nai_rel_e_hists;

  AngDepHistsVec ang_dep_nai_rel_e_hists;

public:

  NaisESpectra(const ElaCuts& cuts, size_t bin_num = 7) :
    bin_num_{bin_num}, logger_{"NaisESpectra"}, cuts_{cuts} {
      assert(cuts_.espri_selector != Espri::both
             && "Can work only on L/R ESPRI separately!");
      init_hists();
    }

  NaisESpectra(const NaisESpectra& other) = delete;
  NaisESpectra(NaisESpectra&&) = delete;
  NaisESpectra operator=(const NaisESpectra&) = delete;
  NaisESpectra operator=(NaisESpectra&&) = delete;

  virtual void setup() {
  }

  virtual void process_event(Evt& evt) {
    if (cuts_.IsSignalEvent(evt)) {
      fill_nai_hists(evt);
    }
  }

  virtual void finalize() {
  }

  ElaCuts cuts() const {
    return cuts_;
  }

  size_t bin_num() {
    return bin_num_;
  }

  double bin_width() {
    return (k_max_theta - k_min_theta) / bin_num_;
  }

  double bin_min_theta(int idx) {
    double minth = k_min_theta + bin_width() * idx;
    return minth;
  }

  double bin_max_theta(int idx) {
    double maxth = k_min_theta + bin_width() * (idx + 1);
    return maxth;
  }

  double bin_center_theta(int idx) {
    return (bin_min_theta(idx) + bin_max_theta(idx)) / 2;
  }

  int bin_idx(Evt& evt) {
    if (evt.p_theta_eff - k_min_theta < 0) {
      return -1;
    }

    size_t idx = (evt.p_theta_eff - k_min_theta) / bin_width();
    if (idx < bin_num_) {
      return idx;
    } else {
      return -1;
    }
  }

private:

  TH1* init_kin_corr_hist(TString title) {
    auto th1 = &s13::ana::make_th1;
    auto hist = th1(100, -8, 8, title, "#theta_{exp} - #theta_{theor}", "Counts");
    return hist;
  }

  TH2* init_2d_kin_corr_hist(TString title) {
    auto th2 = &s13::ana::make_th2;
    auto hist = th2(200, 0, 12, 200, 50, 75, title, "#theta_{He}", "#theta_{p}");
    return hist;
  }

  TH1* init_nai_adc_hist(TString user_title) {
    auto th1 = &s13::ana::make_th1;
    auto title = s13::ana::tstrfmt("NaI ADC data (%s)", user_title.Data());
    auto hist = th1(70, 0, 2000, title, "ADC value [ch]", "Counts");
    return hist;
  }

  TH2* init_nai_adc_2dhist(TString user_title) {
    auto th2 = &s13::ana::make_th2;
    auto title = s13::ana::tstrfmt("Proton #theta vs. NaI ADC (%s)",
                                   user_title.Data());
    auto hist = th2(70, 55, 70, 70, 0, 2000, title,
                    "Proton #theta [lab deg]",
                    "ADC value [ch]");
    return hist;
  }

  TH2* init_nai_e_2dhist(TString user_title) {
    auto th2 = &s13::ana::make_th2;
    auto title = s13::ana::tstrfmt("Proton #theta vs. NaI ADC (%s)",
                                   user_title.Data());
    auto hist = th2(70, 55, 70, 70, 0, 180, title,
                    "Proton #theta [lab deg]",
                    "Proton E [MeV]");
    return hist;
  }

  TH1* init_ecut_bin_hist(double bin_min_theta, double bin_max_theta,
                          TString user_title) {
    auto th1 = &s13::ana::make_th1;
    auto tit = s13::ana::tstrfmt("Rel. E spectrum (%.2f<p_{#theta}<%.2f). %s.",
                                 bin_min_theta, bin_max_theta, user_title.Data());
    auto hist = th1(70, -1, 1, tit, "Rel. E [%/100]", "Counts");
    return hist;
  }

  void init_hists() {
    for (int i = 0; i < s13::gk_espri_num_nais; i++) {
      TString title = s13::ana::tstrfmt("NaI #%i", i);
      nai_adc_hists.push_back(init_nai_adc_2dhist(title));
      nai_e_hists.push_back(init_nai_e_2dhist(title));
      auto hist = init_ecut_bin_hist(k_min_theta, k_max_theta, title);
      nai_rel_e_hists.push_back(hist);
    }

    for (size_t b = 0; b < bin_num_; b++) {
      ang_dep_nai_rel_e_hists.push_back(std::vector<TH1*>());
      for (int i = 0; i < s13::gk_espri_num_nais; i++) {
        TString title = s13::ana::tstrfmt("NaI #%i", i);
        auto hist = init_ecut_bin_hist(bin_min_theta(b),
                                       bin_max_theta(b),
                                       title);
        ang_dep_nai_rel_e_hists[b].push_back(hist);
      }
    }
  }

  int get_nai_idx(Evt& evt) {
    int max_nai_idx = -1;
    double max_nai_adc = 0;
    double* naie_raw = evt.IsLeftEspriEvent() ?
      evt.esl_naie_raw : evt.esr_naie_raw;
    for (int i = 0; i < s13::gk_espri_num_nais; i++) {
      if (naie_raw[i] > max_nai_adc) {
        max_nai_adc = naie_raw[i];
        max_nai_idx = i;
      }
    }
    return max_nai_idx;
  }

  virtual void fill_nai_hists(Evt& evt) {
    if (!is_well_defined(evt.p_theta)
        || !is_well_defined(evt.he_theta)) {
      return;
    }
    double xpos = cuts_.IsLeftEspriEvent(evt) ? evt.esl_xpos : evt.esr_xpos;
    if (std::abs(xpos) > 9000) {
      return;
    }

    double p_theta = evt.p_theta_eff;
    // double e_diff_rel = (evt.p_e - evt.p_e_sim) / evt.p_e_sim;
    Espri espri = evt.IsLeftEspriEvent() ? Espri::left : Espri::right;
    double* naie_raw = cuts_.IsLeftEspriEvent(evt) ? evt.esl_naie_raw : evt.esr_naie_raw;
    double* naie = cuts_.IsLeftEspriEvent(evt) ? evt.esl_naie_cal : evt.esr_naie_cal;
    double nai_rel_e[s13::gk_espri_num_nais];

    for (int i = 0; i < s13::gk_espri_num_nais; i++) {
      nai_rel_e[i] = (naie[i] - evt.p_e_sim) / evt.p_e_sim;
      if (cuts_.espri_e_cut_use_he6_rel_e_corr) {
        nai_rel_e[i] -= he6_e_rel_interp_->eval(p_theta, espri, i);
      }
      if (cuts_.espri_e_cut_use_he4_rel_e_corr) {
        nai_rel_e[i] -= he4_e_rel_interp_->eval(p_theta, espri, i);
      }
    }

    int nai_idx = get_nai_idx(evt);
    if (nai_idx == -1) {
      return;
    }

    nai_adc_hists[nai_idx]->Fill(p_theta, naie_raw[nai_idx]);
    nai_e_hists[nai_idx]->Fill(p_theta, naie[nai_idx]);
    nai_rel_e_hists[nai_idx]->Fill(nai_rel_e[nai_idx]);

    for (size_t b = 0; b < bin_num_; b++) {
      if (p_theta < bin_min_theta(b) || p_theta > bin_max_theta(b)) {
        continue;
      }
      ang_dep_nai_rel_e_hists[b][nai_idx]->Fill(nai_rel_e[nai_idx]);
    }
  }
};


void draw_cut_marker_lines(double cut_width, TH1* hist) {
  s13::ana::draw_marker_line_ysym(cut_width / 2,
                                  hist->GetMinimum(),
                                  hist->GetMaximum());
}

void fit_nai_rel_e_dists(std::shared_ptr<NaisESpectra> alg) {
  for (int i = 0; i < s13::gk_espri_num_nais; i++) {
    if (alg->nai_rel_e_hists[i]->Integral() > 50) {
      alg->nai_rel_e_hists[i]->Fit("gaus", "", "");
    }
  }
}

using NaisFitParams = std::vector< std::vector< std::pair<double, double> > >;

NaisFitParams fit_ang_dep_nai_dists(std::shared_ptr<NaisESpectra> alg) {
  static const char* nai_fit_nan = "nai_fit_nan";
  double notanum = std::nan(nai_fit_nan);
  NaisFitParams fit_params;
  for (size_t b = 0; b < alg->bin_num(); b++) {
    fit_params.push_back({});
    for (int i = 0; i < s13::gk_espri_num_nais; i++) {
      double hist_int = alg->ang_dep_nai_rel_e_hists[b][i]->Integral();
      double mean, sigma;
      if (hist_int > 110) {
        alg->ang_dep_nai_rel_e_hists[b][i]->Fit("gaus", "", "", -0.3, 0.3);
        double* fit_params =
          alg->ang_dep_nai_rel_e_hists[b][i]->GetFunction("gaus")->GetParameters();
        mean = fit_params[1];
        sigma = fit_params[2];
        if (std::abs(mean) > 0.5) {
          mean = 0;
          sigma = 0;
        }
      } else if (hist_int > 10) {
        auto h = std::unique_ptr<TH1>{
          static_cast<TH1*>(alg->ang_dep_nai_rel_e_hists[b][i]->Clone())};
        h->SetAxisRange(-0.4, 0.4, "X");
        mean = h->GetMean();
        sigma = h->GetStdDev();
        // auto nai_hist = alg->ang_dep_nai_rel_e_hists[b][i];
      } else {
        mean = 0;
        sigma = 0;
      }
      fit_params[b].push_back(std::pair<double, double>(mean, sigma));
    }
  }
  return fit_params;
}

template <typename Fn>
void write_nai_fit_params(NaisFitParams& fit_params, TString fname,
                          Fn bin_theta_proc) {
  assert(!fit_params.empty() && fit_params[0].size() == s13::gk_espri_num_nais);
  std::ofstream ofs(fname.Data());
  ofs << "#p_theta,";
  for (int nai_idx = 0; nai_idx < s13::gk_espri_num_nais; nai_idx++) {
    auto str = s13::ana::tstrfmt("nai%i_mean,nai%i_sigma", nai_idx + 1, nai_idx + 1);
    ofs << str;
    if (nai_idx < s13::gk_espri_num_nais - 1) {
      ofs << ",";
    }
  }
  ofs << std::endl;

  for (size_t b = 0; b < fit_params.size(); b++) {
    auto& nais_fp = fit_params[b];
    ofs << bin_theta_proc(b) << ",";
    for (size_t nai_idx = 0; nai_idx < nais_fp.size(); nai_idx++) {
      double mean = std::get<0>(nais_fp[nai_idx]);
      double sigma = std::get<1>(nais_fp[nai_idx]);
      ofs << mean << "," << sigma;
      if (nai_idx < s13::gk_espri_num_nais - 1) {
        ofs << ",";
      }
    }
    ofs << std::endl;
  }
}

TMultiGraph*
make_multi_graph(std::vector<TGraph*> graphs, TString title,
                 TString xaxis_label, TString yaxis_label) {
  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle(TString::Format("%s;%s;%s", title.Data(),
                               xaxis_label.Data(), yaxis_label.Data()));
  for (TGraph* gr : graphs) {
    mg->Add(gr);
  }

  return mg;
}

TGraph*
make_means_tgraph(std::vector<double> p_thetas, std::vector<double> nai_means,
                  TString title) {
  static int line_color = 1;
  static int marker_color = 1;
  static int marker_style = 1;
  TGraph* graph = new TGraph(nai_means.size(), &p_thetas[0], &nai_means[0]);
  graph->SetLineColor(line_color++);
  graph->SetLineWidth(1);
  graph->SetLineStyle(3);
  graph->SetMarkerColor(marker_color++);
  graph->SetMarkerStyle(marker_style++);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitle("Proton #theta [deg]");
  graph->GetYaxis()->SetTitle("dE/E [x100%]");
  graph->SetMinimum(50);
  graph->SetMaximum(75);

  return graph;
}

template <typename Fn>
TMultiGraph*
compute_nai_means_lin_fits(NaisFitParams& fit_params, TString fname,
                            Fn bin_theta_proc) {
  assert(!fit_params.empty() && fit_params[0].size() == s13::gk_espri_num_nais);

  std::vector<TGraph*> mean_graphs;
  for (size_t nai_idx = 0; nai_idx < fit_params[0].size(); nai_idx++) {
    std::vector<double> nai_means;
    std::vector<double> p_thetas;
    for (size_t b = 0; b < fit_params.size(); b++) {
      auto& nais_fp = fit_params[b];
      double mean = std::get<0>(nais_fp[nai_idx]);
      nai_means.push_back(mean);
      p_thetas.push_back(bin_theta_proc(b));
    }
    auto title = s13::ana::tstrfmt("NaI#%i", nai_idx);
    mean_graphs.push_back(make_means_tgraph(p_thetas, nai_means, title));
  }

  std::ofstream ofs(fname);
  ofs << "#nai_id,c(const),k(slope)" << std::endl;
  int fit_line_cr = 1;
  for (size_t i = 0; i < mean_graphs.size(); i++) {
    mean_graphs[i]->Fit("pol1");
    TF1* fit_fn = mean_graphs[i]->GetFunction("pol1");
    fit_fn->SetLineColor(fit_line_cr++);
    assert(fit_fn != nullptr);
    double* fit_params = fit_fn->GetParameters();
    ofs << i << "," << fit_params[0] << "," << fit_params[1] << std::endl;
  }

  auto mgraph = make_multi_graph(mean_graphs, "NaIs E-rel means vs. proton theta",
                                 "Proton #theta [deg]", "dE/E [x100%]");
  return mgraph;
}

void draw_mean_marker_line(Double_t mean, Double_t y1, Double_t y2) {
  TLine *line = new TLine(mean, y1, mean, y2);
  line->SetLineColor(kGreen);
  line->SetLineWidth(3);
  line->Draw();
}

void draw_stdev_marker_lines(Double_t stdev, Double_t mean, Double_t y1, Double_t y2) {
  TLine *line = new TLine(mean - stdev, y1, mean - stdev, y2);
  line->SetLineColor(kRed);
  line->SetLineWidth(3);
  TLine *line2 = new TLine(mean + stdev, y1, mean + stdev, y2);
  line2->SetLineColor(kRed);
  line2->SetLineWidth(3);
  line->Draw();
  line2->Draw();
}

void
draw_nais_hists(s13::ana::RootScript& script, std::shared_ptr<NaisESpectra> alg) {
  fit_nai_rel_e_dists(alg);
  auto bin_theta_proc = [alg](int bin_idx) -> double {
    return alg->bin_center_theta(bin_idx);
  };
  auto fit_params = fit_ang_dep_nai_dists(alg);
  write_nai_fit_params(fit_params,
                       script.GetFilename("_nai_fit_params.csv"),
                       bin_theta_proc);

  script.NewPage(4, 2);
  for (int i = 0; i < s13::gk_espri_num_nais; i++) {
    script.cd();
    alg->nai_adc_hists[i]->Draw("colz");
  }

  script.NewPage(4, 2);
  auto pe_graph =
    s13::ana::make_proton_angle_energy_graph(cli_opts.dataset_type(),
                                             alg->cuts().espri_selector);
  pe_graph->SetLineWidth(2);
  for (int i = 0; i < s13::gk_espri_num_nais; i++) {
    script.cd();
    alg->nai_e_hists[i]->Draw("colz");
    pe_graph->Draw("SAME");
  }

  script.NewPage(4, 2);
  for (int i = 0; i < s13::gk_espri_num_nais; i++) {
    script.cd();
    alg->nai_rel_e_hists[i]->Draw();
  }

  for (size_t b = 0; b < alg->bin_num(); b++) {
    script.NewPage(4, 2);
    for (int i = 0; i < s13::gk_espri_num_nais; i++) {
      script.cd();
      auto hist = alg->ang_dep_nai_rel_e_hists[b][i];
      hist->Draw();
      double mean = std::get<0>(fit_params[b][i]);
      double sigma = std::get<1>(fit_params[b][i]);
      draw_mean_marker_line(mean,
                            hist->GetMinimum(),
                            hist->GetMaximum());
      draw_stdev_marker_lines(3 * sigma, mean,
                              hist->GetMinimum(),
                              hist->GetMaximum());
    }
  }

  auto mgraph =
    compute_nai_means_lin_fits(fit_params,
                               script.GetFilename("_nai_e_rel_lin_fit_params.csv"),
                               bin_theta_proc);
  script.NewPage(1,1).cd();
  mgraph->Draw("ALP");
  gPad->BuildLegend(0.7, 0.15, 0.9, 0.40, "");
}

void nais_e_spectra() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  cuts.use_espri_aang_corr = true;
  auto alg =
    s13::ana::walk_alg2<NaisESpectra>(init_dataset_from_cli(cli_opts),
                                      cli_opts.use_mt(),
                                      cuts, 7);
  auto script_fname =
    s13::ana::tstrfmt("%s_%s", cli_opts.output_file_basename().data(),
                      s13::ana::espri_selector_to_str(cuts.espri_selector).data());
  s13::ana::RootScript script(script_fname);
  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit();
  draw_nais_hists(script, alg);
}

int main(int argc, char** argv) {
  s13::misc::MessageLogger log("main()");
  // s13::ana::SetRootStyle(s13::ana::RootStyle::publication);

  TThread::Initialize();
  cli_opts.require_cuts_conf_opt(true);
  try {
    cli_opts.parse_options(argc, argv);
  } catch (std::exception& ex) {
    log.error("Invalid argument provided: %s", ex.what());
    return -1;
  }

  nais_e_spectra();
}
