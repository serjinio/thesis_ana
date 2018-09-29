

#include <functional>

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


class CutAngDep : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  static constexpr double k_min_theta = s13::gk_cs_range_start;
  static constexpr double k_max_theta = s13::gk_cs_range_end;

  using CutProc = bool (ElaCuts::*)(const Evt& evt) const;
  using HistInitProc = std::function<TH1*(double, double, double)>;
  using AngularSnData = std::vector< std::vector<TH1*> >;

  std::vector<TH1*> e_dists;
  AngularSnData sn_in;
  AngularSnData sn_out;
  AngularSnData sn_tot;
  AngularSnData outphi;

private:

  s13::misc::MessageLogger logger_;

  std::vector<ElaCuts> cuts_range_;

  size_t bin_num_;

  CutProc cut_proc_;
  HistInitProc cut_hist_init_proc_;

  std::shared_ptr<s13::ana::ERelMeanInterp> he6_nais_e_rel_interp_ =
    s13::ana::make_he6_nai_erel_mean_interp();

  std::shared_ptr<s13::ana::ERelMeanInterp> he4_nais_e_rel_interp_ =
    s13::ana::make_he4_nai_erel_mean_interp();

public:

  CutAngDep(const std::vector<ElaCuts>& cuts,
             CutProc cut_proc,
             std::function<TH1*(double, double, double)> hist_init_proc,
             size_t bin_num = s13::gk_cs_bin_number) :
    logger_{"CutAngDep"}, cuts_range_{cuts}, bin_num_{bin_num},
    cut_proc_{cut_proc}, cut_hist_init_proc_{hist_init_proc} {
      init_hists();
  }

  CutAngDep(const CutAngDep&) = delete;
  CutAngDep(CutAngDep&&) = delete;
  CutAngDep operator=(const CutAngDep&) = delete;
  CutAngDep operator=(CutAngDep&&) = delete;

  virtual void setup() {
  }

  virtual void process_event(Evt& evt) {
    assert(!cuts_range_.empty());
    fill_ang_dep_hists(evt);
  }

  virtual void finalize() {
    auto bg_norm_fact_interp =
      s13::ana::make_bg_norm_factor_interp2(DsType::he4,
                                            cuts_range_[0].dsdc_selector);
    for (size_t i = 0; i < outphi.size(); i++) {
      auto& v = outphi.at(i);
      for (auto& h : v) {
        double theta = bin_center_theta(i);
        double phiW = cuts_range_[0].GetPhiCutWidth(theta);
        h->Scale(bg_norm_fact_interp.eval(theta, phiW) - 1);
      }
    }
  }

  const std::vector<ElaCuts>& cuts_range() {
    return cuts_range_;
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

private:

  TH1* init_kin_corr_hist(TString title) {
    auto th1 = &s13::ana::make_th1;
    auto hist = th1(80, -8, 8, title, "#theta_{exp} - #theta_{theor}", "Counts");
    return hist;
  }

  TH2* init_2d_kin_corr_hist(TString title) {
    auto th2 = &s13::ana::make_th2;
    auto hist = th2(200, 0, 12, 200, 50, 75, title, "#theta_{He}", "#theta_{p}");
    return hist;
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

  void init_hists() {
    for (size_t i = 0; i < bin_num_; i++) {
      auto e_hist = cut_hist_init_proc_(bin_min_theta(i),
                                        bin_max_theta(i), 1);
      e_dists.push_back(e_hist);

      sn_in.push_back(std::vector<TH1*>());
      sn_out.push_back(std::vector<TH1*>());
      sn_tot.push_back(std::vector<TH1*>());
      outphi.push_back(std::vector<TH1*>());
      for (size_t c = 0; c < cuts_range_.size(); c++) {
        sn_in[i].push_back(init_kin_corr_hist("IN"));
        sn_out[i].push_back(init_kin_corr_hist("OUT"));
        sn_tot[i].push_back(init_kin_corr_hist("TOT"));
        outphi[i].push_back(init_kin_corr_hist("OUT_phi"));
      }
    }
  }

  void fill_ang_dep_hists(Evt& evt) {
    Double_t he_theta = cuts_range()[0].dsdc_selector == Dsdcs::s1dc
      ? evt.s1dc_theta : evt.fdc0_theta;
    double p_theta = evt.p_theta_eff;
    double he_theta_diff = he_theta - evt.he_theta_theor;

    if (!is_well_defined(p_theta)
        || !is_well_defined(he_theta) ||
        !evt.IsEspriEvent()) {
      return;
    }

    int evt_bin_idx =  bin_idx(evt);
    Espri espri = evt.IsLeftEspriEvent() ? Espri::left : Espri::right;

    if (evt_bin_idx == -1) {
      return;
    }

    // double e_diff = evt.p_e_eff - evt.p_e_sim;
    double e_diff_rel = (evt.p_e - evt.p_e_sim) / evt.p_e_sim;
    if (cuts_range()[0].espri_e_cut_use_he6_rel_e_corr) {
      // NaI crystal mean correction
      e_diff_rel -= he6_nais_e_rel_interp_->eval(evt);
      // Overall E-dist mean correction
      e_diff_rel -= cuts_range()[0].GetHe6EDistMeanInterp(espri)(evt.p_theta_eff);
    }
    if (cuts_range()[0].espri_e_cut_use_he4_rel_e_corr) {
      // NaI crystal mean correction
      e_diff_rel -= he4_nais_e_rel_interp_->eval(evt);
      // Overall E-dist mean correction
      e_diff_rel -= cuts_range()[0].GetHe4EDistMeanInterp(espri)(evt.p_theta_eff);
    }

    for (size_t c = 0; c < cuts_range_.size(); c++) {
      auto loc_cut = cuts_range_[c];
      loc_cut.apply_espri_e_cut = false;
      if (loc_cut.IsSignalEvent(evt)) {
        e_dists[evt_bin_idx]->Fill(e_diff_rel);
        sn_tot[evt_bin_idx][c]->Fill(he_theta_diff);
        if ((loc_cut.*cut_proc_)(evt)) {
          sn_in[evt_bin_idx][c]->Fill(he_theta_diff);
        } else {
          sn_out[evt_bin_idx][c]->Fill(he_theta_diff);
        }
      }

      if (loc_cut.IsBackgroundEvent(evt) && (loc_cut.*cut_proc_)(evt)) {
        outphi[evt_bin_idx][c]->Fill(he_theta_diff);
      }
    }
  }

};

void draw_cut_marker_lines(std::shared_ptr<CutAngDep> alg,
                           int bin_idx,
                           int cut_idx) {
  auto& cuts = alg->cuts_range()[cut_idx];
  TH1* hist = alg->e_dists[bin_idx];
  double cut_hw =
    cuts.GetEDistSigmaInterp(cli_opts.dataset_type())(alg->bin_center_theta(bin_idx)) *
    cuts.espri_e_width_stdev;
  double cut_min_lim = -cut_hw;
  double cut_max_lim = cut_hw;

  s13::ana::draw_marker_line(cut_min_lim,
                             hist->GetMinimum(),
                             cut_min_lim,
                             hist->GetMaximum());
  s13::ana::draw_marker_line(cut_max_lim,
                             hist->GetMinimum(),
                             cut_max_lim,
                             hist->GetMaximum());
}

using GausMean = s13::ana::NamedType<double, struct gaus_mean_>;
using GausSigma = s13::ana::NamedType<double, struct gaus_sigma_>;
using EdistFitParams = std::pair<GausMean, GausSigma>;

std::vector<EdistFitParams>
fit_ecut_dists(std::shared_ptr<CutAngDep> alg) {
  std::vector<EdistFitParams> edist_fit_params;
  for (size_t i = 0; i < alg->e_dists.size(); i++) {
    std::cout << "Stats in fit region: "
              << alg->e_dists[i]->Integral()
              << " events." << std::endl;
    auto hist = alg->e_dists[i];
    double hist_int = hist->Integral();
    GausMean mean(0);
    GausSigma sigma(0);
    if (hist_int > 600) {
      hist->Fit("gaus", "", "", -0.3, 0.3);
      double* fit_params = hist->GetFunction("gaus")->GetParameters();
      mean = GausMean(fit_params[1]);
      sigma = GausSigma(fit_params[2]);
    } else {
      auto h = std::unique_ptr<TH1>{
        static_cast<TH1*>(alg->e_dists[i]->Clone())};
      h->SetAxisRange(-0.4, 0.4, "X");
      mean = GausMean(h->GetMean());
      sigma = GausSigma(h->GetStdDev());
    }
    edist_fit_params.push_back({mean, sigma});
  }
  return edist_fit_params;
}

void write_ecut_fit_params(std::string fname,
                           std::shared_ptr<CutAngDep> alg,
                           std::vector<EdistFitParams> fit_params) {
  std::ofstream ofs(fname);
  ofs << "# p_theta,gaus_mean,gaus_sigma" << std::endl;
  for (size_t i = 0; i < fit_params.size(); ++i) {
    double theta = alg->bin_center_theta(i);
    GausMean mean = std::get<0>(fit_params[i]);
    GausSigma sigma = std::get<1>(fit_params[i]);
    ofs << theta << "," << mean.get() << "," << sigma.get() << std::endl;
  }
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

void draw_angdep_sn(s13::ana::RootScript& script,
                    std::shared_ptr<CutAngDep> alg_ecut,
                    std::shared_ptr<s13::ana::CarbonYields> alg_carbon) {
  auto fit_params = fit_ecut_dists(alg_ecut);
  TString espri_str =
    alg_ecut->cuts_range()[0].espri_selector == Espri::left ?
    "esl" : "esr";
  TString fname = s13::ana::tstrfmt("out/ecut_%s_dist_fit_params.csv",
                                    espri_str.Data());
  write_ecut_fit_params(fname.Data(), alg_ecut, fit_params);

  for (size_t i = 0; i < alg_ecut->e_dists.size(); i++) {
    script.NewPage(alg_ecut->cuts_range().size(), 4);
    for (size_t c = 0; c < alg_ecut->cuts_range().size(); c++) {
      script.cd();
      auto hist = alg_ecut->e_dists[i];
      hist->Draw();
      draw_cut_marker_lines(alg_ecut, i, c);
      GausMean mean = std::get<0>(fit_params[i]);
      GausSigma sigma = std::get<1>(fit_params[i]);
      draw_mean_marker_line(mean.get(),
                            hist->GetMinimum(),
                            hist->GetMaximum());
      draw_stdev_marker_lines(3 * sigma.get(), mean.get(),
                              hist->GetMinimum(),
                              hist->GetMaximum());

      double e_cut_hw = alg_ecut->cuts_range()[c].espri_e_width_stdev *
        alg_ecut->cuts_range()[c].GetEDistSigmaInterp(cli_opts.dataset_type())(alg_ecut->bin_center_theta(i));
      auto cut_width_text =
        s13::ana::tstrfmt("E_rel: +/-%.2f", e_cut_hw);
      s13::ana::draw_text(cut_width_text);

      script.cd();
      auto s1 = s13::ana::draw_thstack({alg_ecut->sn_tot[i][c],
            alg_ecut->sn_out[i][c]},
        "TOT & OUT E cut", "NOSTACK");
      s1->SetMaximum(s1->GetMaximum() * 0.8);
      gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);

      if (cli_opts.dataset_type() == DsType::he6) {
        script.cd();
        alg_carbon->yields_in[i][c]->SetTitle("C IN");
        auto s2 = s13::ana::draw_thstack({alg_ecut->sn_in[i][c],
              alg_carbon->yields_in[i][c]},
          "Aft. cuts vs. carbon BG", "NOSTACK");
        s2->SetMaximum(s1->GetMaximum() * 0.8);
        gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);

        script.cd();
        auto ela_yields = static_cast<TH1*>(alg_ecut->sn_in[i][c]->Clone());
        // std::cout << "Signal integral (with BG): " << ela_yields->Integral() << std::endl;
        // std::cout << "BG_tot integral: " << alg_carbon->yields_tot[i][c]->Integral() << std::endl;
        // std::cout << "BG_in integral: " << alg_carbon->yields_in[i][c]->Integral() << std::endl;
        ela_yields->Add(alg_carbon->yields_in[i][c], -1);
        // std::cout << "Signal integral (without BG): " << ela_yields->Integral() << std::endl;
        ela_yields->SetTitle("Yields w/o C BG");
        ela_yields->Draw("HIST");
        s13::ana::draw_kin_corr_hist_sn(alg_ecut->bin_center_theta(i),
                                        ela_yields,
                                        alg_ecut->cuts_range()[c].dsdc_selector,
                                        alg_ecut->cuts_range()[c].theta_width);
      } else { // he4
        script.cd();
        auto s2 = s13::ana::draw_thstack({alg_ecut->sn_in[i][c],
              alg_ecut->outphi[i][c]},
          "Aft. cuts vs. OUT_phi BG", "NOSTACK");
        s2->SetMaximum(s1->GetMaximum() * 0.8);
        gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);

        script.cd();
        auto ela_yields = static_cast<TH1*>(alg_ecut->sn_in[i][c]->Clone());
        ela_yields->SetTitle("Yields w/o OUT_phi BG");
        ela_yields->Add(alg_ecut->outphi[i][c], -1);
        ela_yields->Draw();
        s13::ana::draw_kin_corr_hist_sn(alg_ecut->bin_center_theta(i),
                                        ela_yields,
                                        alg_ecut->cuts_range()[c].dsdc_selector,
                                        alg_ecut->cuts_range()[c].theta_width);
      }
    }
  }
}

std::shared_ptr<s13::ana::CarbonYields>
compute_carbon_yields(const std::vector<ElaCuts>& cuts_range) {
  auto num_runs = cli_opts.dataset_size();
  int carbon_num_runs = s13::gk_carbon_num_of_runs;
  if (num_runs != -1) {
    carbon_num_runs = 0.1971 * num_runs;
    if (carbon_num_runs == 0) carbon_num_runs = 3;
  }

  auto carbon_ds = init_carbon_segmented_dataset(1, carbon_num_runs);
  auto carbon_yields_alg =
    s13::ana::walk_alg2<s13::ana::CarbonYields>(carbon_ds,
                                                cli_opts.use_mt(),
                                                cuts_range,
                                                &ElaCuts::IsEspriECut);
  return carbon_yields_alg;
}

void ecut_ang_dep() {
  auto cuts_range = s13::io::parse_cuts_range_config(std::fstream(cli_opts.cuts_config()));
  auto vz_alg =
    s13::ana::walk_alg2<CutAngDep>(init_dataset_from_cli(cli_opts),
                                   cli_opts.use_mt(),
                                   cuts_range,
                                   &ElaCuts::IsEspriECut,
                                   &s13::ana::init_ecut_bin_hist);

  auto script_fname =
    s13::ana::tstrfmt("%s_%s_ecut_ang_dep", cli_opts.output_file_basename().data(),
                      s13::ana::espri_selector_to_str(cuts_range[0].espri_selector).data());
  s13::ana::RootScript script(script_fname);
  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit();
  if (cli_opts.dataset_type() == DsType::he6) {
    auto carbon_yields_alg = compute_carbon_yields(cuts_range);
    draw_angdep_sn(script, vz_alg, carbon_yields_alg);
  } else {
    draw_angdep_sn(script, vz_alg, nullptr);
  }
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

  ecut_ang_dep();
  std::cout << "exiting main\n" << std::flush;
}
