

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
using Dsdcs = s13::ana::Dsdcs;
using ElaCuts = s13::ana::ElasticScatteringCuts;


void draw_cut_marker_lines(double cut_width, TH1* hist) {
  double phi_dist_center = s13::ana::ElasticScatteringCuts::k_phi_dist_center_deg;
  double xmin = phi_dist_center - cut_width/2;
  double xmax = phi_dist_center + cut_width/2;
  s13::ana::draw_marker_line(xmin, hist->GetMinimum(),
                             xmin, hist->GetMaximum());
  s13::ana::draw_marker_line(xmax, hist->GetMinimum(),
                             xmax, hist->GetMaximum());
}

void draw_angdep_sn(s13::ana::RootScript& script,
                    std::shared_ptr<s13::ana::CutAngDep> alg_phi_cut,
                    std::shared_ptr<s13::ana::CarbonYields> alg_carbon) {
  for (size_t i = 0; i < alg_phi_cut->cut_dists.size(); i++) {
    script.NewPage(alg_phi_cut->cuts_range().size(), 4);
    for (size_t c = 0; c < alg_phi_cut->cuts_range().size(); c++) {
      script.cd();
      auto& cuts = alg_phi_cut->cuts_range()[c];
      double bin_theta = alg_phi_cut->bin_center_theta(i);
      double phiW = cuts.GetPhiCutWidth(bin_theta);
      std::cout << "phiW: " << phiW << std::endl;
      alg_phi_cut->cut_dists[i]->Draw();
      draw_cut_marker_lines(phiW, alg_phi_cut->cut_dists[i]);
      auto cut_width_text =
        s13::ana::tstrfmt("#Delta #phi: +/-%.2f deg.", phiW);
      s13::ana::draw_text(cut_width_text);

      script.cd();
      auto s1 = s13::ana::draw_thstack({alg_phi_cut->sn_tot[i][c],
            alg_phi_cut->sn_out[i][c]},
        "TOT & OUT phi cut", "NOSTACK");
      s1->SetMaximum(alg_phi_cut->sn_tot[i][c]->GetMaximum());
      gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);

      if (cli_opts.dataset_type() == DsType::he6) {
        script.cd();
        alg_carbon->yields_in[i][c]->SetTitle("C IN");
        auto s2 = s13::ana::draw_thstack({alg_phi_cut->sn_in[i][c],
              alg_carbon->yields_in[i][c]},
          "Aft. cuts vs. carbon BG", "NOSTACK");
        s2->SetMaximum(alg_phi_cut->sn_in[i][c]->GetMaximum());
        gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);

        script.cd();
        auto ela_yields = static_cast<TH1*>(alg_phi_cut->sn_in[i][c]->Clone());
        ela_yields->SetTitle("Yields w/o C BG");
        ela_yields->Add(alg_carbon->yields_in[i][c], -1);
        ela_yields->Draw();
        s13::ana::draw_kin_corr_hist_sn(alg_phi_cut->bin_center_theta(i),
                                        ela_yields,
                                        alg_phi_cut->cuts_range()[c].dsdc_selector,
                                        alg_phi_cut->cuts_range()[c].theta_width);
      } else { // he4
        script.cd();
        auto s2 = s13::ana::draw_thstack({alg_phi_cut->sn_in[i][c],
              alg_phi_cut->outphi[i][c]},
          "Aft. cuts vs. OUT_phi BG", "NOSTACK");
        s2->SetMaximum(alg_phi_cut->sn_in[i][c]->GetMaximum());
        gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);

        script.cd();
        auto ela_yields = static_cast<TH1*>(alg_phi_cut->sn_in[i][c]->Clone());
        ela_yields->SetTitle("Yields w/o OUT_phi BG");
        ela_yields->Add(alg_phi_cut->outphi[i][c], -1);
        ela_yields->Draw();
        s13::ana::draw_kin_corr_hist_sn(alg_phi_cut->bin_center_theta(i),
                                        ela_yields,
                                        alg_phi_cut->cuts_range()[c].dsdc_selector,
                                        alg_phi_cut->cuts_range()[c].theta_width);
      }
    }
  }
}

std::shared_ptr<s13::ana::CarbonYields>
compute_carbon_yields(const std::vector<ElaCuts>& cuts_range) {
  int carbon_num_runs = s13::gk_carbon_num_of_runs;
  auto num_runs = cli_opts.dataset_size();
  if (num_runs != -1) {
    carbon_num_runs = 0.1971 * num_runs;
    if (carbon_num_runs == 0) carbon_num_runs = 1;
  }

  auto carbon_ds = init_carbon_segmented_dataset(1, carbon_num_runs);
  auto carbon_yields_alg =
    s13::ana::walk_alg2<s13::ana::CarbonYields>(carbon_ds,
                                                cli_opts.use_mt(),
                                                cuts_range,
                                                &ElaCuts::IsPhiCut);
  return carbon_yields_alg;
}

void phi_ang_dep() {
  auto cuts_range = s13::io::parse_cuts_range_config(std::fstream(cli_opts.cuts_config()));
  auto compute_cut_dist_evt = [](const s13::ana::ScatteringEvent& evt) -> double {
    return std::abs(evt.p_phi - evt.he_phi);
  };
  auto init_cut_hist = [](double min_theta, double max_theta) -> TH1* {
    return s13::ana::init_phi_bin_hist(min_theta, max_theta, 90);
  };
  auto phi_alg =
    s13::ana::walk_alg2<s13::ana::CutAngDep>(init_dataset_from_cli(cli_opts),
                                             cli_opts.use_mt(),
                                             cuts_range,
                                             &ElaCuts::IsPhiCut,
                                             init_cut_hist,
                                             compute_cut_dist_evt);
  auto script_fname =
    s13::ana::tstrfmt("%s_%s_phi_cut_ang_dep", cli_opts.output_file_basename().data(),
                      s13::ana::espri_selector_to_str(cuts_range[0].espri_selector).data());
  s13::ana::RootScript script(script_fname);
  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit();
  if (cli_opts.dataset_type() == DsType::he6) {
    auto carbon_yields_alg = compute_carbon_yields(cuts_range);
    draw_angdep_sn(script, phi_alg, carbon_yields_alg);
  } else {
    draw_angdep_sn(script, phi_alg, nullptr);
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

  phi_ang_dep();
}
