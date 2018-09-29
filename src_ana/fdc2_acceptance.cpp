
#include "TROOT.h"
#include "TThread.h"
#include "TVector3.h"
#include "TMath.h"

#include "cli.hpp"
#include "init_ds.hpp"
#include "treewalker.hpp"
#include "consts.hpp"
#include "fdc2ypos_heabang.hpp"
#include "drawing.hpp"
#include "cuts_conf.hpp"
#include "rootscript.hpp"


using ElaCuts = s13::ana::ElasticScatteringCuts;
using Evt = s13::ana::ScatteringEvent;
using PtrTreeWalker = std::shared_ptr<s13::ana::ScatteringTreeWalker>;


static s13::misc::CliOptions cli_opts;


std::shared_ptr<s13::ana::Fdc2YposHeAbangAlg>
compute_ypos_relation(ElaCuts cuts) {
  auto alg = s13::ana::make_fdc2_heabang_alg(cuts);
  auto tree_walker = PtrTreeWalker(new s13::ana::ScatteringTreeWalker());
  tree_walker->add(alg);

  auto ds = init_dataset_from_cli(cli_opts);

  if (!cli_opts.use_mt()) {
    tree_walker->walk(*ds[0]);
  } else {
    s13::ana::ParallelTreeWalker ptree_walker(tree_walker, ds);
    ptree_walker.walk();
  }

  return alg;
}

void
save_heabang_fdc2y_fit_params(std::shared_ptr<s13::ana::Fdc2YposHeAbangAlg> fdc2_ypos_alg) {
  fdc2_ypos_alg->fit_ypos_heabang_hists("pol1");
  auto hists = fdc2_ypos_alg->ypos_heabang_hists();
  std::ofstream fit_params_fs("data/heabang_fdc2y_conv_factors.csv");
  fit_params_fs << "# He_aang,Fit_const,Fit_coeff,Fit_const_err,Fit_coeff_err\n";
  int hist_idx = 0;
  for (auto el : hists) {
    s13::ana::FitParams params(2);
    params.SetParams(el->GetFunction("pol1")->GetParameters());
    params.SetErrors(el->GetFunction("pol1")->GetParErrors());
    fit_params_fs << fdc2_ypos_alg->ypos_heabang_hist_aang(hist_idx) << ",";
    fit_params_fs << params;
    hist_idx += 1;
  }
}

void draw_hebang_fdc2y_hists(s13::ana::RootScript& script,
                             std::vector<TH2*> hists) {
  int hists_per_page = hists.size() / 2;
  if (hists_per_page % 2 != 0) {
    hists_per_page += 1;
  }
  size_t hist_counter = 0;
  for (auto el : hists) {
    if (hist_counter == hists.size() / 2 || hist_counter == 0) {
      script.NewPage(hists_per_page / 2, 2).cd(s13::ana::PadSeq::row);
    }

    script.cd();
    el->Draw("col");

    hist_counter += 1;
  }
}

TH2* make_empty_hebang_fdc2y_hist(double he_aang, int nbins = 200) {
  auto hist_title =
    s13::ana::tstrfmt("S1DC_bang vs. FDC2_y (S1DC_aang=%.1f)",
                      he_aang);
  auto hist = s13::ana::make_th2(nbins, -6, 6, nbins, -600, 600,
                                 hist_title,
                                 "S1DC_bang [deg]",
                                 "FDC2 Y [mm]");
  return hist;
}

TH2* make_interp_hebang_fdc2y_hist(double he_aang) {
  s13::ana::HeAbangFdc2YConverter conv("data/heabang_fdc2y_conv_factors.csv");
  auto hist = make_empty_hebang_fdc2y_hist(he_aang, 50);
  hist->SetLineColor(kGreen);
  // hist->SetFillColor(kBlack);
  for (double bang = -6; bang < 6; bang += 0.1) {
    double fdc2_ypos = conv.fdc2_ypos(he_aang, bang);
    hist->Fill(bang, fdc2_ypos, 1000);
  }
  return hist;
}

void draw_interp_hebang_fdc2y_hists(s13::ana::RootScript& script,
                                    double min_he_aang,
                                    double max_he_aang,
                                    double npoints,
                                    std::shared_ptr<s13::ana::Fdc2YposHeAbangAlg> alg) {
  double bin_width = (max_he_aang - min_he_aang) / npoints;
  int nhists_per_page = npoints / 2;
  if (nhists_per_page % 2 != 0) {
    nhists_per_page += 1;
  }
  int hist_counter = 0;
  for (double aang = min_he_aang + bin_width/2;
       aang < max_he_aang; aang += bin_width) {
    if (hist_counter == nhists_per_page - 1 || hist_counter == 0) {
      script.NewPage(nhists_per_page/2, 2).cd(s13::ana::PadSeq::row);
    }
    auto hist = make_interp_hebang_fdc2y_hist(aang);
    script.cd();
    alg->ypos_heabang_hists().at(hist_counter)->Draw("col");
    hist->Draw("same box");
    hist_counter += 1;
  }
}

void draw_fdc2_abang_acceptance_hist(s13::ana::RootScript& script) {
  TH2* hist_abang = s13::ana::make_th2(120, -15, 15, 120, -9, 9,
                                       "FDC2 He a&b angles acceptance",
                                       "He aang [deg.]", "He bang [deg.]");
  TH2* hist_theta_phi = s13::ana::make_th2(60, 0, 15, 1200, -180, 180,
                                           "FDC2 He #theta & #phi angles acceptance",
                                           "He #theta [deg.]", "He #phi [deg.]");
  s13::ana::HeAbangFdc2YConverter conv("data/heabang_fdc2y_conv_factors.csv");

  for (double aang = -9.99; aang < 9.99; aang += 0.25) {
    for (double bang = -8; bang < 8; bang += 0.15) {
      double ypos = conv.fdc2_ypos(aang, bang);
      if (std::abs(ypos) < s13::gk_fdc2_acceptance_height/2) {
        TVector3 he_dir(0, 0, 1);
        he_dir.RotateX(bang * s13::gk_d2r);
        he_dir.RotateY(aang * s13::gk_d2r);

        hist_abang->Fill(aang, bang);

        double phi = (he_dir.Phi() / s13::gk_d2r) - 90;
        if (phi < -180) {
          phi += 360;
        }
        hist_theta_phi->Fill(he_dir.Theta() / s13::gk_d2r, phi);
      }
    }
  }

  script.NewPage(1).cd();
  hist_abang->Draw("box");
  // script.NewPage(1).cd();
  // hist_theta_phi->Draw("box");
}

TVector3 projection_yz(TVector3 vec) {
  vec.SetX(0);
  return vec;
}

TVector3 projection_xz(TVector3 vec) {
  vec.SetY(0);
  return vec;
}

double cm_to_he_lab_angle(double cm_angle) {
  double cm_angle_rad = cm_angle * s13::gk_d2r;
  double lab_he_angle_rad = TMath::Sin(cm_angle_rad);
  lab_he_angle_rad /= TMath::Cos(cm_angle_rad) + 6./1.;
  return lab_he_angle_rad / s13::gk_d2r;
}

void draw_fdc2_lab_acceptance_hist(s13::ana::RootScript& script) {
  TH2* hist_theta_phi = s13::ana::make_th2(60, 0, 15, 90, -180, 180,
                                           "FDC2 He #theta & #phi angles acceptance",
                                           "He #theta [deg.]", "He #phi [deg.]");
  s13::ana::HeAbangFdc2YConverter conv("data/heabang_fdc2y_conv_factors.csv");

  for (double theta = 0; theta < 10; theta += 0.25) {
    for (double phi = -180; phi < 180; phi += 4) {
      TVector3 dir_vec(0, 0, 0);
      dir_vec.SetMagThetaPhi(1, theta * s13::gk_d2r, phi * s13::gk_d2r);
      TVector3 z_dir(0, 0, 1);

      double aang = projection_xz(dir_vec).Angle(z_dir) / s13::gk_d2r;
      double bang = projection_yz(dir_vec).Angle(z_dir) / s13::gk_d2r;
      if (phi < 0) {
        bang = -bang;
      }
      if ((phi > -180 && phi < -90) || (phi > 90 && phi < 180)) {
        aang = -aang;
      }
      double ypos = conv.fdc2_ypos(aang, bang);
      if (std::abs(ypos) < s13::gk_fdc2_acceptance_height/2) {
        double s13_phi = phi - 90;
        if (s13_phi < -180) {
          s13_phi += 360;
        }
        hist_theta_phi->Fill(theta, s13_phi);
      }
    }
  }

  hist_theta_phi->Draw("box");
}

void draw_fdc2_cm_acceptance_hist(s13::ana::RootScript& script) {
  TH2* hist_theta_phi = s13::ana::make_th2(120, 0, 120, 90, -180, 180,
                                           "FDC2 CM #theta & #phi angles acceptance",
                                           "He #theta [deg.]", "He #phi [deg.]");
  s13::ana::HeAbangFdc2YConverter conv("data/heabang_fdc2y_conv_factors.csv");

  for (double theta = 5; theta < 120; theta += 1) {
    for (double phi = -180; phi < 180; phi += 4) {
      // convert to lab angle
      double lab_he_theta = cm_to_he_lab_angle(theta);
      TVector3 dir_vec(0, 0, 0);
      dir_vec.SetMagThetaPhi(1, lab_he_theta * s13::gk_d2r, phi * s13::gk_d2r);
      // then into a & b angles
      TVector3 z_dir(0, 0, 1);
      double aang = projection_xz(dir_vec).Angle(z_dir) / s13::gk_d2r;
      double bang = projection_yz(dir_vec).Angle(z_dir) / s13::gk_d2r;
      if (phi < 0) {
        bang = -bang;
      }
      if ((phi > -180 && phi < -90) || (phi > 90 && phi < 180)) {
        aang = -aang;
      }

      double ypos = conv.fdc2_ypos(aang, bang);
      if (std::abs(ypos) < s13::gk_fdc2_acceptance_height/2) {
        double s13_phi = phi - 90;
        if (s13_phi < -180) {
          s13_phi += 360;
        }
        hist_theta_phi->Fill(theta, s13_phi);
      }
    }
  }

  hist_theta_phi->Draw("box");
}

void fdc2_acceptance() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  // auto fdc2_alg = compute_ypos_relation(cuts);

  s13::ana::RootScript script("fdc2_acceptance");
  gStyle->SetOptStat(1111111);
  // gStyle->SetOptFit();
  // save_heabang_fdc2y_fit_params(fdc2_alg);
  // draw_hebang_fdc2y_hists(script, fdc2_alg->ypos_heabang_hists());
  // draw_interp_hebang_fdc2y_hists(script, -10, 10, 10, fdc2_alg);
  draw_fdc2_abang_acceptance_hist(script);
  draw_fdc2_lab_acceptance_hist(script);
  draw_fdc2_cm_acceptance_hist(script);
}

int main(int argc, char** argv) {
  s13::misc::MessageLogger log("main()");
  gROOT->SetStyle("Modern");
  TThread::Initialize();
  cli_opts.require_cuts_conf_opt(true);
  try {
    cli_opts.parse_options(argc, argv);
  } catch (std::exception& ex) {
    log.error("Invalid argument provided: %s", ex.what());
    return -1;
  }

  fdc2_acceptance();
  std::cout << "exiting main\n" << std::flush;
}
