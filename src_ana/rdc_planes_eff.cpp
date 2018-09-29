/**
 *   \file rdc_planes_eff.cpp
 *   \brief Computes RDC planes efficiency given input rawtree file
 *
 */

#include "TROOT.h"
#include "TPDF.h"
#include "TPaveStats.h"

#include "cli.hpp"
#include "init_ds.hpp"
#include "consts.hpp"
#include "rootscript.hpp"


static s13::misc::CliOptions cli_opts;


class RdcEff {

public:

  RdcEff(TChain* chain, s13::ana::Espri espri_selector) :
    chain_{chain}, espri_selector_{espri_selector} {
      assert(espri_selector_ != s13::ana::Espri::both && "Invalid ESPRI selector: both");
      init_hists();
  };

  void compute() {
    bind_tree_vars();

    for (int i_entry = 0; i_entry < chain_->GetEntries(); i_entry++) {
      clean_tree_vars();
      chain_->GetEntry(i_entry);
      process_event();
    }

    finalize();
  }

  std::vector<TH1*> xplanes_pos;
  std::vector<TH1*> yplanes_pos;
  std::vector<TH1*> xplanes_effs;
  std::vector<TH1*> yplanes_effs;

private:

  TH1* init_pos_hist(TString title, TString axis_name) {
    return s13::ana::make_th1(32, 0, 450, title, axis_name, "Counts");
  }

  TH1* init_xpos_hist(TString title) {
    return init_pos_hist(title, "X [mm]");
  }

  TH1* init_ypos_hist(TString title) {
    return init_pos_hist(title, "Y [mm]");
  }

  void init_hists() {
    for (int i = 0; i < 4; i++) {
      TString espri = espri_selector_ == s13::ana::Espri::left ? "left" : "right";
      TString title = s13::ana::tstrfmt("ESPRI %s X%d position", espri.Data(), i);
      xplanes_pos.push_back(init_xpos_hist(title));
      title = s13::ana::tstrfmt("ESPRI %s Y%d position", espri.Data(), i);
      yplanes_pos.push_back(init_ypos_hist(title));
    }
  }

  void bind_tree_vars() {
    TString vpref = "ESPRI";
    if (espri_selector_ == s13::ana::Espri::left) {
      vpref += "_RDCL";
    } else {
      vpref += "_RDCR";
    }

    chain_->SetBranchAddress(vpref + "_X0wid", &x0wid);
    chain_->SetBranchAddress(vpref + "_X1wid", &x1wid);
    chain_->SetBranchAddress(vpref + "_X2wid", &x2wid);
    chain_->SetBranchAddress(vpref + "_X3wid", &x3wid);
    chain_->SetBranchAddress(vpref + "_Y1wid", &y1wid);
    chain_->SetBranchAddress(vpref + "_Y2wid", &y2wid);
    chain_->SetBranchAddress(vpref + "_Y3wid", &y3wid);
  }

  void clean_tree_vars() {
    x0wid =
      x1wid =
      x2wid =
      x3wid =
      y1wid =
      y2wid =
      y3wid = -99999;
  }

  void process_event() {
    double x0 = 14 * (x0wid - 0.5);
    double x1 = 14 * (x1wid - 0.5);
    double x2 = 14 * (x2wid - 0.0);
    double x3 = 14 * (x3wid - 0.0);
    double y1 = 14 * (y1wid - 0.0);
    double y2 = 14 * (y2wid - 0.0);
    double y3 = 14 * (y3wid - 0.5);

    double x0_extrp = 10 * ((x2 - x3) / 20) + x1;

    // xplanes_pos[1]->Fill(x1);
    xplanes_pos[1]->Fill(x0_extrp);

    xplanes_pos[2]->Fill(x2);
    xplanes_pos[3]->Fill(x3);
    yplanes_pos[2]->Fill(y2);
    yplanes_pos[3]->Fill(y3);

    if (std::abs(x0_extrp - x0) < 1 * 15) {
      xplanes_pos[0]->Fill(x0);
    }

    // if (std::abs(x1wid - x0wid) < 2.5) {
    //   xplanes_pos[0]->Fill(x0);
    // }

    double mean_ywid = (y3wid + y2wid) / 2.;
    if (std::abs(mean_ywid - y1wid) < 1) {
      yplanes_pos[1]->Fill(14 * (y1wid - 0.0));
    }
  }

  TH1* compute_eff_hist(TH1* num, TH1* denom) {
    auto num_clone = (TH1*)num->Clone();
    num_clone->Divide(denom);
    num_clone->SetTitle(s13::ana::tstrfmt("Efficiency (%s)", num_clone->GetTitle()));
    num_clone->SetMaximum(1.6);
    num_clone->SetMinimum(0);
    return num_clone;
  }

  void finalize() {
    xplanes_effs.push_back(compute_eff_hist(xplanes_pos[0], xplanes_pos[1]));
    xplanes_effs.push_back(compute_eff_hist(xplanes_pos[1], xplanes_pos[2]));
    xplanes_effs.push_back(compute_eff_hist(xplanes_pos[2], xplanes_pos[1]));
    xplanes_effs.push_back(compute_eff_hist(xplanes_pos[3], xplanes_pos[2]));
    yplanes_effs.push_back(compute_eff_hist(yplanes_pos[0], yplanes_pos[1]));
    yplanes_effs.push_back(compute_eff_hist(yplanes_pos[1], yplanes_pos[2]));
    yplanes_effs.push_back(compute_eff_hist(yplanes_pos[2], yplanes_pos[1]));
    yplanes_effs.push_back(compute_eff_hist(yplanes_pos[3], yplanes_pos[2]));
  }

  TChain* chain_;
  s13::ana::Espri espri_selector_;

  Int_t x0wid;
  Int_t x1wid;
  Int_t x2wid;
  Int_t x3wid;
  Int_t y1wid;
  Int_t y2wid;
  Int_t y3wid;
};

std::shared_ptr<RdcEff> eval_efficiency(TChain* chain, s13::ana::Espri espri_sel) {
  std::shared_ptr<RdcEff> eff_alg(new RdcEff(chain, espri_sel));
  eff_alg->compute();
  return eff_alg;
}

void draw_rdc_efficiency(s13::ana::RootScript& script,
                         std::shared_ptr<RdcEff> rdc_eff) {
  script.NewPage(2, 4).cd(s13::ana::PadSeq::row);
  gStyle->SetStatY(0.55);
  gStyle->SetStatX(1.0);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.25);
  for (int i = 0; i < 4; i++) {
    script.cd();
    rdc_eff->xplanes_effs[i]->Draw();
    script.cd();
    rdc_eff->yplanes_effs[i]->Draw();
  }

  script.NewPage(2,4).cd(s13::ana::PadSeq::row);
  for (int i = 0; i < 4; i++) {
    script.cd();
    rdc_eff->xplanes_pos[i]->Draw();
    script.cd();
    rdc_eff->yplanes_pos[i]->Draw();
  }
}

void rdc_planes_eff() {
  TChain* chain = tchain_from_file(cli_opts.input_filename(), "ESPRITree");
  auto rdc_eff = eval_efficiency(chain, s13::ana::Espri::left);

  s13::ana::RootScript script(cli_opts.output_file_basename("RDC", "_eff"));
  draw_rdc_efficiency(script, rdc_eff);
}

int main(int argc, char** argv) {
  s13::misc::MessageLogger log("main()");
  // s13::ana::SetRootStyle(s13::ana::RootStyle::publication);

  cli_opts.require_cuts_conf_opt(false);
  try {
    cli_opts.parse_options(argc, argv);
  } catch (std::exception& ex) {
    log.error("Invalid argument provided: %s", ex.what());
    return -1;
  }

  rdc_planes_eff();
  std::cout << "exiting main\n" << std::flush;
}
