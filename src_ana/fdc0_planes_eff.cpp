/**
 *   \file fdc0_planes_eff.cpp
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


class Fdc0Eff {

public:

  Fdc0Eff(TChain* chain, Int_t eval_plane_no) :
    logger_{"Fdc0Eff"}, chain_{chain}, eval_plane_no_{eval_plane_no} {
      init_hists();
  }

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

  TH1* xplane_eff;
  TH1* xplane_pos;
  TH1* xplane_pos_extrp;
  TH1* yplane_eff;
  TH1* yplane_pos;
  TH1* yplane_pos_extrp;

private:

  TH1* init_xpos_hist(TString title) {
    // return s13::ana::make_th1(32, 0, 160, title, "X [mm]", "Counts");
    return s13::ana::make_th1((120-40)/5, 40, 120, title, "X [mm]", "Counts");
  }

  TH1* init_ypos_hist(TString title) {
    // return s13::ana::make_th1(32, 0, 160, title, "Y [mm]", "Counts");
    return s13::ana::make_th1((120-40)/5, 40, 120, title, "Y [mm]", "Counts");
  }

  void init_hists() {
    for (int i = 0; i < 4; i++) {
      TString title = s13::ana::tstrfmt("FDC0 X%d position", i+1);
      xplanes_pos.push_back(init_xpos_hist(title));
      title = s13::ana::tstrfmt("FDC0 Y%d position", i+1);
      yplanes_pos.push_back(init_ypos_hist(title));
    }

    xplane_pos = init_xpos_hist("FDC0 Xn position");
    yplane_pos = init_ypos_hist("FDC0 Yn position");
    xplane_pos_extrp = init_xpos_hist("FDC0 Xn position extrp.");
    yplane_pos_extrp = init_ypos_hist("FDC0 Yn position extrp.");
  }

  void bind_tree_vars() {
    TString vpref = "FDC0";
    chain_->SetBranchAddress(vpref + "_x1wid", &x1wid);
    chain_->SetBranchAddress(vpref + "_x2wid", &x2wid);
    chain_->SetBranchAddress(vpref + "_x3wid", &x3wid);
    chain_->SetBranchAddress(vpref + "_x4wid", &x4wid);
    chain_->SetBranchAddress(vpref + "_y1wid", &y1wid);
    chain_->SetBranchAddress(vpref + "_y2wid", &y2wid);
    chain_->SetBranchAddress(vpref + "_y3wid", &y3wid);
    chain_->SetBranchAddress(vpref + "_y4wid", &y4wid);
  }

  void clean_tree_vars() {
    x1wid =
      x2wid =
      x3wid =
      x4wid =
      y1wid =
      y2wid =
      y3wid =
      y4wid = -99999;
  }

  void compute_plane_eff_values(double x, double y, double x_extrp, double y_extrp) {
    xplane_pos->Fill(x);
    xplane_pos_extrp->Fill(x_extrp);
    // if (std::abs(x - x_extrp) <= 1.0 * 5) {
    // }

    yplane_pos->Fill(y);
    yplane_pos_extrp->Fill(y_extrp);
    // if (std::abs(y - y_extrp) <= 1.0 * 5) {
    // }
  }

  void process_event() {
    double xx[4];
    double yy[4];
    double xx_extrp[4];
    double yy_extrp[4];
    for (int i = 0; i < 4; i++) {
      xx[i] = -9999; yy[i] = -9999; xx_extrp[i] = -9999; yy_extrp[i] = -9999;
    }

    // logger_.debug("X1wid: %d", x1wid);
    xx[0] = 5 * (x1wid - 0.0);
    xx[1] = 5 * (x2wid - 0.5);
    xx[2] = 5 * (x3wid - 0.0);
    xx[3] = 5 * (x4wid - 0.5);
    yy[0] = 5 * (y1wid - 0.0);
    yy[1] = 5 * (y2wid - 0.5);
    yy[2] = 5 * (y3wid - 0.0);
    yy[3] = 5 * (y4wid - 0.5);

    xplanes_pos[0]->Fill(xx[0]);
    xplanes_pos[1]->Fill(xx[1]);
    xplanes_pos[2]->Fill(xx[2]);
    xplanes_pos[3]->Fill(xx[3]);
    yplanes_pos[0]->Fill(yy[0]);
    yplanes_pos[1]->Fill(yy[1]);
    yplanes_pos[2]->Fill(yy[2]);
    yplanes_pos[3]->Fill(yy[3]);

    xx_extrp[0] = 2 * ((xx[1] - xx[3]) / 2) + xx[2];
    xx_extrp[1] = 2 * ((xx[0] - xx[2]) / 2) + xx[3];
    xx_extrp[2] = 2 * ((xx[3] - xx[1]) / 2) + xx[0];
    xx_extrp[3] = 2 * ((xx[2] - xx[0]) / 2) + xx[1];
    yy_extrp[0] = 2 * ((yy[1] - yy[3]) / 2) + yy[2];
    yy_extrp[1] = 2 * ((yy[0] - yy[2]) / 2) + yy[3];
    yy_extrp[2] = 2 * ((yy[3] - yy[1]) / 2) + yy[0];
    yy_extrp[3] = 2 * ((yy[2] - yy[0]) / 2) + yy[1];
    compute_plane_eff_values(xx[eval_plane_no_ - 1],
                             yy[eval_plane_no_ - 1],
                             xx_extrp[eval_plane_no_ - 1],
                             yy_extrp[eval_plane_no_ - 1]);
  }

  TH1* compute_eff_hist(TH1* num, TH1* denom) {
    auto num_clone = (TH1*)num->Clone();
    num_clone->Divide(denom);
    num_clone->SetTitle(s13::ana::tstrfmt("Efficiency (%s)", num_clone->GetTitle()));
    // num_clone->SetMaximum(1.5);
    num_clone->SetMinimum(0.0);
    return num_clone;
  }

  void finalize() {
    xplanes_effs.push_back(compute_eff_hist(xplanes_pos[0], xplanes_pos[2]));
    xplanes_effs.push_back(compute_eff_hist(xplanes_pos[1], xplanes_pos[3]));
    xplanes_effs.push_back(compute_eff_hist(xplanes_pos[2], xplanes_pos[0]));
    xplanes_effs.push_back(compute_eff_hist(xplanes_pos[3], xplanes_pos[1]));
    yplanes_effs.push_back(compute_eff_hist(yplanes_pos[0], yplanes_pos[2]));
    yplanes_effs.push_back(compute_eff_hist(yplanes_pos[1], yplanes_pos[3]));
    yplanes_effs.push_back(compute_eff_hist(yplanes_pos[2], yplanes_pos[0]));
    yplanes_effs.push_back(compute_eff_hist(yplanes_pos[3], yplanes_pos[1]));
    xplane_eff = compute_eff_hist(xplane_pos, xplane_pos_extrp);
    yplane_eff = compute_eff_hist(yplane_pos, yplane_pos_extrp);
  }

  s13::misc::MessageLogger logger_;

  TChain* chain_;

  Int_t eval_plane_no_;

  Int_t x1wid;
  Int_t x2wid;
  Int_t x3wid;
  Int_t x4wid;
  Int_t y1wid;
  Int_t y2wid;
  Int_t y3wid;
  Int_t y4wid;
};

std::shared_ptr<Fdc0Eff> eval_efficiency(TChain* chain, int eval_plane_no) {
  std::shared_ptr<Fdc0Eff> eff_alg(new Fdc0Eff(chain, eval_plane_no));
  eff_alg->compute();
  return eff_alg;
}

void draw_fdc0_efficiency(s13::ana::RootScript& script,
                         std::shared_ptr<Fdc0Eff> fdc0_eff) {
  script.NewPage(4, 2).cd(s13::ana::PadSeq::column);
  for (int i = 0; i < 4; i++) {
    script.cd();
    fdc0_eff->xplanes_effs[i]->Draw();
    // gPad->SetLogy();
    script.cd();
    fdc0_eff->yplanes_effs[i]->Draw();
    // gPad->SetLogy();
  }

  script.NewPage(4,2).cd(s13::ana::PadSeq::column);
  for (int i = 0; i < 4; i++) {
    script.cd();
    fdc0_eff->xplanes_pos[i]->Draw();
    // gPad->SetLogy();
    script.cd();
    fdc0_eff->yplanes_pos[i]->Draw();
    // gPad->SetLogy();
  }

  script.NewPage(2,2).cd(s13::ana::PadSeq::column);
  gStyle->SetStatY(0.95);
  gStyle->SetStatX(1.0);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.2);

  script.cd();
  fdc0_eff->xplane_pos->Draw();
  gPad->SetLogy();
  script.cd();
  fdc0_eff->xplane_pos_extrp->Draw();
  gPad->SetLogy();
  script.cd();
  fdc0_eff->yplane_pos->Draw();
  gPad->SetLogy();
  script.cd();
  fdc0_eff->yplane_pos_extrp->Draw();
  gPad->SetLogy();

  script.NewPage(2,1).cd(s13::ana::PadSeq::column);
  gStyle->SetStatY(0.50);
  gStyle->SetStatX(1.0);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.15);

  script.cd();
  fdc0_eff->xplane_eff->Draw();
  script.cd();
  fdc0_eff->yplane_eff->Draw();
}

void fdc0_planes_eff() {
  int eval_plane_no = 2;
  TChain* chain = tchain_from_file(cli_opts.input_filename(), "DSDCsTree");
  auto fdc0_eff = eval_efficiency(chain, eval_plane_no);

  s13::ana::RootScript script(cli_opts.output_file_basename("FDC0", "_eff"));
  gStyle->SetOptStat(1111111);
  draw_fdc0_efficiency(script, fdc0_eff);
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

  fdc0_planes_eff();
  std::cout << "exiting main\n" << std::flush;
}
