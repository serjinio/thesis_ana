/**
 *   \file s1dc_planes_eff.cpp
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


class S1dcEff {

public:

  S1dcEff(TChain* chain, Int_t eval_plane_no) :
    logger_{"S1dcEff"}, chain_{chain}, eval_plane_no_{eval_plane_no} {
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
    return s13::ana::make_th1((432-72)/12, 72, 432, title, "X [mm]", "Counts");
  }

  TH1* init_ypos_hist(TString title) {
    // return s13::ana::make_th1(32, 0, 160, title, "Y [mm]", "Counts");
    return s13::ana::make_th1((216-36)/12, 36, 216, title, "Y [mm]", "Counts");
  }

  void init_hists() {
    for (int i = 0; i < 4; i++) {
      TString title = s13::ana::tstrfmt("S1DC X%d position", i+1);
      xplanes_pos.push_back(init_xpos_hist(title));
      title = s13::ana::tstrfmt("S1DC Y%d position", i+1);
      yplanes_pos.push_back(init_ypos_hist(title));
    }

    xplane_pos = init_xpos_hist("S1DC Xn position");
    yplane_pos = init_ypos_hist("S1DC Yn position");
    xplane_pos_extrp = init_xpos_hist("S1DC Xn position extrp.");
    yplane_pos_extrp = init_ypos_hist("S1DC Yn position extrp.");
  }

  void bind_tree_vars() {
    TString vpref = "S1MWDC";
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

  void fill_eval_plane_hists(double x, double y, double x_extrp, double y_extrp) {
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
    int eval_plane_idx = eval_plane_no_ - 1;

    for (int i = 0; i < 4; i++) {
      xx[i] = -9999; yy[i] = -9999; xx_extrp[i] = -9999; yy_extrp[i] = -9999;
    }

    xx[0] = 12 * ((x1wid - 0.0) - 2);
    xx[1] = 12 * ((x2wid - 0.5) - 2);
    xx[2] = 12 * ((x3wid - 0.5) + 2);
    xx[3] = 12 * ((x4wid - 0.0) + 2);
    yy[0] = 12 * (y1wid - 0.0);
    yy[1] = 12 * (y2wid - 0.5);
    yy[2] = 12 * (y3wid - 0.0);
    yy[3] = 12 * (y4wid - 0.5);

    xplanes_pos[0]->Fill(xx[0]);
    xplanes_pos[1]->Fill(xx[1]);
    xplanes_pos[2]->Fill(xx[2]);
    xplanes_pos[3]->Fill(xx[3]);
    yplanes_pos[0]->Fill(yy[0]);
    yplanes_pos[1]->Fill(yy[1]);
    yplanes_pos[2]->Fill(yy[2]);
    yplanes_pos[3]->Fill(yy[3]);

    xx_extrp[0] = ((69. + 79.) * ((xx[1] - xx[2]) / (69. + 79.))) + xx[3];
    xx_extrp[1] = ((69. + 79.) * ((xx[0] - xx[3]) / (69. + 79.))) + xx[2];
    xx_extrp[2] = ((69. + 79.) * ((xx[3] - xx[0]) / (69. + 79.))) + xx[1];
    xx_extrp[3] = ((69. + 79.) * ((xx[2] - xx[1]) / (69. + 79.))) + xx[0];
    yy_extrp[0] = ((59. * 2) * ((yy[1] - yy[3]) / (49. * 2))) + yy[2];
    yy_extrp[1] = ((49. * 2) * ((yy[0] - yy[2]) / (59. * 2))) + yy[3];
    yy_extrp[2] = ((59. * 2) * ((yy[3] - yy[1]) / (49. * 2))) + yy[0];
    yy_extrp[3] = ((49. * 2) * ((yy[2] - yy[0]) / (59. * 2))) + yy[1];
    fill_eval_plane_hists(xx[eval_plane_idx],
                          yy[eval_plane_idx],
                          xx_extrp[eval_plane_idx],
                          yy_extrp[eval_plane_idx]);
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

std::shared_ptr<S1dcEff> eval_efficiency(TChain* chain, int eval_plane_no) {
  std::shared_ptr<S1dcEff> eff_alg(new S1dcEff(chain, eval_plane_no));
  eff_alg->compute();
  return eff_alg;
}

void draw_s1dc_efficiency(s13::ana::RootScript& script,
                         std::shared_ptr<S1dcEff> s1dc_eff) {
  // script.NewPage(4, 2).cd(s13::ana::PadSeq::column);
  // for (int i = 0; i < 4; i++) {
  //   script.cd();
  //   s1dc_eff->xplanes_effs[i]->Draw();
  //   // gPad->SetLogy();
  //   script.cd();
  //   s1dc_eff->yplanes_effs[i]->Draw();
  //   // gPad->SetLogy();
  // }

  script.NewPage(4,2).cd(s13::ana::PadSeq::column);
  for (int i = 0; i < 4; i++) {
    script.cd();
    s1dc_eff->xplanes_pos[i]->Draw();
    // gPad->SetLogy();
    script.cd();
    s1dc_eff->yplanes_pos[i]->Draw();
    // gPad->SetLogy();
  }

  script.NewPage(2,2).cd(s13::ana::PadSeq::column);
  gStyle->SetStatY(0.95);
  gStyle->SetStatX(1.0);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.2);

  script.cd();
  s1dc_eff->xplane_pos->Draw();
  gPad->SetLogy();
  script.cd();
  s1dc_eff->xplane_pos_extrp->Draw();
  gPad->SetLogy();
  script.cd();
  s1dc_eff->yplane_pos->Draw();
  gPad->SetLogy();
  script.cd();
  s1dc_eff->yplane_pos_extrp->Draw();
  gPad->SetLogy();

  script.NewPage(1,2).cd(s13::ana::PadSeq::column);
  gStyle->SetStatY(0.50);
  gStyle->SetStatX(0.95);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.20);

  script.cd();
  s1dc_eff->xplane_eff->Draw();
  script.cd();
  s1dc_eff->yplane_eff->Draw();
}

void s1dc_planes_eff() {
  int eval_plane_no = 4;
  TChain* chain = tchain_from_file(cli_opts.input_filename(), "DSDCsTree");
  auto s1dc_eff = eval_efficiency(chain, eval_plane_no);

  s13::ana::RootScript script(cli_opts.output_file_basename("S1DC", "_eff"));
  gStyle->SetOptStat(1111111);
  draw_s1dc_efficiency(script, s1dc_eff);
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

  s1dc_planes_eff();
  std::cout << "exiting main\n" << std::flush;
}
