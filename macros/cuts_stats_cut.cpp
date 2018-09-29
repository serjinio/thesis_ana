#include <iostream>
#include "TChain.h"
#include <string.h>

#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

#include "common.hpp"
#include "init_4he_ds.hpp"
#include "elacuts.hpp"


using namespace s13::ana;


void inline print_pdf_first(TCanvas& canv) {
  TString output_filename = TString::Format("out/%s.pdf(", __FILENAME__);
  canv.Print(output_filename, "pdf");
}

void inline print_pdf(TCanvas& canv) {
  TString output_filename = TString::Format("out/%s.pdf", __FILENAME__);
  canv.Print(output_filename, "pdf");
}

void inline print_pdf_last(TCanvas& canv) {
  TString output_filename = TString::Format("out/%s.pdf)", __FILENAME__);
  canv.Print(output_filename, "pdf");
}


void draw_phi_corr(TChain& chain, float phi_width = 0, TString draw_opts = "colz") {
  static int hist_num = 1;
  ElasticScatteringCuts cuts;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_vertexZ_cut = true;
  cuts.apply_phi_cut = false;
  cuts.apply_theta_cut = false;

  hdraw(chain, cuts.GetTotalCutName() + "_" + std::to_string(hist_num),
        "abs(p_phi*57.3-s1dc_phi*57.3)",
        "(300,150,210)",
        cuts.GetTotalCut(),
        "abs(p_{#phi} - He_{#phi}) (" + cuts.GetTotalCutTitle() + ")",
        "p_{#phi} - He_{#phi} [deg]",
        "Counts");
  draw_marker_line(177.5 - -phi_width/2., 0,
                   177.5 - -phi_width/2., 80);
  draw_marker_line(177.5 - phi_width/2., 0,
                   177.5 - phi_width/2., 80);
  hist_num += 1;
}

void draw_vertz_corr(TChain& chain, Espri espri, float vertz_width = 0,
                     TString draw_opts = "colz") {
  assert(espri != Espri::both);
  static int hist_num = 1;
  ElasticScatteringCuts cuts;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_vertexZ_cut = false;
  cuts.apply_phi_cut = true;
  cuts.apply_theta_cut = false;

  TString espri_det = espri == Espri::left ? "esl_" : "esr_";

  hdraw(g_chain_up, "vertz_" + std::to_string(hist_num),
        espri_det + "vertex_zpos",
        "(200,-150,150)",
        cuts.GetTotalCut() && "",
        "Vertex Z (" + cuts.GetTotalCutTitle() + ")",
        "Z [mm]",
        "Counts");
  Double_t zpos_offset = espri == Espri::left ? \
    ElasticScatteringCuts::k_esl_zpos_offset : \
    ElasticScatteringCuts::k_esr_zpos_offset;
  draw_marker_line(-vertz_width/2. + zpos_offset, 0, 
                   -vertz_width/2. + zpos_offset, 35);
  draw_marker_line(vertz_width/2. + zpos_offset, 0, 
                   vertz_width/2. + zpos_offset, 35);
  hist_num += 1;
}

void draw_tgtxy_corr(TChain& chain, float vertexXY_radius) {
  static int hist_num = 1;
  ElasticScatteringCuts cuts;
  cuts.apply_vertexXY_cut = false;
  cuts.apply_vertexZ_cut = true;
  cuts.apply_phi_cut = true;
  cuts.apply_theta_cut = false;

  hdraw(g_chain_up, "tgtXY_" + std::to_string(hist_num),
        "tgt_up_ypos:tgt_up_xpos",
        "(400,-40,40,400,-40,40)",
        cuts.GetTotalCut(),
        "Vertex X&Y (" + cuts.GetTotalCutTitle() + ")",
        "X [mm]",
        "Y [mm]");
  gPad->SetLogz();
  draw_marker_circle(-2, -2, vertexXY_radius);

  hist_num += 1;
}

void draw_theta_corr(TChain& chain, TCut cut,
                     TString title, int yaxis_limit = -1, TString draw_opts = "colz") {
  static int hist_num = 1;
  TCut p_angle_acceptance_cut = "p_theta*57.3>55 && p_theta*57.3<70";
  TString hist_name = TString::Format("thcorr_%i", hist_num);

  hdraw(chain, hist_name,
        "s1dc_theta*57.3 - he_theta_theor*57.3",
        "(200,-8,8)",
        cut && p_angle_acceptance_cut,
        title,
        "Fragment: #theta_{exp} - #theta_{theor} [lab. deg.]",
        "Counts", draw_opts);

  if (yaxis_limit != -1) {
    TH1* hist = static_cast<TH1*>(gPad->GetPrimitive(hist_name));
    hist->GetYaxis()->SetRangeUser(0, yaxis_limit);
  }
  hist_num += 1;
}

void draw_phi_cut_range(TCanvas& canv, ElasticScatteringCuts& cuts,
                        TChain& chain, int phi_width_initial, int phi_step) {
  canv.Clear(); canv.Divide(4,3);
  int canv_pad_counter = 1;
  for (cuts.phi_width = phi_width_initial;
       cuts.phi_width < phi_width_initial + phi_step * 4;
       cuts.phi_width += phi_step) {
    TString phi_cut_title = TString::Format("4He-p: #phi cut (|#phi|<%.1f)", cuts.phi_width/2);
    TString not_phi_cut_title = TString::Format("4He-p: #phi cut (|#phi|>%.1f)", cuts.phi_width/2);
    canv.cd(canv_pad_counter);
    draw_phi_corr(g_chain_total, cuts.phi_width);
    canv.cd(canv_pad_counter + 4);
    draw_theta_corr(g_chain_total, cuts.GetTotalCut() && cuts.GetPhiCut(),
                    phi_cut_title, 200);
    canv.cd(canv_pad_counter + 8);
    draw_theta_corr(g_chain_total, cuts.GetTotalCut() && !cuts.GetPhiCut(),
                    not_phi_cut_title, 200);
    canv_pad_counter += 1;
  }
}

void draw_vertz_cut_range(TCanvas& canv, ElasticScatteringCuts& cuts,
                          TChain& chain, int vertexZ_width_initial, int vertexZ_step,
                          Espri espri) {
  canv.Clear(); canv.Divide(4,3);
  int canv_pad_counter = 1;
  for (cuts.vertexZ_width = vertexZ_width_initial;
       cuts.vertexZ_width < vertexZ_width_initial + vertexZ_step * 4;
       cuts.vertexZ_width += vertexZ_step) {
    TString vertexZ_cut_title = TString::Format("4He-p: #vertexZ cut (|vertexZ|<%.1f)",
                                                cuts.vertexZ_width/2);
    TString not_vertexZ_cut_title = TString::Format("4He-p: #vertexZ cut (|vertexZ|>%.1f)",
                                                    cuts.vertexZ_width/2);
    canv.cd(canv_pad_counter);
    draw_vertz_corr(g_chain_total, espri, cuts.vertexZ_width);
    canv.cd(canv_pad_counter + 4);
    draw_theta_corr(g_chain_total, cuts.GetTotalCut() && cuts.GetVertexZCut(),
                    vertexZ_cut_title, 80);
    canv.cd(canv_pad_counter + 8);
    draw_theta_corr(g_chain_total, cuts.GetTotalCut() && !cuts.GetVertexZCut(),
                    not_vertexZ_cut_title, 80);
    canv_pad_counter += 1;
  }
}

void draw_tgtxy_cut_range(TCanvas& canv, ElasticScatteringCuts& cuts,
                          TChain& chain, int vertexXY_initial, int vertexXY_step) {
  canv.Clear(); canv.Divide(4,3);
  int canv_pad_counter = 1;
  for (cuts.vertexXY_radius = vertexXY_initial;
       cuts.vertexXY_radius < vertexXY_initial + vertexXY_step * 4;
       cuts.vertexXY_radius += vertexXY_step) {
    TString vertexXY_cut_title = TString::Format("4He-p: vert-XY cut (|vert-XY|<%.1f)",
                                                cuts.vertexXY_radius);
    TString not_vertexXY_cut_title = TString::Format("4He-p: vert-XY cut (|vert-XY|>%.1f)",
                                                    cuts.vertexXY_radius);
    canv.cd(canv_pad_counter);
    draw_tgtxy_corr(g_chain_total, cuts.vertexXY_radius);
    canv.cd(canv_pad_counter + 4);
    draw_theta_corr(g_chain_total, cuts.GetTotalCut() && cuts.GetVertexXYCut(),
                    vertexXY_cut_title, 200);
    canv.cd(canv_pad_counter + 8);
    draw_theta_corr(g_chain_total, cuts.GetTotalCut() && !cuts.GetVertexXYCut(),
                    not_vertexXY_cut_title, 200);
    canv_pad_counter += 1;
  }
}

void cuts_stats_cut() {
  init_dataset();

  TString output_root_filename = TString::Format("out/%s.root", __FILENAME__);
  TFile hist_out(output_root_filename, "RECREATE");
  TCanvas c1("c1");
  gStyle->SetOptStat(1111111);

  ElasticScatteringCuts cuts;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_theta_cut = false;
  cuts.apply_phi_cut = false;
  cuts.apply_vertexZ_cut = true;
  cuts.apply_hodf_he6_cut = false;
  cuts.espri_selector = Espri::both;

  c1.Clear(); c1.Divide(4,1);

  c1.cd(1); draw_phi_corr(g_chain_total);
  c1.cd(2); draw_theta_corr(g_chain_total, cuts.GetTotalCut(),
                            "4He-p: NO #phi cut");
  c1.cd(3); draw_theta_corr(g_chain_total, cuts.GetTotalCut() && cuts.GetPhiCut(),
                            "4He-p: #phi cut");
  c1.cd(4); draw_theta_corr(g_chain_total, cuts.GetTotalCut() && !cuts.GetPhiCut(),
                            "4He-p: NOT #phi cut");
  print_pdf_first(c1);

  // draw_phi_cut_range(c1, cuts, g_chain_total, 4, 8);
  // print_pdf(c1);
  // draw_phi_cut_range(c1, cuts, g_chain_total, 32, 8);
  // print_pdf(c1);

  cuts = ElasticScatteringCuts();
  cuts.apply_vertexXY_cut = true;
  cuts.apply_theta_cut = false;
  cuts.apply_phi_cut = true;
  cuts.apply_vertexZ_cut = false;
  cuts.apply_hodf_he6_cut = false;
  cuts.espri_selector = Espri::left;

  draw_vertz_cut_range(c1, cuts, g_chain_total, 20, 20, Espri::left);
  print_pdf(c1);
  draw_vertz_cut_range(c1, cuts, g_chain_total, 120, 40, Espri::left);
  print_pdf(c1);

  cuts.espri_selector = Espri::right;
  draw_vertz_cut_range(c1, cuts, g_chain_total, 20, 20, Espri::right);
  print_pdf(c1);
  draw_vertz_cut_range(c1, cuts, g_chain_total, 120, 40, Espri::right);
  print_pdf(c1);

  cuts = ElasticScatteringCuts();
  cuts.apply_vertexXY_cut = false;
  cuts.apply_theta_cut = false;
  cuts.apply_phi_cut = true;
  cuts.apply_vertexZ_cut = true;
  cuts.apply_hodf_he6_cut = false;
  cuts.espri_selector = Espri::both;

  draw_tgtxy_cut_range(c1, cuts, g_chain_total, 3, 6);
  print_pdf_last(c1);

  hist_out.Write();
  hist_out.Close();
}
