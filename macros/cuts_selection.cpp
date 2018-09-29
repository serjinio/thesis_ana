
#include "TChain.h"
#include "TGraph.h"

#include "init_4he_ds.hpp"
#include "elacuts.hpp"
#include "anapow.hpp"


using namespace s13::ana;


TCut vertex_zpos_cut{"es_vertex_zpos > -45 && es_vertex_zpos < 45"};


TCut target_cut{"pow(tgt_up_xpos + 2, 2) + pow(tgt_up_ypos + 2, 2) < 36"};
TCut phi_corr_cut_1d{"abs(abs(p_phi*57.3-s1dc_phi*57.3) - 180) < 2.5"};

TCut theta_cut_1d2{"abs(s1dc_theta*57.3 - he_theta_theor*57.3) < 3"};

TCut proton_1d_cut{"p_theta*57.3>55 && p_theta*57.3<70"};

TCut bw_proton_angles = "p_theta*57.3 > 67 && p_theta*57.3 < 70";
TCut fw_proton_angles = "p_theta*57.3 > 55 && p_theta*57.3 < 58";


void hdraw(TTree& tree, TString name, TString draw_cmd,
           TString binning, TCut cuts = "", TString title = "",
           TString xaxis_title = "", TString yaxis_title = "",
           TString draw_opts = "colz") {
  TString hstr = TString::Format("%s >>%s%s", draw_cmd.Data(), name.Data(),
                                 binning.Data());
  tree.Draw(hstr, cuts, draw_opts);
  TH1 *hist = static_cast<TH1*>(gPad->GetPrimitive(name));
  if (yaxis_title != "") {
    hist->GetYaxis()->SetTitle(yaxis_title);
  }
  if (xaxis_title != "") {
    hist->GetXaxis()->SetTitle(xaxis_title);
  }
  if (title != "") {
    hist->SetTitle(title);
  }
}

void draw_marker_line(Double_t x1, Double_t y1, Double_t x2, Double_t y2) {
  TLine *line = new TLine(x1, y1, x2, y2);
  line->SetLineColor(kRed);
  line->Draw();
}

void draw_marker_circle(Double_t center_x, Double_t center_y, Double_t radius) {
  TEllipse* ellipse = new TEllipse(center_x, center_y, radius, radius);
  ellipse->SetLineColor(kRed);
  ellipse->SetFillStyle(0);
  ellipse->SetFillColor(0);
  ellipse->Draw();
}

void draw_theta_corr(TCanvas& canv, Espri espri,
                     float init_theta_cut_width = 4, float theta_cut_step = 2,
                     bool use_vertXY_cut = true, bool use_phi_cut = false,
                     bool use_vertZ_cut = false) {
  ElasticScatteringCuts cuts{};
  cuts.apply_vertexXY_cut = true;
  cuts.apply_vertexZ_cut = true;
  cuts.apply_phi_cut = true;
  cuts.apply_theta_cut = true;

  cuts.theta_width = init_theta_cut_width;
  cuts.espri_selector = espri;

  canv.Clear();
  canv.Divide(4,3);
  bool first_iter = true;
  for (int i = 1; i < 5; i++) {
    if (!first_iter)
      cuts.theta_width = init_theta_cut_width + (theta_cut_step * (i-2));

    // 2d theta correlation with cut markers
    canv.cd(i);
    cuts.apply_theta_cut = false;
    hdraw(g_chain_up, cuts.GetTotalCutName() + "_" + std::to_string(i),
          "p_theta*57.3:s1dc_theta*57.3 - he_theta_theor*57.3",
          "(200,-10,10,200,50,75)",
          cuts.GetTotalCut(),
          cuts.GetTotalCutTitle(),
          "Fragment #theta angle [lab. deg.]",
          "Proton #theta angle [lab. deg.]", "colz");
    gPad->SetLogz();
    draw_marker_line(-cuts.theta_width/2.,50,-cuts.theta_width/2.,75);
    draw_marker_line(cuts.theta_width/2.,50,cuts.theta_width/2.,75);
    draw_marker_line(-10,67,10.,67);
    draw_marker_line(-10,70,10.,70);
    draw_marker_line(-10,55,10.,55);
    draw_marker_line(-10,58,10.,58);

    //1d theta correlation FW proton angles
    canv.cd(i + 4);
    if (!first_iter)
      cuts.apply_theta_cut = true;
    hdraw(g_chain_up, cuts.GetTotalCutName() + "_fw",
          "s1dc_theta*57.3 - he_theta_theor*57.3",
          "(200,-10,10)",
          cuts.GetTotalCut() && fw_proton_angles,
          cuts.GetTotalCutTitle() + " FW",
          "Fragment #theta angle [lab. deg.]",
          "Counts");

    //1d theta correlation BW proton angles
    canv.cd(i + 8);
    cuts.apply_theta_cut = true;
    hdraw(g_chain_up, cuts.GetTotalCutName() + "_bw",
          "s1dc_theta*57.3 - he_theta_theor*57.3",
          "(200,-10,10)",
          cuts.GetTotalCut() && bw_proton_angles,
          cuts.GetTotalCutTitle() + " BW",
          "Fragment #theta angle [lab. deg.]",
          "Counts");

    first_iter = false;
  }
}


void draw_phi_corr(TCanvas& canv, Espri espri,
                   float init_phi_cut_width = 4, float phi_cut_step = 2,
                   bool use_vertXY_cut = true, bool use_phi_cut = false,
                   bool use_vertZ_cut = false) {
  //assert(espri != Espri::both);

  ElasticScatteringCuts cuts{};
  cuts.apply_vertexXY_cut = true;
  cuts.apply_vertexZ_cut = true;
  cuts.apply_phi_cut = true;
  cuts.apply_theta_cut = false;

  cuts.phi_width = init_phi_cut_width;
  cuts.espri_selector = espri;

  canv.Clear();
  canv.Divide(4,3);
  bool first_iter = true;
  for (int i = 1; i < 5; i++) {
    if (!first_iter)
      cuts.phi_width = init_phi_cut_width + (phi_cut_step * (i-2));

    // 2d theta correlation with cut markers
    canv.cd(i);
    cuts.apply_phi_cut = false;
    hdraw(g_chain_up, cuts.GetTotalCutName() + "_" + std::to_string(i),
          "abs(p_phi*57.3-s1dc_phi*57.3)",
          "(300,150,210)",
          cuts.GetTotalCut(),
          "abs(p_{#phi} - He_{#phi}) (" + cuts.GetTotalCutTitle() + ")",
          "p_{#phi} - He_{#phi} [deg]",
          "Counts");
    canv.Update();
    gPad->Update();
    draw_marker_line(177.5 - -cuts.phi_width/2., gPad->GetFrame()->GetY1(),
                     177.5 - -cuts.phi_width/2., gPad->GetFrame()->GetY2());
    draw_marker_line(177.5 - cuts.phi_width/2., gPad->GetFrame()->GetY1(),
                     177.5 - cuts.phi_width/2., gPad->GetFrame()->GetY2());
    //gPad->SetLogy();

    if (!first_iter)
      cuts.apply_phi_cut = true;
    //1d theta correlation FW proton angles
    canv.cd(i + 4);
    hdraw(g_chain_up, cuts.GetTotalCutName() + "_fw",
          "s1dc_theta*57.3 - he_theta_theor*57.3",
          "(200,-10,10)",
          cuts.GetTotalCut() && fw_proton_angles,
          cuts.GetTotalCutTitle() + " FW",
          "Fragment #theta angle [lab. deg.]",
          "Counts");

    //1d theta correlation BW proton angles
    canv.cd(i + 8);
    hdraw(g_chain_up, cuts.GetTotalCutName() + "_bw",
          "s1dc_theta*57.3 - he_theta_theor*57.3",
          "(200,-10,10)",
          cuts.GetTotalCut() && bw_proton_angles,
          cuts.GetTotalCutTitle() + " BW",
          "Fragment #theta angle [lab. deg.]",
          "Counts");

    first_iter = false;
  }
}


void draw_tgtr_corr(TCanvas& canv, Espri espri,
                   float init_tgt_cut_r = 3, float tgt_cut_step = 3,
                   bool use_vertXY_cut = true, bool use_phi_cut = false,
                   bool use_vertZ_cut = false) {
  ElasticScatteringCuts cuts{};
  cuts.apply_vertexXY_cut = true;
  cuts.apply_vertexZ_cut = true;
  cuts.apply_phi_cut = true;
  cuts.apply_theta_cut = false;

  cuts.vertexXY_radius = init_tgt_cut_r;
  cuts.espri_selector = espri;

  canv.Clear();
  canv.Divide(4,3);
  bool first_iter = true;
  for (int i = 1; i < 5; i++) {
    if (!first_iter)
      cuts.vertexXY_radius = init_tgt_cut_r + (tgt_cut_step * (i-2));

    // 2d theta correlation with cut markers
    canv.cd(i);
    cuts.apply_vertexXY_cut = false;
    hdraw(g_chain_up, cuts.GetTotalCutName() + "_" + std::to_string(i),
          "tgt_up_ypos:tgt_up_xpos",
          "(400,-40,40,400,-40,40)",
          cuts.GetTotalCut(),
          "Vertex X&Y (" + cuts.GetTotalCutTitle() + ")",
          "X [mm]",
          "Y [mm]");
    gPad->SetLogz();
    draw_marker_circle(-2, -2, cuts.vertexXY_radius);
    //draw_marker_circle(-2, -2, cuts.vertexXY_radius - 2);

    if (!first_iter)
      cuts.apply_vertexXY_cut = true;
    // 1d theta correlation FW proton angles
    canv.cd(i + 4);
    hdraw(g_chain_up, cuts.GetTotalCutName() + "_fw",
          "s1dc_theta*57.3 - he_theta_theor*57.3",
          "(200,-10,10)",
          cuts.GetTotalCut() && fw_proton_angles,
          cuts.GetTotalCutTitle() + " FW",
          "Fragment #theta angle [lab. deg.]",
          "Counts");

    // 1d theta correlation BW proton angles
    canv.cd(i + 8);
    hdraw(g_chain_up, cuts.GetTotalCutName() + "_bw",
          "s1dc_theta*57.3 - he_theta_theor*57.3",
          "(200,-10,10)",
          cuts.GetTotalCut() && bw_proton_angles,
          cuts.GetTotalCutTitle() + " BW",
          "Fragment #theta angle [lab. deg.]",
          "Counts");

    first_iter = false;
  }
}


void draw_vertz_corr(TCanvas& canv, Espri espri,
                     float init_vertz_cut_width = 8, float vertz_cut_step = 8,
                     bool use_vertXY_cut = true, bool use_phi_cut = false,
                     bool use_vertZ_cut = false) {
  assert(espri != Espri::both);

  ElasticScatteringCuts cuts{};
  cuts.apply_vertexXY_cut = true;
  cuts.apply_vertexZ_cut = true;
  cuts.apply_phi_cut = true;
  cuts.apply_theta_cut = false;

  cuts.vertexZ_width = init_vertz_cut_width;
  cuts.espri_selector = espri;

  canv.Clear();
  canv.Divide(4,3);
  bool first_iter = true;
  for (int i = 1; i < 5; i++) {
    if (!first_iter)
      cuts.vertexZ_width = init_vertz_cut_width + (vertz_cut_step * (i-2));

    TString espri_det = espri == Espri::left ? "esl_" : "esr_";

    // 2d theta correlation with cut markers
    canv.cd(i);
    cuts.apply_vertexZ_cut = false;
    hdraw(g_chain_up, cuts.GetTotalCutName() + "_" + std::to_string(i),
          espri_det + "vertex_zpos",
          "(200,-150,150)",
          cuts.GetTotalCut() && "",
          "Vertex Z (" + cuts.GetTotalCutTitle() + ")",
          "Z [mm]",
          "Counts");
    canv.Update();
    Double_t zpos_offset = espri == Espri::left ? \
      ElasticScatteringCuts::k_esl_zpos_offset : \
      ElasticScatteringCuts::k_esr_zpos_offset;
    draw_marker_line(-cuts.vertexZ_width/2. + zpos_offset, 0,
                     -cuts.vertexZ_width/2. + zpos_offset, gPad->GetFrame()->GetY2());
    draw_marker_line(cuts.vertexZ_width/2. + zpos_offset, 0,
                     cuts.vertexZ_width/2. + zpos_offset, gPad->GetFrame()->GetY2());
    if (!first_iter)
      cuts.apply_vertexZ_cut = true;
    //1d theta correlation FW proton angles
    canv.cd(i + 4);
    hdraw(g_chain_up, cuts.GetTotalCutName() + "_fw",
          "s1dc_theta*57.3 - he_theta_theor*57.3",
          "(200,-10,10)",
          cuts.GetTotalCut() && fw_proton_angles,
          cuts.GetTotalCutTitle() + " FW",
          "Fragment #theta angle [lab. deg.]",
          "Counts");

    //1d theta correlation BW proton angles
    canv.cd(i + 8);
    hdraw(g_chain_up, cuts.GetTotalCutName() + "_bw",
          "s1dc_theta*57.3 - he_theta_theor*57.3",
          "(200,-10,10)",
          cuts.GetTotalCut() && bw_proton_angles,
          cuts.GetTotalCutTitle() + " BW",
          "Fragment #theta angle [lab. deg.]",
          "Counts");

    first_iter = false;
  }
}


void cuts_selection() {
  init_dataset();

  TFile hist_out("out/cuts_selection.root", "RECREATE");
  TCanvas c1("c1", "asymmetry graph", 200,10,700,500);
  gStyle->SetOptStat(1111111);

  s13::ana::ScatteringAsymetryAlg asymmetry_alg{
    g_chain_up, g_chain_down, "4he", "-(S13: 4He asymmetry)", 5, 55, 70};

  draw_theta_corr(c1, Espri::left, 1, 2);
  c1.Print("out/cuts_selection.pdf(", "pdf");
  draw_theta_corr(c1, Espri::left, 9, 2);
  c1.Print("out/cuts_selection.pdf", "pdf");
  // draw_theta_corr(c1, Espri::right, 1, 2);
  // c1.Print("out/cuts_selection.pdf", "pdf");
  // draw_theta_corr(c1, Espri::right, 9, 2);
  // c1.Print("out/cuts_selection.pdf", "pdf");

  draw_phi_corr(c1, Espri::left, 4, 4);
  c1.Print("out/cuts_selection.pdf", "pdf");
  draw_phi_corr(c1, Espri::left, 20, 4);
  c1.Print("out/cuts_selection.pdf", "pdf");
  // draw_phi_corr(c1, Espri::right, 4, 4);
  // c1.Print("out/cuts_selection.pdf", "pdf");
  // draw_phi_corr(c1, Espri::right, 20, 4);
  // c1.Print("out/cuts_selection.pdf", "pdf");

  draw_tgtr_corr(c1, Espri::left, 4, 8);
  c1.Print("out/cuts_selection.pdf", "pdf");
  // draw_tgtr_corr(c1, Espri::right, 3, 3);
  // c1.Print("out/cuts_selection.pdf", "pdf");

  draw_vertz_corr(c1, Espri::left, 20, 40);
  c1.Print("out/cuts_selection.pdf", "pdf");
  draw_vertz_corr(c1, Espri::left, 180, 40);
  c1.Print("out/cuts_selection.pdf)", "pdf");
  // draw_vertz_corr(c1, Espri::right, 20, 40);
  // c1.Print("out/cuts_selection.pdf", "pdf");
  // draw_vertz_corr(c1, Espri::right, 180, 40);
  // c1.Print("out/cuts_selection.pdf)", "pdf");

  hist_out.Write();
  hist_out.Close();
}
