
#include "TChain.h"
#include "TGraph.h"

#include "ay_init_ds.hpp"
#include "anapow.hpp"


using namespace s13::ana;


TCut vertex_zpos_cut{"es_vertex_zpos > -45 && es_vertex_zpos < 45"};


TCut target_cut{"pow(tgt_up_xpos + 2, 2) + pow(tgt_up_ypos + 2, 2) < 36"};
TCut theta_corr_cut_1d_left_env = "89.22 - p_theta*57.3 - 4 < 2.65*s1dc_theta*57.3";
TCut theta_corr_cut_1d_right_env = "89.22 - p_theta*57.3 + 4 > 3.12*s1dc_theta*57.3";
TCut theta_cut_1d = theta_corr_cut_1d_left_env && theta_corr_cut_1d_right_env;


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

void ay_4he_phi_cut() {
  init_dataset();

  TFile hist_out("out/ay_4he_phi_cut.root", "RECREATE");
  TCanvas c1("c1", "asymmetry graphs", 200,10,700,500);
  gStyle->SetOptStat(1111111);

  s13::ana::ScatteringAsymetryAlg asymmetry_alg{
    g_chain_up, g_chain_down, "4he", "-(S13: 4He asymmetry)", 5, 55, 70};

  c1.Clear();
  c1.Divide(3,3);
  for (int w = 5; w < 50; w += 5) {
    c1.cd((w+1) / 5);
    TCut phi_corr_cut_1d{
      TString::Format("abs(abs(p_phi*57.3-s1dc_phi*57.3) - 180) < %i", w).Data()};
    TString title = TString::Format("tgt(R<6 mm), "
                                    "#phi(+/-%.2f deg), #theta(+/-4 deg)", w/2.);
    hdraw(g_chain_up, TString::Format("phi_cut_w_%i", w), "p_theta*57.3:s1dc_theta*57.3",
          "(200,0.1,20,200,50,75)",
          "triggers[5]==1" && target_cut &&
          phi_corr_cut_1d && theta_cut_1d,
          title,
          "Fragment #theta angle [lab. deg.]",
          "Proton #theta angle [lab. deg.]", "colz");
    gPad->SetLogz();
  }
  c1.Print("out/ay_4he_phi_cut.pdf(", "pdf");

  c1.Clear();
  c1.Divide(3,3);
  for (int w = 5; w < 50; w += 5) {
    std::cout << "Computing asymmetry with phi cut width " << w << \
      " deg..." << std::endl;
    asymmetry_alg.SetThetaCutWidth(8);
    asymmetry_alg.SetPhiCutWidth(w);
    asymmetry_alg.SetTargetCutRadius(6);
    asymmetry_alg.Run();

    TString title = TString::Format("tgt(R<6 mm), "
                                    "#phi(+/-%.2f deg), #theta(+/-4 deg)", w/2.);
    std::vector<TGraph*> graphs = {Build4HeAyGraph(),
                                   BuildAyGraph(asymmetry_alg.GetResult(), title)};
    auto multi_graph = PrepareAyMultiGraph(title, graphs);

    c1.cd((w+1) / 5);
    multi_graph->Draw("ACP");
    //c1.BuildLegend(0.15, 0.67, 0.5, 0.88);
    //c1.BuildLegend();
  }
  c1.Print("out/ay_4he_phi_cut.pdf)", "pdf");

  hist_out.Write();
  hist_out.Close();
}
