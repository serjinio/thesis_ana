
#include "TChain.h"
#include "TGraph.h"

#include "ay_init_ds.hpp"
#include "anapow.hpp"


using namespace s13::ana;


TCut vertex_zpos_cut{"es_vertex_zpos > -45 && es_vertex_zpos < 45"};
static TCut esl_vertex_zpos_cut{"esl_vertex_zpos > -7-13 && esl_vertex_zpos < -7+13"};
static TCut esr_vertex_zpos_cut{"esr_vertex_zpos > -9-18 && esr_vertex_zpos < -9+18"};


TCut phi_corr_cut_1d{"abs(abs(p_phi*57.3-s1dc_phi*57.3) - 180) < 2.5"};
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

void ay_4he_r_tgt_cut() {
  init_dataset();

  TFile hist_out("out/ay_4he_r_tgt_cut.root", "RECREATE");
  TCanvas c1("c1", "asymmetry graph", 200,10,700,500);
  gStyle->SetOptStat(1111111);

  s13::ana::ScatteringAsymetryAlg asymmetry_alg{
    g_chain_up, g_chain_down, "4he", "-(S13: 4He asymmetry)", 5, 55, 70};

  c1.Clear();
  c1.Divide(3,3);
  for (int r = 3; r < 12; r++) {
    c1.cd(r-2);
    TCut target_cut{
      TString::Format("pow(tgt_up_xpos + 2, 2) + pow(tgt_up_ypos + 2, 2) < %.2f", 
                      std::pow(r, 2))};
    TString title = TString::Format("tgt(R<%i mm), "
                                    "#phi(+/-2.5 deg), #theta(+/-4 deg)", r);
    hdraw(g_chain_up, TString::Format("tgt_cut_r_%i", r), "p_theta*57.3:s1dc_theta*57.3",
          "(200,0.1,20,200,50,75)",
          "triggers[5]==1" && target_cut && 
          phi_corr_cut_1d && theta_cut_1d,
          title,
          "Fragment #theta angle [lab. deg.]",
          "Proton #theta angle [lab. deg.]", "colz");
    gPad->SetLogz();
  }
  c1.Print("out/ay_4he_r_tgt_cut.pdf(", "pdf");

  c1.Clear();
  c1.Divide(3,3);
  for (int r = 3; r < 12; r++) {
    std::cout << "Computing asymmetry with tgt R cut " << r << " mm..." << std::endl;
    asymmetry_alg.SetThetaCutWidth(8);
    asymmetry_alg.SetPhiCutWidth(5);
    asymmetry_alg.SetTargetCutRadius(r);
    asymmetry_alg.Run();

    TString title = TString::Format("tgt(R<%i mm), "
                                    "#phi(+/-2.5 deg), #theta(+/-4 deg)", r);
    std::vector<TGraph*> graphs = {Build4HeAyGraph(),
                                   BuildAyGraph(asymmetry_alg.GetResult(), title)};
    auto multi_graph = PrepareAyMultiGraph(title, graphs);

    c1.cd(r-2);
    multi_graph->Draw("ACP");
    //c1.BuildLegend(0.15, 0.67, 0.5, 0.88);
    //c1.BuildLegend();
  }
  c1.Print("out/ay_4he_r_tgt_cut.pdf)", "pdf");

  hist_out.Write();
  hist_out.Close();
}
