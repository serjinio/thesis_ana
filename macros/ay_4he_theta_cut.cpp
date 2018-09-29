
#include "TChain.h"
#include "TGraph.h"

#include "ay_init_ds.hpp"
#include "anapow.hpp"


using namespace s13::ana;


TCut vertex_zpos_cut{"es_vertex_zpos > -45 && es_vertex_zpos < 45"};


TCut target_cut{"pow(tgt_up_xpos + 2, 2) + pow(tgt_up_ypos + 2, 2) < 36"};
TCut phi_corr_cut_1d{"abs(abs(p_phi*57.3-s1dc_phi*57.3) - 180) < 2.5"};

TCut theta_cut_1d2{"abs(s1dc_theta*57.3 - he_theta_theor*57.3) < 3"};

TCut proton_1d_cut{"p_theta*57.3>55 && p_theta*57.3<70"};


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

void ay_4he_theta_cut() {
  init_dataset();

  TFile hist_out("out/ay_4he_theta_cut.root", "RECREATE");
  TCanvas c1("c1", "asymmetry graph", 200,10,700,500);
  gStyle->SetOptStat(1111111);

  s13::ana::ScatteringAsymetryAlg asymmetry_alg{
    g_chain_up, g_chain_down, "4he", "-(S13: 4He asymmetry)", 5, 55, 70};

  c1.Clear();
  c1.Divide(2,2);
  c1.cd(1);
  hdraw(g_chain_up, "esl_theta_cut", "esl_p_theta*57.3:s1dc_theta*57.3 - he_theta_theor*57.3",
        "(200,-10,10,200,50,75)",
        "triggers[5]==1" && target_cut && 
        phi_corr_cut_1d,
        "EPRI left: #theta exp. - #theta theoretical",
        "Fragment #theta angle [lab. deg.]",
        "Proton #theta angle [lab. deg.]", "colz");
  // TLine *line1 = new TLine(-3,50,-3,75);
  // line1->SetLineColor(kRed);
  // line1->Draw();
  // TLine *line2 = new TLine(3,50,3,75);
  // line2->SetLineColor(kRed);
  // line2->Draw();
  gPad->SetLogz();
  c1.cd(2);
  hdraw(g_chain_up, "esr_theta_cut", "esr_p_theta*57.3:s1dc_theta*57.3 - he_theta_theor*57.3",
        "(200,-10,10,200,50,75)",
        "triggers[5]==1" && target_cut && 
        phi_corr_cut_1d && theta_cut_1d2,
        "ESPRI right: #theta exp. - #theta theoretical",
        "Fragment #theta angle [lab. deg.]",
        "Proton #theta angle [lab. deg.]", "colz");
  gPad->SetLogz();
  c1.cd(3);
  hdraw(g_chain_up, "esl_theta_cut_1d", "s1dc_theta*57.3 - he_theta_theor*57.3",
        "(200,-10,10)",
        "triggers[5]==1" && target_cut && 
        phi_corr_cut_1d && proton_1d_cut,
        "EPRI left: #theta exp. - #theta theoretical 1D",
        "Fragment #theta angle [lab. deg.]",
        "Counts", "colz");
  c1.cd(4);
  hdraw(g_chain_up, "esr_theta_cut_1d", "s1dc_theta*57.3 - he_theta_theor*57.3",
        "(200,-10,10)",
        "triggers[5]==1" && target_cut && 
        phi_corr_cut_1d && theta_cut_1d2 && proton_1d_cut,
        "ESPRI right: #theta exp. - #theta theoretical 1D",
        "Fragment #theta angle [lab. deg.]",
        "Counts", "colz");
  c1.Print("out/ay_4he_theta_cut.pdf(", "pdf");

  c1.Clear();
  c1.Divide(3,3);
  for (int w = 2; w < 20; w += 2) {
    c1.cd((w+1) / 2);
    TCut theta_cut_1d3 = TString::Format("abs(s1dc_theta*57.3 - "
                                         "he_theta_theor*57.3) < %.2f",
                                         w/2.).Data();
    TString title = TString::Format("tgt(R<6 mm), "
                                    "#phi(+/-2.5 deg), #theta(+/-%.2f deg)", w/2.);
    hdraw(g_chain_up, TString::Format("theta_cut_w%i", w), 
          "p_theta*57.3:s1dc_theta*57.3",
          "(200,0.1,20,200,50,75)",
          "triggers[5]==1" && target_cut && 
          phi_corr_cut_1d && theta_cut_1d3,
          title,
          "Fragment #theta angle [lab. deg.]",
          "Proton #theta angle [lab. deg.]", "colz");
    gPad->SetLogz();
  }
  c1.Print("out/ay_4he_theta_cut.pdf", "pdf");

  c1.Clear();
  c1.Divide(3,3);
  for (int w = 2; w < 20; w += 2) {
    c1.cd((w+1) / 2);
    TCut theta_cut_1d3 = TString::Format("abs(s1dc_theta*57.3 - "
                                         "he_theta_theor*57.3) < %.2f",
                                         w/2.).Data();
    TString title = TString::Format("tgt(R<6 mm), "
                                    "#phi(+/-2.5 deg), #theta(+/-%.2f deg)", w/2.);
    hdraw(g_chain_up, TString::Format("theta_cut_1d_w%i", w), 
          "s1dc_theta*57.3 - he_theta_theor*57.3",
          "(200,-10,10)",
          "triggers[5]==1" && target_cut && 
          phi_corr_cut_1d && theta_cut_1d3 && proton_1d_cut,
          title,
          "Fragment #theta angle [lab. deg.]",
          "Proton #theta angle [lab. deg.]", "colz");
    gPad->SetLogz();
  }
  c1.Print("out/ay_4he_theta_cut.pdf)", "pdf");

  // c1.Clear();
  // c1.Divide(3,3);
  // for (int w = 4; w < 40; w += 4) {
  //   std::cout << "Computing asymmetry with theta cut width " << w << \
  //     " deg..." << std::endl;
  //   asymmetry_alg.SetThetaCutWidth(w);
  //   asymmetry_alg.SetPhiCutWidth(5);
  //   asymmetry_alg.SetTargetCutRadius(6);
  //   asymmetry_alg.Run();

  //   TString title = TString::Format("tgt(R<6 mm), "
  //                                   "#phi(+/-2.5 deg), #theta(+/-%.2f deg)", w/2.);
  //   std::vector<TGraph*> graphs = {Build4HeAyGraph(),
  //                                  BuildAyGraph(asymmetry_alg.GetResult(), title)};
  //   auto multi_graph = PrepareAyMultiGraph(title, graphs);

  //   c1.cd((w+1) / 4);
  //   multi_graph->Draw("ACP");
  //   //c1.BuildLegend(0.15, 0.67, 0.5, 0.88);
  //   //c1.BuildLegend();
  // }
  // c1.Print("out/ay_4he_theta_cut.pdf)", "pdf");

  hist_out.Write();
  hist_out.Close();
}
