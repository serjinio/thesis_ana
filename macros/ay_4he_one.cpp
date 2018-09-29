
#include "TChain.h"
#include "TGraph.h"

#include "ay_init_ds.hpp"
#include "anapow.hpp"


using namespace s13::ana;


static TCut g_p_angle_sn_cut = "p_theta*57.3>65 && p_theta*57.3<68";

static TCut g_phi_corr_cut_1d = "abs(abs(p_phi*57.3-s1dc_phi*57.3) - 180) < 15";
static TCut g_target_cut{"tgt_up_xpos*tgt_up_xpos + tgt_up_ypos*tgt_up_ypos < 100"};
static TCut g_vertex_zpos_cut{"es_vertex_zpos > -45 && es_vertex_zpos < 45"};


ScatteringYield compute_scattering_asymmetry(TString name, TString title) {
  s13::ana::ScatteringAsymetryAlg asymmetry_alg{
    g_chain_up, g_chain_down, name, title, 6, 55, 70};
  
  asymmetry_alg.Run();
  
  return asymmetry_alg.GetResult();
}

void ay_4he_one() {
  init_dataset();

  TFile hist_out("out/ay_4he_one.root", "RECREATE");
  TCanvas c1("c1", "asymmetry graph", 200,10,700,500);
  //gStyle->SetOptStat(1111);

  c1.Clear();
  c1.Divide(1);

  s13::ana::ScatteringAsymetryAlg asymmetry_alg{
    g_chain_up, g_chain_down, "4he", "-(S13: 4He asymmetry)", 5, 55, 70};
  asymmetry_alg.SetThetaCutWidth(8);
  asymmetry_alg.SetPhiCutWidth(5);
  asymmetry_alg.SetTargetCutRadius(6);
  asymmetry_alg.Run();

  TString title = TString::Format("tgt(R<6 mm), "
                                  "#phi(+/-2.5 deg), #theta(+/-4 )");
  std::vector<TGraph*> graphs = {Build4HeAyGraph(),
                                 BuildAyGraph(asymmetry_alg.GetResult(), title)};
  auto multi_graph = PrepareAyMultiGraph(title, graphs);

  multi_graph->Draw("ACP");
  c1.BuildLegend(0.15, 0.67, 0.5, 0.88);
  //c1.BuildLegend();
  c1.Print("out/ay_4he_one.pdf", "pdf");

  hist_out.Write();
  hist_out.Close();
}
