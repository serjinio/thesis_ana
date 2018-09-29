
#include "TChain.h"
#include "TGraph.h"

#include "ay_init_ds.hpp"
#include "anapow.hpp"


using namespace s13::ana;

static TCut g_vertex_zpos_cut{"es_vertex_zpos > -45 && es_vertex_zpos < 45"};
static TCut g_p_angle_sn_cut = "p_theta*57.3>65 && p_theta*57.3<68";


TCut phi_corr_cut_1d = "abs(abs(p_phi*57.3-s1dc_phi*57.3) - 180) < 10";
TCut theta_corr_cut_1d_left_env = "89.22 - p_theta*57.3 - 4 < 2.65*s1dc_theta*57.3";
TCut theta_corr_cut_1d_right_env = "89.22 - p_theta*57.3 + 4 > 3.12*s1dc_theta*57.3";
TCut theta_cut_1d = theta_corr_cut_1d_left_env && theta_corr_cut_1d_right_env;

TCut g_target_cut{"pow(tgt_up_xpos + 2, 2) + pow(tgt_up_ypos + 2, 2) < 64"};


void ay_4he_tgt_cut_search() {
  init_dataset();

  Tfile hist_out("out/ay_4he_tgt_cut_search.root", "RECREATE");
  TCanvas c1("c1", "Graph example", 200,10,700,500);
  //gStyle->SetOptStat(1111);

  c1.Clear();
  c1.Divide(4,4);

  s13::ana::ScatteringAsymetryAlg asymmetry_alg{
    g_chain_up, g_chain_down, "4he", "S13: 4He asymmetry", 3, 55, 70};

  // mesh seach on target cut
  const int k_mesh_step = 4;
  int hist_counter = 1;
  for (int i_y = 3; i_y > -1; i_y--) {
    int ymin = k_mesh_step * i_y - 8 - 2;
    int ymax = k_mesh_step * i_y - 8 - 2 + k_mesh_step;
    for (int i_x = 0; i_x < 4; i_x++) {
      int xmin = k_mesh_step * i_x - 8 - 2;
      int xmax = k_mesh_step * i_x - 8 - 2 + k_mesh_step;
      
      std::cout << "Computing asymmetry with target cut " << 
        "(xmin; ymin) - (xmax; ymax): " << \
        "(" << xmin << "; " << ymin << ") - (" << xmax << "; " << ymax << ")" << \
        std::endl;
      asymmetry_alg.SetThetaCutWidth(8);
      asymmetry_alg.SetPhiCutWidth(5);
      asymmetry_alg.SetTargetCut(xmin, ymin, xmax, ymax);
      asymmetry_alg.Run();

      TString title = TString::Format("Tgt cut: %i,%i - %i,%i", xmin, ymin, xmax, ymax);
      std::vector<TGraph*> graphs = {Build4HeAyGraph(), 
                                     BuildAyGraph(asymmetry_alg.GetResult(), title)};
      auto multi_graph = PrepareAyMultiGraph(title, graphs);

      c1.cd(hist_counter);
      multi_graph->Draw("ACP");
      //c1.BuildLegend(0.15, 0.67, 0.5, 0.88);
      //c1.BuildLegend();
      hist_counter++;
    }
  }
  c1.Print("out/ay_4he_tgt_cut_search.pdf", "pdf");

  hist_out.Write();
  hist_out.Close();
}
