#include <iostream>

#include "TChain.h"


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


void phi_correlation_up_down() {
  TChain chain_up{"scattree"};
  TChain chain_down{"scattree"};

  // for (int i = 272; i < 297; i++){
  //for (int i = 272; i < 292; i++){
  for (int i = 272; i < 282; i++){
    TString name;
    name.Form("scattree_rootfiles/scattree_run%i.root", i);
    std::cout << "Adding file to chain_up: " << name << std::endl;
    chain_up.Add(name);
  }
  // for (int i = 297; i < 316; i++){
  //for (int i = 297; i < 316; i++){
  for (int i = 297; i < 306; i++){
    TString name;
    name.Form("scattree_rootfiles/scattree_run%i.root", i);
    std::cout << "Adding file to chain_down: " << name << std::endl;
    chain_down.Add(name);
  }


  TCut p_angle_sn_cut = "p_theta*57.3>65 && p_theta*57.3<68";
  TCut p_angle_sn_cut2 = "p_theta*57.3>58 && p_theta*57.3<61";

  TCut phi_corr_cut_1d = "abs(abs(p_phi*57.3-s1dc_phi*57.3) - 180) < 15";
  TCut target_cut{"tgt_up_xpos*tgt_up_xpos + tgt_up_ypos*tgt_up_ypos < 100"};
  TCut rect_target_cut{"tgt_up_xpos > -2 && tgt_up_xpos < 2 && tgt_up_ypos > -2 && tgt_up_ypos < 2"}

  TCut vertex_zpos_cut{"es_vertex_zpos > -45 && es_vertex_zpos < 45"};

  TCut theta_corr_cut_1d = "abs(89.22 - p_theta*57.3 - 2.86*s1dc_theta*57.3) < 5";

  TCut theta_corr_cut_1d_left_env = "89.22 - p_theta*57.3 - 4 < 2.65*s1dc_theta*57.3";
  TCut theta_corr_cut_1d_right_env = "89.22 - p_theta*57.3 + 4 > 3.12*s1dc_theta*57.3";


  TFile hist_out("out/phi_correlation_up_down.root", "RECREATE");
  TCanvas c1("c1");
  gStyle->SetOptStat(1111111);


  // target cut
  c1.Clear();
  c1.Divide(2,2);
  c1.cd(1);
  hdraw(chain_up, "tgt_hit_pol_up", "tgt_up_ypos:tgt_up_xpos",
        "(200,-40,40,200,-40,40)",
        "triggers[5]==1",
        "Target hit position - upstream extrapolation, pol. up",
        "X hit position [mm]",
        "Y hit position [mm]", "colz");
  c1.cd(2);
  hdraw(chain_up, "tgt_hit_pol_up_cut", "tgt_up_ypos:tgt_up_xpos",
        "(200,-40,40,200,-40,40)",
        "triggers[5]==1" && target_cut,
        "Target hit position - upstream extrapolation, pol. up, cut",
        "X hit position [mm]",
        "Y hit position [mm]", "colz");
  c1.cd(3);
  hdraw(chain_down, "tgt_hit_pol_down", "tgt_up_ypos:tgt_up_xpos",
        "(200,-40,40,200,-40,40)",
        "triggers[5]==1",
        "Target hit position - upstream extrapolation, pol. down",
        "X hit position [mm]",
        "Y hit position [mm]", "colz");
  c1.cd(4);
  hdraw(chain_down, "tgt_hit_pol_down_cut", "tgt_up_ypos:tgt_up_xpos",
        "(200,-40,40,200,-40,40)",
        "triggers[5]==1" && target_cut,
        "Target hit position - upstream extrapolation, pol. down, cut",
        "X hit position [mm]",
        "Y hit position [mm]", "colz");
  c1.Print("out/phi_correlation_up_down.pdf(", "pdf");

  // vert-Z cut 1D, pol. up
  c1.Clear();
  c1.Divide(2,2);
  c1.cd(1);
  hdraw(chain_up, "esl_vertz_pos_pol_up", "esl_vertex_zpos",
        "(200,-90,90)",
        "triggers[5]==1" && target_cut,
        "Vertex Z-pos - ESPRI left, pol. up",
        "Z-pos [mm]",
        "Counts", "colz");
  c1.cd(2);
  hdraw(chain_up, "esr_vertz_pos_pol_up", "esr_vertex_zpos",
        "(200,-90,90)",
        "triggers[5]==1" && target_cut,
        "Vertex Z-pos - ESPRI right, pol. up",
        "Z-pos [mm]",
        "Counts", "colz");
  c1.cd(3);
  hdraw(chain_up, "esl_vertz_pos_pol_up_cut", "esl_vertex_zpos",
        "(200,-90,90)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut,
        "Vertex Z-pos - ESPRI left, pol. up, cut",
        "Z-pos [mm]",
        "Counts", "colz");
  c1.cd(4);
  hdraw(chain_up, "esr_vertz_pos_pol_up_cut", "esr_vertex_zpos",
        "(200,-90,90)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut,
        "Vertex Z-pos - ESPRI right, pol. up, cut",
        "Z-pos [mm]",
        "Counts", "colz");
  c1.Print("out/phi_correlation_up_down.pdf", "pdf");

  // vert-Z cut 1D, pol. down
  c1.Clear();
  c1.Divide(2,2);
  c1.cd(1);
  hdraw(chain_down, "esl_vertz_pos_pol_down", "esl_vertex_zpos",
        "(200,-90,90)",
        "triggers[5]==1" && target_cut,
        "Vertex Z-pos - ESPRI left, pol. down",
        "Z-pos [mm]",
        "Counts", "colz");
  c1.cd(2);
  hdraw(chain_down, "esr_vertz_pos_pol_down", "esr_vertex_zpos",
        "(200,-90,90)",
        "triggers[5]==1" && target_cut,
        "Vertex Z-pos - ESPRI right, pol. down",
        "Z-pos [mm]",
        "Counts", "colz");
  c1.cd(3);
  hdraw(chain_down, "esl_vertz_pos_pol_down_cut", "esl_vertex_zpos",
        "(200,-90,90)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut,
        "Vertex Z-pos - ESPRI left, pol. down, cut",
        "Z-pos [mm]",
        "Counts", "colz");
  c1.cd(4);
  hdraw(chain_down, "esr_vertz_pos_pol_up_cut", "esr_vertex_zpos",
        "(200,-90,90)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut,
        "Vertex Z-pos - ESPRI right, pol. down, cut",
        "Z-pos [mm]",
        "Counts", "colz");
  c1.Print("out/phi_correlation_up_down.pdf", "pdf");

  // phi corr. cut 1D
  c1.Clear();
  c1.Divide(2,2);
  c1.cd(1);
  hdraw(chain_up, "polar_corr_pol_up_no_phi_cut", "p_theta*57.3:s1dc_theta*57.3",
        "(200,0.1,20,200,50,75)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut,
        "Polar angle correlation, pol. up - tgt, vertZ cuts",
        "Fragment #theta angle [lab. deg.]",
        "Proton #theta angle [lab. deg.]", "colz");
  c1.cd(2);
  hdraw(chain_up, "polar_corr_phi_cut", "p_theta*57.3:s1dc_theta*57.3",
        "(200,0.1,20,200,50,75)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut && phi_corr_cut_1d,
        "Polar angle correlation, pol. up - tgt, vertZ, #phi cuts",
        "Fragment #theta angle [lab. deg.]",
        "Proton #theta angle [lab. deg.]", "colz");
  c1.cd(3);
  hdraw(chain_down, "polar_corr_pol_down_no_phi_cut", "p_theta*57.3:s1dc_theta*57.3",
        "(200,0.1,20,200,50,75)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut,
        "Polar angle correlation, pol. down - tgt, vertZ cuts",
        "Fragment #theta angle [lab. deg.]",
        "Proton #theta angle [lab. deg.]", "colz");
  c1.cd(4);
  hdraw(chain_down, "polar_corr_pol_down_phi_cut", "p_theta*57.3:s1dc_theta*57.3",
        "(200,0.1,20,200,50,75)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut && phi_corr_cut_1d,
        "Polar angle correlation, pol. down - tgt, vertZ, #phi cuts",
        "Fragment #theta angle [lab. deg.]",
        "Proton #theta angle [lab. deg.]", "colz");
  c1.Print("out/phi_correlation_up_down.pdf", "pdf");

  // pol up & down corr.
  c1.Clear();
  c1.Divide(2);
  c1.cd(1);
  hdraw(chain_up, "polar_up", "p_theta*57.3:s1dc_theta*57.3",
        "(200,0,20,200,50,75)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut &&
        phi_corr_cut_1d,
        "Polar angle correlation - pol. up runs",
        "Fragment #theta angle [lab. deg.]",
        "Proton #theta angle [lab. deg.]", "colz");
  c1.cd(2);
  hdraw(chain_down, "polar_down", "p_theta*57.3:s1dc_theta*57.3",
        "(200,0,20,200,50,75)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut &&
        phi_corr_cut_1d,
        "Polar angle correlation - pol. down runs",
        "Fragment #theta angle [lab. deg.]",
        "Proton #theta angle [lab. deg.]", "colz");
  c1.Print("out/phi_correlation_up_down.pdf", "pdf");

  // S/N
  c1.Clear();
  c1.Divide(1,2);
  c1.cd(1);
  hdraw(chain_up, "polar_up_sn", "s1dc_theta*57.3",
        "(200,0,20)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut &&
        phi_corr_cut_1d && p_angle_sn_cut,
        "Polar angle cut S/N - pol. up runs",
        "Fragment #theta angle [lab. deg.]",
        "Counts", "colz");
  c1.cd(2);
  hdraw(chain_down, "polar_down_sn", "s1dc_theta*57.3",
        "(200,0,20)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut &&
        phi_corr_cut_1d && p_angle_sn_cut,
        "Polar angle cut S/N - pol. down runs",
        "Fragment #theta angle [lab. deg.]",
        "Counts", "colz");
  c1.Print("out/phi_correlation_up_down.pdf", "pdf");

  // pol up corr. with theta corr cut
  c1.Clear();
  c1.Divide(2);
  c1.cd(1);
  hdraw(chain_up, "polar_up_cut", "p_theta*57.3:s1dc_theta*57.3",
        "(200,0,20,200,50,75)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut &&
        phi_corr_cut_1d && 
        theta_corr_cut_1d_left_env && theta_corr_cut_1d_right_env,
        "Polar angle correlation - pol. up runs",
        "Fragment #theta angle [lab. deg.]",
        "Proton #theta angle [lab. deg.]", "colz");
  c1.cd(2);
  hdraw(chain_down, "polar_down_cut", "p_theta*57.3:s1dc_theta*57.3",
        "(200,0,20,200,50,75)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut &&
        phi_corr_cut_1d && 
        theta_corr_cut_1d_left_env && theta_corr_cut_1d_right_env,
        "Polar angle correlation - pol. down runs",
        "Fragment #theta angle [lab. deg.]",
        "Proton #theta angle [lab. deg.]", "colz");
  c1.Print("out/phi_correlation_up_down.pdf", "pdf");

  // pol up: left & right ESPRIs
  c1.Clear();
  c1.Divide(2,2);
  c1.cd(1);
  hdraw(chain_up, "polar_up_esl", "esl_p_theta*57.3:s1dc_theta*57.3",
        "(200,0,20,200,50,75)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut &&
        phi_corr_cut_1d,
        "Polar angle correlation - pol. up, ESPRI left",
        "Fragment #theta angle [lab. deg.]",
        "Proton #theta angle [lab. deg.]", "colz");
  c1.cd(2);
  hdraw(chain_up, "polar_up_esr", "esr_p_theta*57.3:s1dc_theta*57.3",
        "(200,0,20,200,50,75)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut &&
        phi_corr_cut_1d,
        "Polar angle correlation - pol. up, ESPRI right",
        "Fragment #theta angle [lab. deg.]",
        "Proton #theta angle [lab. deg.]", "colz");
  c1.cd(3);
  hdraw(chain_down, "polar_down_esl", "esl_p_theta*57.3:s1dc_theta*57.3",
        "(200,0,20,200,50,75)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut &&
        phi_corr_cut_1d,
        "Polar angle correlation - pol. down, ESPRI left",
        "Fragment #theta angle [lab. deg.]",
        "Proton #theta angle [lab. deg.]", "colz");
  c1.cd(4);
  hdraw(chain_down, "polar_down_esr", "esr_p_theta*57.3:s1dc_theta*57.3",
        "(200,0,20,200,50,75)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut &&
        phi_corr_cut_1d,
        "Polar angle correlation - pol. down, ESPRI right",
        "Fragment #theta angle [lab. deg.]",
        "Proton #theta angle [lab. deg.]", "colz");
  c1.Print("out/phi_correlation_up_down.pdf", "pdf");

  // pol up: left ESPRI S/N
  c1.Clear();
  c1.Divide(4,2);
  for (int i = 55; i < 68; i += 2) {
    c1.cd((i - 53) / 2);
    TCut sn_cut{TString::Format("esl_p_theta*57.3>%i && p_theta*57.3<%i", i, i + 2)};
    hdraw(chain_up, TString::Format("polar_up_esl_sn_%i", (i - 53) / 2), "s1dc_theta*57.3",
          "(200,0,20)",
          "triggers[5]==1" && target_cut && vertex_zpos_cut &&
          phi_corr_cut_1d && sn_cut,
          TString::Format("Polar angle cut S/N - pol. up, ESPRI left (#theta_p : %i - %i deg.)", i, i + 2),
          "Fragment #theta angle [lab. deg.]",
          "Counts", "colz");
  }
  c1.Print("out/phi_correlation_up_down.pdf", "pdf");

  // pol up: right ESPRI S/N
  c1.Clear();
  c1.Divide(4,2);
  for (int i = 55; i < 68; i += 2) {
    c1.cd((i - 53) / 2);
    TCut sn_cut{TString::Format("esr_p_theta*57.3>%i && p_theta*57.3<%i", i, i + 2)};
    hdraw(chain_up, TString::Format("polar_up_esr_sn_%i", (i - 53) / 2), "s1dc_theta*57.3",
          "(200,0,20)",
          "triggers[5]==1" && target_cut && vertex_zpos_cut &&
          phi_corr_cut_1d && sn_cut,
          TString::Format("Polar angle cut S/N - pol. up, ESPRI right (#theta_p : %i - %i deg.)", i, i + 2),
          "Fragment #theta angle [lab. deg.]",
          "Counts", "colz");
  }
  c1.Print("out/phi_correlation_up_down.pdf", "pdf");

  // pol down: left & right ESPRIs
  c1.Clear();
  c1.Divide(2);
  c1.cd(1);
  hdraw(chain_down, "polar_up_esl", "esl_p_theta*57.3:s1dc_theta*57.3",
        "(200,0,20,200,50,75)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut &&
        phi_corr_cut_1d,
        "Polar angle correlation - pol. down, ESPRI left",
        "Fragment #theta angle [lab. deg.]",
        "Proton #theta angle [lab. deg.]", "colz");
  c1.cd(2);
  hdraw(chain_down, "polar_up_esr", "esr_p_theta*57.3:s1dc_theta*57.3",
        "(200,0,20,200,50,75)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut &&
        phi_corr_cut_1d,
        "Polar angle correlation - pol. down, ESPRI right",
        "Fragment #theta angle [lab. deg.]",
        "Proton #theta angle [lab. deg.]", "colz");
  c1.Print("out/phi_correlation_up_down.pdf", "pdf");

  // pol down: left ESPRI S/N
  c1.Clear();
  c1.Divide(4,2);
  for (int i = 55; i < 68; i += 2) {
    c1.cd((i - 53) / 2);
    TCut sn_cut{TString::Format("esl_p_theta*57.3>%i && p_theta*57.3<%i", i, i + 2)};
    hdraw(chain_down, TString::Format("polar_up_esl_sn_%i", (i - 53) / 2), "s1dc_theta*57.3",
          "(200,0,20)",
          "triggers[5]==1" && target_cut && vertex_zpos_cut &&
          phi_corr_cut_1d && sn_cut,
          TString::Format("Polar angle cut S/N - pol. down, ESPRI left (#theta_p : %i - %i deg.)", i, i + 2),
          "Fragment #theta angle [lab. deg.]",
          "Counts", "colz");
  }
  c1.Print("out/phi_correlation_up_down.pdf", "pdf");

  // pol down: right ESPRI S/N
  c1.Clear();
  c1.Divide(4,2);
  for (int i = 55; i < 68; i += 2) {
    c1.cd((i - 53) / 2);
    TCut sn_cut{TString::Format("esr_p_theta*57.3>%i && p_theta*57.3<%i", i, i + 2)};
    hdraw(chain_down, TString::Format("polar_up_esr_sn_%i", (i - 53) / 2), "s1dc_theta*57.3",
          "(200,0,20)",
          "triggers[5]==1" && target_cut && vertex_zpos_cut &&
          phi_corr_cut_1d && sn_cut,
          TString::Format("Polar angle cut S/N - pol. down, ESPRI right (#theta_p : %i - %i deg.)", i, i + 2),
          "Fragment #theta angle [lab. deg.]",
          "Counts", "colz");
  }
  c1.Print("out/phi_correlation_up_down.pdf)", "pdf");

  // locus angle determination
  // c1.Clear();
  // c1.Divide(2,2);
  // c1.cd(1);
  // hdraw(chain_up, "polar_up_sn_10", "s1dc_theta*57.3",
  //       "(200,0,20)",
  //       "triggers[5]==1" && target_cut && vertex_zpos_cut &&
  //       phi_corr_cut_1d && p_angle_sn_cut,
  //       "Polar angle cut S/N - pol. up runs",
  //       "Fragment #theta angle [lab. deg.]",
  //       "Counts", "colz");
  // c1.cd(3);
  // c1.cd(3);
  // hdraw(chain_up, "polar_up_sn_11", "s1dc_theta*57.3",
  //       "(200,0,20)",
  //       "triggers[5]==1" && target_cut && vertex_zpos_cut &&
  //       phi_corr_cut_1d && p_angle_sn_cut2,
  //       "Polar angle cut S/N - pol. up runs",
  //       "Fragment #theta angle [lab. deg.]",
  //       "Counts", "colz");
  // c1.cd(2);
  // hdraw(chain_down, "polar_down_sn_20", "s1dc_theta*57.3",
  //       "(200,0,20)",
  //       "triggers[5]==1" && target_cut && vertex_zpos_cut &&
  //       phi_corr_cut_1d && p_angle_sn_cut,
  //       "Polar angle cut S/N - pol. down runs",
  //       "Fragment #theta angle [lab. deg.]",
  //       "Counts", "colz");
  // c1.cd(4);
  // hdraw(chain_down, "polar_down_sn_21", "s1dc_theta*57.3",
  //       "(200,0,20)",
  //       "triggers[5]==1" && target_cut && vertex_zpos_cut &&
  //       phi_corr_cut_1d && p_angle_sn_cut2,
  //       "Polar angle cut S/N - pol. down runs",
  //       "Fragment #theta angle [lab. deg.]",
  //       "Counts", "colz");
  // c1.Print("out/phi_correlation_up_down.pdf", "pdf");

  hist_out.Write();
  hist_out.Close();
}
