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


void beam_timing() {
    TChain chain{"scattree"};

    //for (int i = 272; i < 283; i++){
    for (int i = 272; i < 276; i++){
        TString name;
        name.Form("scattree_rootfiles/scattree_run%i.root", i);
        std::cout << "Adding file to chain: " << name << std::endl;
        chain.Add(name);
    }

    
    TCut target_cut{"tgt_up_xpos*tgt_up_xpos + tgt_up_ypos*tgt_up_ypos < 64"};
    TCut vertex_zpos_cut{"es_vertex_zpos > -30 && es_vertex_zpos < 30"};
    TCut esl_phi_corr_cut = "abs(-1.3 * s1dc_phi*57.3 - 25 - esl_p_phi*57.3) < 10";
    TCut esr_phi_corr_cut = "abs(-1.35 * s1dc_phi*57.3 + 33 - esr_p_phi*57.3) < 10";

    TCut beam_timing_cut{"f7_t > 680 && f7_t < 720"};

    TFile hist_out("out/phi_angle_correlation.root", "RECREATE");
    TCanvas *c1 = new TCanvas();
    gStyle->SetOptStat(1111111);
 
    c1->Clear();
    c1->Divide(1,2);
    c1->cd(1);
    hdraw(chain, "f7_t", "f7_t", "(400,-200,1200)", "triggers[5]==1", 
          "F7 timing", "Timing [ch]", "Counts", "colz");
    gPad->SetLogy();
    c1->cd(2);
    hdraw(chain, "f7_q", "f7_q", "(400,-200,1200)", "triggers[5]==1", 
          "F7 dE", "dE [ch]", "Counts", "colz");
    gPad->SetLogy();
    c1->Print("out/beam_timing.pdf(", "pdf");

    c1->Clear();
    c1->Divide(1,2);
    c1->cd(1);
    hdraw(chain, "f7_t_cut", "f7_t", "(400,-200,1200)", 
          "triggers[5]==1" && beam_timing_cut, 
          "F7 timing", "Timing [ch]", "Counts", "colz");
    gPad->SetLogy();
    c1->cd(2);
    hdraw(chain, "f7_q_cut", "f7_q", "(400,-200,1200)", 
          "triggers[5]==1" && beam_timing_cut, 
          "F7 dE", "dE [ch]", "Counts", "colz");
    gPad->SetLogy();
    c1->Print("out/beam_timing.pdf", "pdf");

    c1->Clear();
    c1->Divide(1,2);
    c1->cd(1);
    hdraw(chain, "sbt1_t", "sbt1_t", "(400,-200,1200)", "triggers[5]==1", 
          "SBT1 timing", "Timing [ch]", "Counts");
    gPad->SetLogy();
    c1->cd(2);
    hdraw(chain, "sbt1_q", "sbt1_q", "(400,-200,1200)", "triggers[5]==1", 
          "SBT1 dE", "dE [ch]", "Counts");
    gPad->SetLogy();
    c1->Print("out/beam_timing.pdf", "pdf");

    c1->Clear();
    c1->Divide(1,2);
    c1->cd(1);
    hdraw(chain, "sbt2_t", "sbt2_t", "(400,-200,1200)", "triggers[5]==1", 
          "SBT2 timing", "Timing [ch]", "Counts");
    gPad->SetLogy();
    c1->cd(2);
    hdraw(chain, "sbt2_q", "sbt2_q", "(400,-200,1200)", "triggers[5]==1", 
          "SBT2 dE", "dE [ch]", "Counts");
    gPad->SetLogy();
    c1->Print("out/beam_timing.pdf", "pdf");

    c1->Clear();
    c1->Divide(1,2);
    c1->cd(1);
    hdraw(chain, "sbv_t", "sbv_t", "(400,-200,1200)", "triggers[5]==1", 
          "SBV timing", "Timing [ch]", "Counts");
    gPad->SetLogy();
    c1->cd(2);
    hdraw(chain, "sbv_q", "sbv_q", "(400,-200,1200)", "triggers[5]==1", 
          "SBV dE", "dE [ch]", "Counts");
    gPad->SetLogy();
    c1->Print("out/beam_timing.pdf", "pdf");

    c1->Clear();
    c1->Divide(1);
    c1->cd(1);
    hdraw(chain, "polar1", "p_theta*57.3:s1dc_theta*57.3",
          "(200,0,20,200,50,75)",
          "triggers[5]==1" && target_cut && vertex_zpos_cut &&
          (esr_phi_corr_cut || esl_phi_corr_cut),
          "Polar angle correlation: tgt, vertZ, $phi cuts",
          "Fragment #theta angle [lab. deg.]",
          "Proton #theta angle [lab. deg.]", "colz");
    c1->Print("out/beam_timing.pdf", "pdf");

    c1->Clear();
    c1->Divide(1);
    c1->cd(1);
    hdraw(chain, "polar2", "p_theta*57.3:s1dc_theta*57.3",
          "(200,0,20,200,50,75)",
          "triggers[5]==1" && target_cut && vertex_zpos_cut &&
          (esr_phi_corr_cut || esl_phi_corr_cut) && beam_timing_cut,
          "Polar angle correlation: tgt, vertZ, #phi, btiming cuts",
          "Fragment #theta angle [lab. deg.]",
          "Proton #theta angle [lab. deg.]", "colz");
    c1->Print("out/beam_timing.pdf", "pdf");

    c1->Clear();
    c1->Divide(1);
    c1->cd(1);
    hdraw(chain, "polar3", "p_theta*57.3:s1dc_theta*57.3",
          "(200,0,20,200,50,75)",
          "triggers[5]==1" && target_cut && beam_timing_cut,
          "Polar angle correlation: tgt, btiming cuts",
          "Fragment #theta angle [lab. deg.]",
          "Proton #theta angle [lab. deg.]", "colz");
    c1->Print("out/beam_timing.pdf)", "pdf");

    
    hist_out.Write();
    hist_out.Close();
}
