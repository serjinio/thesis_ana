#include <iostream>

#include "TChain.h"


void beam_tracking_dcs() {
    TChain chain{"scattree"};

    for (int i = 272; i < 274; i++){
        TString name;
        name.Form("scattree_rootfiles/scattree_run%i.root", i);
        std::cout << "Adding file to chain: " << name << std::endl;
        chain.Add(name);
    }

    chain.Print();

    TCanvas c1("c1");
    c1.Divide(3,2);

    c1.cd(1);
    gPad->SetLogz();
    chain.Draw("bdc1_ypos:bdc1_xpos >>hbdc1(400,-40,40,400,-40,40)",
               "triggers[5]==1", "colz");
    c1.cd(2);
    gPad->SetLogz();
    chain.Draw("bdc2_ypos:bdc2_xpos >>hbdc2(400,-40,40,400,-40,40)",
               "triggers[5]==1", "colz");
    c1.cd(3);
    gPad->SetLogz();
    chain.Draw("tgt_up_ypos:tgt_up_xpos >>hbdc3(400,-40,40,400,-40,40)",
               "triggers[5]==1", "colz");
    c1.cd(4);
    gPad->SetLogz();
    chain.Draw("s1dc_ypos:s1dc_xpos >>hdsdc1(400,-40,40,400,-40,40)",
               "triggers[5]==1", "colz");
    c1.cd(5);
    gPad->SetLogz();
    chain.Draw("fdc0_ypos:fdc0_xpos >>hdsdc2(400,-40,40,400,-40,40)",
               "triggers[5]==1", "colz");
    c1.cd(6);
    gPad->SetLogz();
    chain.Draw("tgt_dwn_ypos:tgt_dwn_xpos >>hdsdc3(400,-40,40,400,-40,40)",
               "triggers[5]==1", "colz");
    c1.Print("out/beam_tracking_dcs.pdf(", "pdf");

    c1.Clear();
    c1.Divide(3,2);

    c1.cd(1);
    //gPad->SetLogy();
    chain.Draw("bdc1_xpos >>hbeam1(100,-40,40)",
               "triggers[5]==1", "colz");
    c1.cd(2);
    //gPad->SetLogy();
    chain.Draw("bdc2_xpos >>hbeam3(100,-40,40)",
               "triggers[5]==1", "colz");
    c1.cd(3);
    //gPad->SetLogy();
    chain.Draw("tgt_up_xpos >>hbeam5(100,-40,40)",
               "triggers[5]==1", "colz");
    c1.cd(4);
    //gPad->SetLogy();
    chain.Draw("bdc1_ypos >>hbeam2(100,-40,40)",
               "triggers[5]==1", "colz");
    c1.cd(5);
    //gPad->SetLogy();
    chain.Draw("bdc2_ypos >>hbeam4(100,-40,40)",
               "triggers[5]==1", "colz");
    c1.cd(6);
    //gPad->SetLogy();
    chain.Draw("tgt_up_ypos >>hbeam6(100,-40,40)",
               "triggers[5]==1", "colz");
    c1.Print("out/beam_tracking_dcs.pdf", "pdf");

    c1.Clear();
    c1.Divide(3,2);

    c1.cd(1);
    //gPad->SetLogy();
    chain.Draw("s1dc_xpos >>hbeam7(100,-40,40)",
               "triggers[5]==1", "colz");
    c1.cd(2);
    //gPad->SetLogy();
    chain.Draw("fdc0_xpos >>hbeam8(100,-40,40)",
               "triggers[5]==1", "colz");
    c1.cd(3);
    //gPad->SetLogy();
    chain.Draw("tgt_dwn_xpos >>hbeam9(100,-40,40)",
               "triggers[5]==1", "colz");
    c1.cd(4);
    //gPad->SetLogy();
    chain.Draw("s1dc_ypos >>hbeam10(100,-40,40)",
               "triggers[5]==1", "colz");
    c1.cd(5);
    //gPad->SetLogy();
    chain.Draw("fdc0_ypos >>hbeam11(100,-40,40)",
               "triggers[5]==1", "colz");
    c1.cd(6);
    //gPad->SetLogy();
    chain.Draw("tgt_dwn_ypos >>hbeam12(100,-40,40)",
               "triggers[5]==1", "colz");
    c1.Print("out/beam_tracking_dcs.pdf", "pdf");

    c1.Clear();
    c1.Divide(2,2);

    c1.cd(1);
    gPad->SetLogy();
    chain.Draw("tgt_up_aang >>hbeam13(100,-0.4,0.4)",
               "triggers[5]==1", "colz");
    c1.cd(2);
    gPad->SetLogy();
    chain.Draw("tgt_up_bang >>hbeam14(100,-0.4,0.4)",
               "triggers[5]==1", "colz");
    c1.cd(3);
    gPad->SetLogy();
    chain.Draw("tgt_dwn_aang >>hbeam15(100,-0.4,0.4)",
               "triggers[5]==1", "colz");
    c1.cd(4);
    gPad->SetLogy();
    chain.Draw("tgt_dwn_bang >>hbeam16(100,-0.4,0.4)",
               "triggers[5]==1", "colz");
    c1.Print("out/beam_tracking_dcs.pdf", "pdf");

    c1.Clear();
    c1.Divide(2,2);

    c1.cd(1);
    gPad->SetLogy();
    chain.Draw("tgt_up_theta*57.3 >>hbeam17(100,-20,20)",
               "triggers[5]==1", "colz");
    c1.cd(2);
    gPad->SetLogy();
    chain.Draw("tgt_up_phi*57.3 >>hbeam18(100,-20,20)",
               "triggers[5]==1", "colz");
    c1.cd(3);
    gPad->SetLogy();
    chain.Draw("tgt_dwn_theta*57.3 >>hbeam19(100,-20,20)",
               "triggers[5]==1", "colz");
    c1.cd(4);
    gPad->SetLogy();
    chain.Draw("tgt_dwn_phi*57.3 >>hbeam20(100,-20,20)",
               "triggers[5]==1", "colz");
    c1.Print("out/beam_tracking_dcs.pdf", "pdf");

    c1.Clear();
    c1.Divide(2,2);

    c1.cd(1);
    gPad->SetLogy();
    chain.Draw("s1dc_aang >>hbeam21(100,-0.4,0.4)",
               "triggers[5]==1", "colz");
    c1.cd(2);
    gPad->SetLogy();
    chain.Draw("s1dc_bang >>hbeam22(100,-0.4,0.4)",
               "triggers[5]==1", "colz");
    c1.cd(3);
    gPad->SetLogy();
    chain.Draw("fdc0_aang >>hbeam23(100,-0.4,0.4)",
               "triggers[5]==1", "colz");
    c1.cd(4);
    gPad->SetLogy();
    chain.Draw("fdc0_bang >>hbeam24(100,-0.4,0.4)",
               "triggers[5]==1", "colz");
    c1.Print("out/beam_tracking_dcs.pdf", "pdf");

    c1.Clear();
    c1.Divide(2,2);

    c1.cd(1);
    gPad->SetLogy();
    chain.Draw("s1dc_theta*57.3 >>hbeam22(100,-20,20)",
               "triggers[5]==1", "colz");
    c1.cd(2);
    gPad->SetLogy();
    chain.Draw("s1dc_phi*57.3 >>hbeam23(100,-20,20)",
               "triggers[5]==1", "colz");
    c1.cd(3);
    gPad->SetLogy();
    chain.Draw("fdc0_theta*57.3 >>hbeam24(100,-20,20)",
               "triggers[5]==1", "colz");
    c1.cd(4);
    gPad->SetLogy();
    chain.Draw("fdc0_phi*57.3 >>hbeam25(100,-20,20)",
               "triggers[5]==1", "colz");
    c1.Print("out/beam_tracking_dcs.pdf)", "pdf");
}
