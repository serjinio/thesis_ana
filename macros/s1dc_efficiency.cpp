#include <iostream>

#include "TChain.h"


void s1dc_efficiency() {
    TChain chain{"scattree"};

    for (int i = 272; i < 274; i++){
        TString name;
        name.Form("scattree_rootfiles/scattree_run%i.root", i);
        std::cout << "Adding file to chain: " << name << std::endl;
        chain.Add(name);
    }

    chain.Print();

    TFile f_hist_out("out/s1dc_efficiency.root", "RECREATE");
    TCanvas c1("c1");
    gStyle->SetOptStat(1111111);

    c1.Divide(2);

    c1.cd(1);
    //gPad->SetLogz();
    chain.Draw("s1dc_xpos >>s1dc_eff1(400,-80,80)",
               "triggers[5]==1 && s1dc_ypos < 4 && s1dc_ypos > -4", "colz");
    c1.cd(2);
    chain.Draw("fdc0_xpos >>s1dc_eff2(400,-80,80)",
               "triggers[5]==1 && fdc0_ypos < 2 && fdc0_ypos > -2", "colz");
    c1.Print("out/s1dc_efficiency.pdf(", "pdf");
    //f_hist_out.Write();

    c1.Clear();
    c1.Divide(2);

    c1.cd(1);
    chain.Draw("s1dc_ypos:s1dc_xpos >>s1dc_eff10(400,-80,80,400,-80,80)", "",
               "colz");
    c1.cd(2);
    chain.Draw("fdc0_ypos:fdc0_xpos >>s1dc_eff11(400,-80,80,400,-80,80)", "",
               "colz");
    c1.Print("out/s1dc_efficiency.pdf)", "pdf");

    f_hist_out.Write();
    f_hist_out.Close();
}
