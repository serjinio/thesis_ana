#include <iostream>

#include "TChain.h"


void proton_e_de() {
    TChain chain{"scattree"};

    for (int i = 272; i < 274; i++){
        TString name;
        name.Form("scattree_rootfiles/scattree_run%i.root", i);
        std::cout << "Adding file to chain: " << name << std::endl;
        chain.Add(name);
    }

    TCanvas c1("c1");
    c1.Divide(2,2);

    c1.cd(1);
    gPad->SetLogy();
    chain.Draw("esl_xpos >> hrdc1(200, -250, 250)", "triggers[5]==1", "colz");
    c1.cd(2);
    gPad->SetLogy();
    chain.Draw("esl_ypos >> hrdc2(200, -250, 250)", "triggers[5]==1", "colz");
    c1.cd(3);
    gPad->SetLogy();
    chain.Draw("esr_xpos >> hrdc3(200, -250, 250)", "triggers[5]==1", "colz");
    c1.cd(4);
    gPad->SetLogy();
    chain.Draw("esr_ypos >> hrdc4(200, -250, 250)", "triggers[5]==1", "colz");
    c1.Print("out/rdc_data.pdf(", "pdf");

    c1.Clear();
    c1.Divide(2,4);

    c1.cd(1);
    gPad->SetLogy();
    chain.Draw("esl_aang >> hrdc1(200, -250, 250a)", "triggers[5]==1", "colz");

}
