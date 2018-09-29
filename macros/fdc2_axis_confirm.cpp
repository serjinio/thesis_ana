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


void fdc2_axis_confirm() {

  auto ds = s13::ana::MakeDataset({351});
  TTree *tt = ds.GetTree();


  TCanvas c1("c1");
  gStyle->SetOptStat(1111111);
  
  c1.Clear();
  c1.Divide(2);
  c1.cd(1);
  hdraw(*tt, "fdc2y_hodfu", "FDC2_YPos:HodoPlaU_Qraw", "(200,500,3000,200,-400,400)",
        "", "HODF up PMTs - FDC2 Ypos correlation", "HODF PMTs Q [ch]", "FDC2 Ypos [mm]");
  gPad->SetLogz();
  c1.cd(2);
  hdraw(*tt, "fdc2y_hodfd", "FDC2_YPos:HodoPlaD_Qraw", "(200,500,3000,200,-400,400)",
        "", "HODF down PMTs - FDC2 Ypos correlation", "HODF PMTs Q [ch]", "FDC2 Ypos [mm]");
  gPad->SetLogz();
  c1.Print("out/fdc2_axis_confirm.pdf", "pdf");

}
