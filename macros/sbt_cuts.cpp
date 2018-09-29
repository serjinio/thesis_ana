{
    auto ds = s13::ana::MakeDataset({273,274,275,276});
    TTree *tt = ds.GetTree();

    TCut main_bunch_cut = "SBT1Pla_TCal > 358 && SBT1Pla_TCal < 375";

    TCanvas c1("c1");

    c1.Divide(2);
    c1.cd(1);
    tt->Draw("SBT1Pla_TCal >>hist_sbt1(400,0,500)", "", "colz");
    gPad->SetLogy();

    c1.cd(2);
    tt->Draw("SBT1Pla_TCal >>hist_sbt2(400,0,500)", main_bunch_cut, "colz");
    gPad->SetLogy();
}
