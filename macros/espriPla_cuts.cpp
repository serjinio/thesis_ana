{

    auto ds = s13::ana::MakeDataset({273,274,275,276});
    TTree *tt = ds.GetTree();

    TCut espriPlaR_low_thresh_cut = "ESPRIPlaR_QCal>250";
    TCut espriPlaL_low_thresh_cut = "ESPRIPlaL_QCal>250";

    TCanvas c1("c1");
    c1.Divide(2,2);

    c1.cd(1);
    tt->Draw("ESPRIPlaR_QCal >>h1(400,0,4000)", "", "colz");
    gPad->SetLogy();

    c1.cd(2);
    tt->Draw("ESPRIPlaR_QCal >>h2(400,0,4000)", espriPlaR_low_thresh_cut, "colz");
    gPad->SetLogy();

    c1.cd(3);
    tt->Draw("ESPRIPlaL_QCal >>h3(400,0,4000)", "", "colz");
    gPad->SetLogy();

    c1.cd(4);
    tt->Draw("ESPRIPlaL_QCal >>h4(400,0,4000)", espriPlaL_low_thresh_cut, "colz");
    gPad->SetLogy();

}
