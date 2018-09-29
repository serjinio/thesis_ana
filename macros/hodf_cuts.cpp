{
    auto ds = s13::ana::MakeDataset({273,274,275,276});
    TTree *tt = ds.GetTree();

    TCut main_bunch_cut = "F7Pla_TCal > 670 && F7Pla_TCal < 700";

    // hodoscope cuts
    TCut hodo_cut = "HodoPla_QCal>0.2";
    TCutG* hodo_gcut = (TCutG*)
        (new TFile("cuts/hodf_4he_cut.root"))->Get("CUTG");
    hodo_gcut->SetName("hodo_gcut");
    TCut he4_left_scat_hodo_cut = "(HodoPla_QCal[1]>13||HodoPla_QCal[2]>13||HodoPla_QCal[3]>13||HodoPla_QCal[4]>13||HodoPla_QCal[5]>13||HodoPla_QCal[6]>13||HodoPla_QCal[7]>13||HodoPla_QCal[8]>13||HodoPla_QCal[9]>13)";
// maybe for right scattering, not sure
    TCut he4_right_scat_hodo_cut = "HodoPla_QCal[11]>13||HodoPla_QCal[12]>13||HodoPla_QCal[13]>13";
// ?
    TCut espri_min_e = "enai>100";

    TCanvas c1("c1");

    // c1.Divide(5,5);
    // for (int i = 0; i < 24; i++) {
    //     c1.cd(i + 1);
    //     tt->Draw(TString::Format("HodoPlaU_Qraw[%i] >> hist_hodo%i(80, 0, 30)",
    //                            i, i));
    // }

    c1.Divide(2);
    c1.cd(1);
    tt->Draw("HodoPla_QCal:HodoID >>hist_hodo50(24,0,24,400,0,45)",
             main_bunch_cut, "colz");
    gPad->SetLogz();

    c1.cd(2);
    tt->Draw("HodoPla_QCal:HodoPla_TCal >>hist_hodo51(400,-300,-50,400,0,45)",
             main_bunch_cut, "colz");
    gPad->SetLogz();
}
