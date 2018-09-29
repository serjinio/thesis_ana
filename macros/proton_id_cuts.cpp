{

    auto ds = s13::ana::MakeDataset({273,274,275,276,277,278,279,280});
    TTree *tt = ds.GetTree();

    TCut main_bunch_cut = "SBT1Pla_TCal > 340 && SBT1Pla_TCal < 380";
    TCut target_cut = "TargetUp_YPos*TargetUp_YPos + TargetUp_XPos*TargetUp_XPos < 36";

    TCut hodo_cut = "HodoPla_QCal>3.0";
    TCutG* hodo_gcut = (TCutG*)
        (new TFile("cuts/hodf_4he_cut.root"))->Get("CUTG");
    hodo_gcut->SetName("hodo_gcut");
    TCut he4_left_scat_hodo_cut = "(HodoPla_QCal[1]>13||HodoPla_QCal[2]>13||HodoPla_QCal[3]>13||HodoPla_QCal[4]>13||HodoPla_QCal[5]>13||HodoPla_QCal[6]>13||HodoPla_QCal[7]>13||HodoPla_QCal[8]>13||HodoPla_QCal[9]>13)";
    //not sure needs checking
    TCut he4_right_scat_hodo_cut = "HodoPla_QCal[11]>13||HodoPla_QCal[12]>13||HodoPla_QCal[13]>13";

    TCut espriPlaR_low_thresh_cut = "ESPRIPlaR_QCal>250";
    TCut espriPlaL_low_thresh_cut = "ESPRIPlaL_QCal>250";

    TCutG* espriPlaR_proton_gcut = (TCutG*)
        (new TFile("cuts/ESPRIPlaR_proton_cut.root"))->Get("CUTG");
    espriPlaR_proton_gcut->SetName("espriPlaR_proton_gcut");
    TCutG* espriPlaL_proton_gcut = (TCutG*)
        (new TFile("cuts/ESPRIPlaL_proton_cut.root"))->Get("CUTG");
    espriPlaL_proton_gcut->SetName("espriPlaL_proton_gcut");


    TCanvas c1("c1");
    c1.Divide(2,2);

    c1.cd(1);
    tt->Draw("ESPRIPlaR_QCal:ESPRIPlaR_TCal >>h3(600,6000,10000,600,0,4500)",
             "", "colz");
    gPad->SetLogz();

    c1.cd(2);
    tt->Draw("ESPRIPlaR_QCal:ESPRIPlaR_TCal >>h4(600,6000,10000,600,0,4500)",
             main_bunch_cut && target_cut && "TriggerBit==5" &&
             //"espriPlaR_proton_gcut" &&
             "hodo_gcut" &&
             espriPlaR_low_thresh_cut,
             "colz");
    gPad->SetLogz();

    c1.cd(3);
    tt->Draw("ESPRIPlaL_QCal:ESPRIPlaL_TCal >>h5(600,6000,10000,600,0,4500)",
             "", "colz");
    gPad->SetLogz();

    c1.cd(4);
    tt->Draw("ESPRIPlaL_QCal:ESPRIPlaL_TCal >>h6(600,6000,10000,600,0,4500)",
             main_bunch_cut && target_cut && "TriggerBit==5" &&
             //"espriPlaL_proton_gcut" &&
             "hodo_gcut" &&
             espriPlaL_low_thresh_cut,
             "colz");
    gPad->SetLogz();
}
