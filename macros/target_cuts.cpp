{

    auto ds = s13::ana::MakeDataset({273});
    TTree *tt = ds.GetTree();

    TCut main_bunch_cut = "SBT1Pla_TCal > 360 && SBT1Pla_TCal < 370";
    TCut target_cut = "TargetUp_YPos*TargetUp_YPos + TargetUp_XPos*TargetUp_XPos < 36";

    // angle vs. position correlation cut for RDC
    // it should select only events on RDC originating from the target region
    // TCutG* proton0_rdc_gcut = (TCutG*)
    //     (new TFile("cuts/proton_epla0_tof0.root"))->Get("CUTG");
    // proton0_rdc_gcut->SetName("proton0_rdc_gcut");
    // TCutG* proton1_rdc_gcut = (TCutG*)
    //     (new TFile("cuts/proton_epla1_tof1.root"))->Get("CUTG");
    // proton1_rdc_gcut->SetName("proton1_rdc_gcut");

    // proton & fragment angle correlation cut
    TCut p_6He_cut = "theta[1]*180/3.14+s1v[0]*0.25 > 80 && theta[1]*180/3.14+s1v[0]*0.25 < 87";
    // TCutG* phi1_phi2_corr_gcut = (TCutG*)
    //     (new TFile("cuts/phi1_phi2_corr.root"))->Get("CUTG");
    // phi1_phi2_corr_gcut->SetName("phi1_phi2_corr_gcut");

    TCanvas c1("c1");
    c1.Divide(2,1);

    c1.cd(1);
    c1.SetLogz();
    tt->Draw("TargetUp_YPos:TargetUp_XPos >>h1(400,-40,40,400,-40,40)", "", "colz");
    gPad->SetLogz();
    h1->SetTitle("Target image, no cuts");

    c1.cd(2);
    c1.SetLogz();
    tt->Draw("TargetUp_YPos:TargetUp_XPos >>h2(400,-40,40,400,-40,40)",
             main_bunch_cut && target_cut, "colz");
    gPad->SetLogz();
    h2->SetTitle("Target image, cut: R < 6 mm");

}
