{
    auto ds = s13::ana::MakeDataset({273,274,275,276});
    TTree *tt = ds.GetTree();

    TCut target_cut = "TargetUp_YPos*TargetUp_YPos + TargetUp_XPos*TargetUp_XPos < 36";

    tt->Draw("ESPRI_RDCL_aAng:ESPRI_RDCL_XPos >>h(400,-250,250,400,-1500,1500)",
             "TriggerBit==5" && target_cut, "colz");
    gPad->SetLogz();
}
