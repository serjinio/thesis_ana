{

    auto ds = s13::ana::MakeDataset({273});
    TTree *tt = ds.GetTree();

    TCanvas *c1 = new TCanvas();
    c1->Divide(1);

    tt->Draw("TargetUp_XPos - TargetDwn_XPos >>h(200,-80,80)", "TriggerBit==5", "colz");
    c1->Print("out/data_check1/TargetPos_deltas.pdf(", "pdf");

    tt->Draw("TargetUp_YPos - TargetDwn_YPos >>h(200,-80,80)", "TriggerBit==5", "colz");
    c1->Print("out/data_check1/TargetPos_deltas.pdf)", "pdf");
}
