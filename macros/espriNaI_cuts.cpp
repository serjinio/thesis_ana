{

  auto ds = s13::ana::MakeDataset({273});
  TTree *tt = ds.GetTree();

  TCut espriPlaR_low_thresh_cut = "ESPRIPlaR_QCal>250";
  TCut espriPlaL_low_thresh_cut = "ESPRIPlaL_QCal>250";

  TCanvas *c1 = new TCanvas();
  c1->Divide(1);

  int num_nai = 28;
  for (int i = 0; i < num_nai; i++) {
    c1->cd(i + 1 - 7);
    TString hist_fmt, hist_title;
    hist_fmt.Form("ESPRI_NaIPla_QCal[%i] >>NaI_PMT%i(400,0,4000)", i, i+1);
    hist_title.Form("NaI PMT #%i - E", i+1);
    cout << "Drawing: " << hist_fmt << "..." << endl << flush;

    tt->Draw(hist_fmt, "TriggerBit==5", "colz");
    gPad->SetLogy();
    gPad->SetTitle(hist_title);

    if (i == 0) { // open file
      c1->Print("out/data_check1/ESPRI_NaI_E.pdf(", "pdf");
    } else if (i == num_nai - 1) {  // close file
      c1->Print("out/data_check1/ESPRI_NaI_E.pdf)", "pdf");
    } else { // print new page
      c1->Print("out/data_check1/ESPRI_NaI_E.pdf", "pdf");
    }
  }
}
