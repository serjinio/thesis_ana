#include <iostream>
#include "TChain.h"
#include <string.h>

#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)


TChain g_chain_up{"scattree"};
TChain g_chain_down{"scattree"};
TChain g_chain_total{"scattree"};
TChain g_chain_carbon{"scattree"};

TChain* init_dataset() {
  for (int i = 153; i < 154; i++) {
    TString name;
    name.Form("scattree_rootfiles/scattree_run%i.root", i);
    std::cout << "Adding file to chain: " << name << std::endl;
    g_chain_total.Add(name);
  }

  return &g_chain_total;
}

TH1* hdraw(TTree& tree, TString name, TString draw_cmd,
           TString binning, TCut cuts, TString title,
           TString xaxis_title, TString yaxis_title,
           TString draw_opts = "") {
  TString hstr = TString::Format("%s >>%s%s", draw_cmd.Data(), name.Data(),
                                 binning.Data());
  tree.Draw(hstr, cuts, draw_opts);
  TH1 *hist = static_cast<TH1*>(gPad->GetPrimitive(name));
  assert(hist != nullptr && "Histogram plotting failed!");
  if (yaxis_title != "") {
    hist->GetYaxis()->SetTitle(yaxis_title);
  }
  if (xaxis_title != "") {
    hist->GetXaxis()->SetTitle(xaxis_title);
  }
  if (title != "") {
    hist->SetTitle(title);
  }

  return hist;
}

void inline print_pdf_first(TCanvas& canv) {
  TString output_filename = TString::Format("out/%s.pdf(", __FILENAME__);
  canv.Print(output_filename, "pdf");
}

void inline print_pdf(TCanvas& canv) {
  TString output_filename = TString::Format("out/%s.pdf", __FILENAME__);
  canv.Print(output_filename, "pdf");
}

void inline print_pdf_last(TCanvas& canv) {
  TString output_filename = TString::Format("out/%s.pdf)", __FILENAME__);
  canv.Print(output_filename, "pdf");
}

void draw_beam_profile(TCanvas &c1) {
  TCut vertexXY_cut =
    TString::Format("triggers[1] == 1").Data();
  TString title = TString::Format("Vertex X & Y at pol. target");

  c1.Divide(1, 1);
  c1.cd(1);
  hdraw(g_chain_total, "tgtXY_pol_target",
        "tgt_up_ypos:tgt_up_xpos",
        "(200,-40,40, 200, -40, 40)",
        vertexXY_cut,
        title,
        "X [mm]",
        "Y [mm]",
        "col");
  gPad->SetLogz();
  print_pdf_first(c1);

  c1.Clear();
  c1.Divide(1, 1);
  c1.cd(1);
  hdraw(g_chain_total, "tgtX_pol_target",
        "tgt_up_xpos",
        "(60,-40,40)",
        vertexXY_cut,
        title,
        "X [mm]",
        "Counts");
  gPad->SetLogz();
  print_pdf(c1);

  c1.Clear();
  c1.Divide(1, 1);
  c1.cd(1);
  hdraw(g_chain_total, "tgtY_pol_target",
        "tgt_up_ypos",
        "(60,-40,40)",
        vertexXY_cut,
        title,
        "Y [mm]",
        "Counts");
  gPad->SetLogz();
  print_pdf_last(c1);
}

void count_incident_particles(TCanvas &c1, Double_t vertexXY_radius) {
  TCut vertexXY_cut =
    TString::Format("triggers[1] == 1 && pow(tgt_up_xpos + 2, 2) "
                    "+ pow(tgt_up_ypos + 2, 2) < %.2f",
                    std::pow(vertexXY_radius, 2)).Data();
  TString title = TString::Format("Vertex X&Y pol. target (R<%.2f)",
                                  vertexXY_radius);

  c1.Clear();
  c1.Divide(2, 2);
  c1.cd(1);

  hdraw(g_chain_total, "tgtXY_pol_target",
        "tgt_up_ypos:tgt_up_xpos",
        "(100,-40,40, 100, -40, 40)",
        vertexXY_cut,
        title,
        "X [mm]",
        "Y [mm]",
        "col");
  gPad->SetLogz();
  gPad->SetLineWidth(3);

  c1.cd(2);
  hdraw(g_chain_total, "tgtX_pol_target",
        "tgt_up_xpos",
        "(100,-40,40)",
        vertexXY_cut,
        title,
        "X [mm]",
        "Counts");
  gPad->SetLogz();
  gPad->SetLineWidth(3);

  c1.cd(3);
  gPad->SetLineWidth(3);
  hdraw(g_chain_total, "tgtY_pol_target",
        "tgt_up_ypos",
        "(100,-40,40)",
        vertexXY_cut,
        title,
        "X [mm]",
        "Counts");
  gPad->SetLogz();
}

void beam_intens_eval() {
  gROOT->SetStyle("Pub");
  gStyle->SetLabelFont(2, "xyz");
  gStyle->SetTitleFont(2, "xyz");
  gStyle->SetLegendFont(3);
  gStyle->SetLabelSize(0.06, "xyz");
  gStyle->SetTitleSize(0.06, "xyz");
  gStyle->SetPadLeftMargin(0.27);
  gStyle->SetPadBottomMargin(0.20);
  gStyle->SetTitleOffset(1.9, "y");
  gStyle->SetLineWidth(2);
  gStyle->SetOptStat(1111111);
  init_dataset();
  // init_carbon_dataset();

  TString output_root_filename = TString::Format("out/%s.root", __FILENAME__);
  TFile hist_out(output_root_filename, "RECREATE");
  TCanvas c1("c1");

  // count_incident_particles(c1, 4);
  // print_pdf_first(c1);

  // count_incident_particles(c1, 6);
  // print_pdf(c1);

  // count_incident_particles(c1, 8);
  // print_pdf(c1);

  // count_incident_particles(c1, 10);
  // print_pdf(c1);

  // count_incident_particles(c1, 120);
  // print_pdf_last(c1);

  draw_beam_profile(c1);
 }
