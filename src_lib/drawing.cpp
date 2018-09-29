
#include "drawing.hpp"


////////////////////////////////////////////////////////////
// histogramming common functions
////////////////////////////////////////////////////////////


void
s13::ana::draw_marker_line(Double_t x1, Double_t y1, Double_t x2, Double_t y2,
                           Double_t line_width) {
  TLine *line = new TLine(x1, y1, x2, y2);
  line->SetLineColor(kRed);
  line->SetLineWidth(line_width);
  line->Draw();
}

void
s13::ana::draw_marker_line_ysym(Double_t abs_x, Double_t y1, Double_t y2,
                                Double_t line_width) {
  draw_marker_line(-abs_x, y1, -abs_x, y2, line_width);
  draw_marker_line(abs_x, y1, abs_x, y2, line_width);
}

void
s13::ana::draw_marker_circle(Double_t center_x, Double_t center_y, Double_t radius,
                             Double_t width) {
  TEllipse* ellipse = new TEllipse(center_x, center_y, radius, radius);
  ellipse->SetLineColor(kRed);
  ellipse->SetLineWidth(width);
  ellipse->SetFillStyle(0);
  ellipse->SetFillColor(0);
  ellipse->Draw();
}

TH1* s13::ana::hdraw(TTree& tree, TString name, TString draw_cmd,
                     TString binning, TCut cuts, TString title,
                     TString xaxis_title, TString yaxis_title,
                     TString draw_opts) {
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

TH1*
s13::ana::draw_theta_corr(TChain& chain, const TCut& cuts,
                          TString title, TString draw_opts, bool draw2d,
                          int bin_num, TString name) {
  static int hist_num = 1;
  TCut p_angle_acceptance_cut = "p_theta_eff*57.3>55 && p_theta_eff*57.3<70";
  TString hist_name = name == "" ? TString::Format("thcorr_%i", hist_num) : name;
  if (title == "") {
    title = hist_name;
  }

  if (draw2d) {
    TString limits = TString::Format("(%i,0,15,%i,55,70)", bin_num, bin_num);
    hdraw(chain, hist_name,
          "p_theta_eff*57.3:s1dc_theta*57.3",
          limits,
          cuts,
          title,
          "#theta_{6He} [lab. deg.]",
          "#theta_p [lab. deg.]", draw_opts);
  } else {
    if (name == "") {
      hist_name += "_1d";
    }
    TString limits = TString::Format("(%i,-8,8)", bin_num);
    hdraw(chain, hist_name,
          "s1dc_theta*57.3 - he_theta_theor*57.3",
          limits,
          cuts && p_angle_acceptance_cut,
          title + " 1D",
          "Fragment: #theta_{exp} - #theta_{theor} [lab. deg.]",
          "Counts", draw_opts);
  }

  hist_num += 1;
  TH1 *hist = static_cast<TH1*>(gPad->GetPrimitive(hist_name));

  return hist;
}

TH1* s13::ana::draw_theta_corr_2d(TChain& chain, const TCut& cuts,
                                  TString title, TString draw_opts,
                                  int bin_num, TString name) {
  return s13::ana::draw_theta_corr(chain, cuts, title, draw_opts, true, bin_num, name);
}


TH1* s13::ana::draw_theta_corr_1d(TChain& chain, const TCut& cuts,
                                  TString title, TString draw_opts,
                                  int bin_num, TString name) {
  return s13::ana::draw_theta_corr(chain, cuts, title, draw_opts, false, bin_num, name);
}

/*
  Function build a TGraph for a given scattering yield object
*/
TGraph*
s13::ana::BuildAyGraph(ScatteringYield yield, TString title) {
  Double_t* xs = new Double_t[yield.Size()];
  Double_t* ys = new Double_t[yield.Size()];
  std::copy(yield.ArgsRangeStartIter(), yield.ArgsRangeEndIter(), xs);
  std::copy(yield.ValsRangeStartIter(), yield.ValsRangeEndIter(), ys);

  TGraph* graph = new TGraph(yield.Size(), xs, ys);
  graph->SetLineColor(2);
  graph->SetLineWidth(4);
  graph->SetMarkerColor(4);
  graph->SetMarkerStyle(21);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitle("Lab angle [deg]");
  graph->GetYaxis()->SetTitle("Asymmetry");
  //      graph->SetMinimum(-1.2);
  //      graph->SetMaximum(1.2);

  return graph;
};

/*
  Function build a TGraph for a given scattering yield object
*/
TGraph*
s13::ana::BuildAsymmetryGraphErrors(ScatteringYield yield,
                                    TString title) {
  Double_t* xs = new Double_t[yield.Size()];
  Double_t* ys = new Double_t[yield.Size()];
  Double_t* yerr = new Double_t[yield.Size()];
  std::copy(yield.ArgsRangeStartIter(), yield.ArgsRangeEndIter(), xs);
  std::copy(yield.ValsRangeStartIter(), yield.ValsRangeEndIter(), ys);
  std::copy(yield.ErrsRangeStartIter(), yield.ErrsRangeEndIter(), yerr);

  TGraphErrors* graph = new TGraphErrors(yield.Size(), xs, ys, nullptr, yerr);
  graph->SetLineColor(kGreen);
  graph->SetLineWidth(2);
  graph->SetMarkerColor(kRed);
  graph->SetMarkerStyle(21);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitle("Lab angle [deg]");
  graph->GetYaxis()->SetTitle("Asymmetry");
  // graph->SetMinimum(-1.2);
  // graph->SetMaximum(1.2);

  return graph;
};

/*
  Function build a TGraph for a given scattering yield object
*/
TGraph*
s13::ana::BuildTGraphErrors(ScatteringYield yield,
                            TString title,
                            TString xaxis_title,
                            TString yaxis_title) {
  Double_t* xs = new Double_t[yield.Size()];
  Double_t* ys = new Double_t[yield.Size()];
  Double_t* yerr = new Double_t[yield.Size()];
  std::copy(yield.ArgsRangeStartIter(), yield.ArgsRangeEndIter(), xs);
  std::copy(yield.ValsRangeStartIter(), yield.ValsRangeEndIter(), ys);
  std::copy(yield.ErrsRangeStartIter(), yield.ErrsRangeEndIter(), yerr);

  TGraphErrors* graph = new TGraphErrors(yield.Size(), xs, ys, nullptr, yerr);
  graph->SetLineColor(kGreen);
  graph->SetLineWidth(2);
  graph->SetMarkerColor(kRed);
  graph->SetMarkerStyle(22);
  graph->SetMarkerSize(1);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitle(xaxis_title);
  graph->GetYaxis()->SetTitle(yaxis_title);

  return graph;
}

TGraph*
s13::ana::BuildCsGraphFromData(const s13::io::CsvColumnsVector& csv_data,
                               std::string title) {
  s13::misc::MessageLogger log("BuildCsGraphFromData()");
  static std::vector<int> line_colors =
    {kBlue, kBlack, kMagenta, kRed, kCyan, kGreen};
  static size_t line_color_idx = 0;
  static size_t line_style = 2;

  int rows_num = csv_data.at(0).size();
  double xes[rows_num], ys[rows_num];
  // log.debug("Loading data for %s...", title.c_str());
  for (int i = 0; i < rows_num; i++) {
    xes[i] = csv_data.at(0).at(i);
    ys[i] = csv_data.at(1).at(i);
    // log.debug("\tthetacm: %.4f; CS: %.4f", xes[i], ys[i]);
  }

  TGraph* gr = new TGraph(rows_num, xes, ys);
  gr->SetLineWidth(1);
  gr->SetLineColor(line_colors.at(line_color_idx));
  gr->SetLineStyle(line_style++);
  // gr->SetMarkerColor(line_colors.at(line_colors.size() - 1 - line_color_idx));
  gr->SetTitle(title.c_str());
  gr->GetXaxis()->SetTitle("CM angle [deg]");
  gr->GetYaxis()->SetTitle("CS [mbarns]");

  ++line_color_idx;
  if (line_color_idx == line_colors.size()) {
    line_color_idx = 0;
  }

  return gr;
}

/*
  Function returns a graph of Ay from Moss paper
*/
TGraph* s13::ana::Build4HeAyGraph() {
  static const double p_lab_theta[] = {55., 57., 59., 61., 64, 67, 69};
  static const double ay[] = {-0.9, -1.0, -0.97, -0.86, -0.65, -0.32, -0.09};
  static const int num_points = 7;

  TGraph* graph = new TGraph(num_points, p_lab_theta, ay);
  graph->SetLineColor(3);
  graph->SetLineWidth(4);
  graph->SetMarkerColor(5);
  graph->SetMarkerStyle(22);
  graph->SetTitle("4He A_y (Moss)");
  graph->GetXaxis()->SetTitle("Lab angle [deg]");
  graph->GetYaxis()->SetTitle("A_y");
  graph->SetMinimum(-1.2);
  graph->SetMaximum(1.2);

  return graph;
}

/*
  Function prepares multi-graph object for asymmetry / Ay plotting
*/
TMultiGraph*
s13::ana::PrepareAyMultiGraph(TString title,
                              std::vector<TGraph*> graphs) {
  TMultiGraph* mg = new TMultiGraph();
  mg->SetMinimum(-1.2);
  mg->SetMaximum(1.2);
  mg->SetTitle(TString::Format("%s;%s;%s", title.Data(),
                               "Lab angle [deg.]", "Asymmetry / A_y"));
  for (TGraph* gr : graphs) {
    mg->Add(gr);
  }

  return mg;
}

void
s13::ana::AddHe6TheorCses(TMultiGraph* mg) {
  std::vector<s13::io::CsvColumnsVector> datas;
  datas.push_back(s13::io::load_csv("data/cs_karataglidis.csv", 2, 2));
  // datas.push_back(s13::io::load_csv("data/cs_nh_karataglidis.csv", 2, 2));
  // datas.push_back(s13::io::load_csv("data/cs_toyokawa_g-matrix.csv", 2, 2));
  // datas.push_back(s13::io::load_csv("data/cs_Elster.csv", 2, 2));
  // datas.push_back(s13::io::load_csv("data/cs_weppner_elsther.csv", 2, 2));
  // datas.push_back(s13::io::load_csv("data/cs_Kaki1.csv", 2, 2));
  datas.push_back(s13::io::load_csv("data/cs_kaki_paper_tmav1.csv", 2, 2));
  // datas.push_back(s13::io::load_csv("data/cs_Kaki2.csv", 2, 2));
  datas.push_back(s13::io::load_csv("data/cs_kaki_paper_tmav2.csv", 2, 2));
  // datas.push_back(s13::io::load_csv("data/cs_Kaki3.csv", 2, 2));
  datas.push_back(s13::io::load_csv("data/cs_kaki_paper_tmav3.csv", 2, 2));
  // datas.push_back(s13::io::load_csv("data/he6_theoretical_cs_s13_proposal.csv", 2, 1));

  std::vector<std::string> theor_cses_titles =
    {
      "Karataglidis",
      // "Karataglidis NH",
      // "Toyokawa g-matrix",
      // "Elsther",
      // "Weppner, Elsther"
      "Kaki 1",
      "Kaki 2",
      "Kaki 3",
      // "S13 prop."
    };

  int i = 0;
  for (auto& el : datas) {
    TGraph* graph = s13::ana::BuildCsGraphFromData(el, theor_cses_titles.at(i));
    mg->Add(graph);
    i++;
  }
}

TH1* s13::ana::init_ecut_bin_hist(double bin_min_theta, double bin_max_theta,
                                  double xlim) {
  auto th1 = &s13::ana::make_th1;
  auto tit = s13::ana::tstrfmt("Rel. E spectrum (%.2f<p_{#theta}<%.2f)",
                               bin_min_theta, bin_max_theta);
  auto hist = th1(120, -xlim, xlim, tit, "Rel. E [%/100]", "Counts");
  return hist;
}

TH1* s13::ana::init_phi_bin_hist(double bin_min_theta, double bin_max_theta,
                                 double xlim) {
  auto th1 = &s13::ana::make_th1;
  auto tit = s13::ana::tstrfmt("#Delta #phi spectrum (%.2f<p_{#theta}<%.2f)",
                               bin_min_theta, bin_max_theta);
  auto hist = th1(120, 180-xlim, 180+xlim, tit, "#Delta #phi [deg]", "Counts");
  return hist;
}

TH1* s13::ana::init_vertR_bin_hist(double bin_min_theta, double bin_max_theta,
                                 double rlim) {
  auto th1 = &s13::ana::make_th1;
  auto tit = s13::ana::tstrfmt("Vertex R spectrum (%.2f<p_{#theta}<%.2f)",
                               bin_min_theta, bin_max_theta);
  auto hist = th1(120, 0, rlim, tit, "Vertex R [deg]", "Counts");
  return hist;
}
