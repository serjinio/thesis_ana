/**
 *   \file common.cpp
 *   \brief Sources of common functions
 *
 */

#include "TROOT.h"
#include "TPad.h"
#include "TStyle.h"
#include "TF1.h"
#include "TPaveStats.h"

#include "common.hpp"


std::string s13::ana::espri_selector_to_str(s13::ana::Espri espri_sel) {
  switch (espri_sel) {
  case s13::ana::Espri::left:
    return "left";
  case s13::ana::Espri::right:
    return "right";
  case s13::ana::Espri::both:
    return "both";
  default:
    throw std::invalid_argument("Invalid s13::ana::Espri enum value!");
  }
}

std::string
s13::ana::ds_type_to_str(DatasetType ds_type) {
  switch (ds_type) {
  case DatasetType::he6:
    return "he6";
  case DatasetType::he4:
    return "he4";
  case DatasetType::carbon:
    return "carbon";
  default:
    throw std::invalid_argument("Invalid s13::ana::DatasetType enum value!");
  }
}

std::string
s13::ana::pol_dir_to_str(s13::ana::PolarizationDirection pol_dir) {
  switch (pol_dir) {
  case PolarizationDirection::up:
    return "pol_up";
  case PolarizationDirection::down:
    return "pol_down";
  case PolarizationDirection::both:
    return "pol_both";
  default:
    throw std::invalid_argument("Invalid s13::ana::PolarizationDirection enum value!");
  }
}

TH1* s13::ana::make_th1(int bins_no, double min_x, double max_x,
                        TString title,
                        TString xlabel, TString ylabel) {
  static int th1_hist_num = 0;
  TString name = s13::ana::tstrfmt("th1_%d", th1_hist_num);
  TH1* retval = new TH1F(name, title, bins_no, min_x, max_x);
  retval->GetXaxis()->SetTitle(xlabel);
  retval->GetYaxis()->SetTitle(ylabel);
  th1_hist_num += 1;
  return retval;
}

TH2* s13::ana::make_th2(int nbins_x, double min_x, double max_x,
                        int nbins_y, double min_y, double max_y,
                        TString title,
                        TString xlabel, TString ylabel) {
  static int th2_hist_num = 0;
  TString name = s13::ana::tstrfmt("th2_%d", th2_hist_num);
  TH2* retval = new TH2F(name, title,
                         nbins_x, min_x, max_x,
                         nbins_y, min_y, max_y);
  retval->GetXaxis()->SetTitle(xlabel);
  retval->GetYaxis()->SetTitle(ylabel);
  th2_hist_num += 1;
  return retval;
}

THStack* s13::ana::make_thstack(std::vector<TH1*> hists, TString title) {
  assert(hists.size() != 0);
  static int hstack_id = 0;
  THStack *hs = new THStack(s13::ana::tstrfmt("hs_%i", hstack_id), title);

  for (auto el : hists) {
    hs->Add(el);
  }

  ++hstack_id;
  return hs;
}

THStack* s13::ana::draw_thstack(std::vector<TH1*> hists, TString title, TString draw_opt) {
  static std::vector<int>
    colors = {kBlue, kRed, kGreen, kCyan, kOrange, kPink, kBlack, kMagenta};
  assert(hists.size() <= colors.size());

  for (size_t i = 0; i < hists.size(); ++i) {
    hists.at(i)->SetLineColor(colors.at(i));
  }

  auto hstack = make_thstack(hists, title);
  hstack->Draw(draw_opt);
  hstack->GetXaxis()->SetTitle(hists.at(0)->GetXaxis()->GetTitle());
  hstack->GetYaxis()->SetTitle(hists.at(0)->GetYaxis()->GetTitle());

  return hstack;
}

TPaveText* s13::ana::draw_text(TString text, double x1, double y1, double x2, double y2,
                               TString draw_opts) {
  TPaveText *pt = new TPaveText(x1, y1, x2, y2, "NDC");
  pt->AddText(text);
  pt->Draw(draw_opts);
  return pt;
}

std::istream& s13::ana::operator>>(std::istream& str, CSVRow& columns) {
  columns.readNextRow(str);
  return str;
}

std::ostream& s13::ana::operator<<(std::ostream& os, const FitParams& obj)
{
  // write params
  for (auto iter = obj.params.begin(); iter != obj.params.end(); iter++) {
    os << *iter << ",";
  }
  // write errors
  for (auto iter = obj.errors.begin(); iter != obj.errors.end(); iter++) {
    os << *iter << ",";
  }
  // object writing finished, end the line
  os << std::endl;

  return os;
}

std::istream& s13::ana::operator>>(std::istream& is, FitParams& obj)
{
  s13::misc::MessageLogger logger("[FitParams] operator>>()");

  std::string line;
  std::getline(is, line);
  logger.debug("Got a line from input stream: %s", line.c_str());

  std::stringstream lineStream(line);
  std::string cell;
  obj.params.clear();
  obj.errors.clear();

  std::vector<std::string> columns;
  while(std::getline(lineStream, cell, ',')) {
    columns.push_back(cell);
  }

  if (columns.size() % 2 != 0) {
    is.setstate(std::ios::failbit);
    std::cerr << "ERROR: Invalid size of input columns for FitParams deserialization." <<
      std::endl;
  }

  obj.SetParamsNum(columns.size() / 2);
  logger.debug("Params number set to %d", columns.size() / 2);
  for (size_t i = 0; i < columns.size(); i++) {
    Double_t data_value = std::stod(columns[i]);
    if (i + 1 <= columns.size() / 2) {
      obj.params.push_back(data_value);
    } else {
      obj.errors.push_back(data_value);
    }
  }

  return is;
}

s13::io::CsvColumnsVector
s13::io::load_csv(std::string filename, size_t columns_number,
                  size_t skip_rows) {
  s13::misc::MessageLogger log("load_csv()");
  std::ifstream ifile(filename);
  s13::io::CsvColumnsVector columns;

  for (size_t i = 0; i < columns_number; i++) {
    std::vector<double> vec;
    columns.push_back(vec);
  }

  log.debug("Loading data from %s...", filename.c_str());
  int rows_counter = 0;
  for (s13::ana::CSVIterator iter(ifile, skip_rows);
       iter != s13::ana::CSVIterator(); iter++) {
    auto csv_row = *iter;
    assert(csv_row.size() == columns_number &&
           "Wrong column number in input file!");

    for (size_t i = 0; i < columns_number; i++) {
      double val = std::stod(csv_row[i]);
      // log.debug("setting col#%i to %.6f", i, val);
      columns.at(i).push_back(val);
    }

    ++rows_counter;
  }
  log.debug("Total number of rows loaded from file: %d", rows_counter);

  return columns;
}

/*
  Fits passed 1D hist with a gaussian.
  Returns: FitParams
*/
s13::ana::FitParams s13::ana::fit_hist(TH1* hist) {
  s13::ana::FitParams fit_params{3};
  const Double_t *fit_params_raw;
  const Double_t *fit_errors_raw;

  if (hist->GetEntries() == 0) {
    std::cout << "WARN: Empty histogram passed for fitting!" << std::endl;
    const Double_t zero_params[] = {0,0,0};
    fit_params.SetParams(zero_params);
    fit_params.SetErrors(zero_params);
    return fit_params;
  }

  hist->Fit("gaus", "W");
  fit_params_raw = hist->GetFunction("gaus")->GetParameters();
  fit_errors_raw = hist->GetFunction("gaus")->GetParErrors();

  fit_params.SetParams(fit_params_raw);
  fit_params.SetErrors(fit_errors_raw);
  return fit_params;
}

/*
  Loads Ay Moss columns for p-4He@200MeV.
*/
std::unique_ptr<ROOT::Math::Interpolator>
s13::ana::LoadHe4MossAy(std::string filename) {
  std::ifstream ifile(filename);
  std::vector<Double_t> angles, ays;

  for (s13::ana::CSVIterator iter(ifile, 2);
       iter != s13::ana::CSVIterator(); iter++) {
    auto csv_row = *iter;
    assert(csv_row.size() == 2 && "Invalid input file format");

    Double_t angle = std::stod(csv_row[0]) * D2R;
    Double_t ay = std::stod(csv_row[1]);

    angles.push_back(angle);
    ays.push_back(ay);
  }

  auto interp = new ROOT::Math::Interpolator{
    static_cast<unsigned int>(angles.size()),
    ROOT::Math::Interpolation::kLINEAR};
  interp->SetData(angles, ays);

  using interp_ptr_t = std::unique_ptr<ROOT::Math::Interpolator>;
  return interp_ptr_t(interp);
}

/*
  Loads Ay Moss columns for p-4He@200MeV.
*/
std::unique_ptr<ROOT::Math::Interpolator>
s13::ana::LoadHe4MossLabCs(std::string filename) {
  std::ifstream ifile(filename);
  std::vector<Double_t> angles, cses;

  for (s13::ana::CSVIterator iter(ifile, 2);
       iter != s13::ana::CSVIterator(); iter++) {
    auto csv_row = *iter;
    assert(csv_row.size() == 5 && "Invalid input file format");

    Double_t angle = std::stod(csv_row[0]) * D2R;
    Double_t cs = std::stod(csv_row[4]);

    angles.push_back(angle);
    cses.push_back(cs);
  }

  auto interp = new ROOT::Math::Interpolator{
    static_cast<unsigned int>(angles.size()),
    ROOT::Math::Interpolation::kLINEAR};
  interp->SetData(angles, cses);

  using interp_ptr_t = std::unique_ptr<ROOT::Math::Interpolator>;
  return interp_ptr_t(interp);
}

/*
  Loads Ay Moss columns for p-4He@200MeV.
*/
std::unique_ptr<ROOT::Math::Interpolator>
s13::ana::LoadHe4MossCmsCs(std::string filename) {
  std::ifstream ifile(filename);
  std::vector<Double_t> angles, cses;

  for (s13::ana::CSVIterator iter(ifile, 2);
       iter != s13::ana::CSVIterator(); iter++) {
    auto csv_row = *iter;
    assert(csv_row.size() == 2 && "Invalid input file format");

    Double_t angle = std::stod(csv_row[0]) * D2R;
    Double_t cs = std::stod(csv_row[1]);

    angles.push_back(angle);
    cses.push_back(cs);
  }

  auto interp = new ROOT::Math::Interpolator{
    static_cast<unsigned int>(angles.size()),
    ROOT::Math::Interpolation::kLINEAR};
  interp->SetData(angles, cses);

  using interp_ptr_t = std::unique_ptr<ROOT::Math::Interpolator>;
  return interp_ptr_t(interp);
}

std::unique_ptr<ROOT::Math::Interpolator>
s13::ana::LoadRdcEtoThetaConv(std::string filename) {
  std::ifstream ifile(filename);
  std::vector<Double_t> energies, thetas;

  for (s13::ana::CSVIterator iter(ifile, 1);
       iter != s13::ana::CSVIterator(); iter++) {
    auto csv_row = *iter;
    assert(csv_row.size() == 2 && "Invalid input file format");

    Double_t e = std::stod(csv_row[0]);
    Double_t a = std::stod(csv_row[1]);

    energies.push_back(e);
    thetas.push_back(a);
  }

  using interp_t = ROOT::Math::Interpolator;
  using interp_ptr_t = std::unique_ptr<interp_t>;
  interp_ptr_t interp(new interp_t(energies, thetas));
  return interp;
}

/*
  Function builds a TGraph with Moss CS columns for p-4He@200Mev.

  Params:
  filename: Name of the CSV file with Moss columns.
*/
TGraph* s13::ana::BuildMossLabCsGraph(const std::string &filename) {
  const int k_num_points = 30;
  const Double_t k_start_angle = 55.1;
  const Double_t k_finish_angle = 69.9;
  const Double_t k_step = (k_finish_angle - k_start_angle) / k_num_points;

  auto moss_cs = s13::ana::LoadHe4MossLabCs(filename);
  Double_t* angles = new Double_t[k_num_points];
  Double_t* cses = new Double_t[k_num_points];
  int arr_idx = 0;
  for (Double_t angle = k_start_angle; angle < k_finish_angle; angle += k_step) {
    Double_t angle_rad = angle * D2R;
    angles[arr_idx] = angle;
    cses[arr_idx] = moss_cs->Eval(angle_rad);

    std::cout << "angle: " << angle << "; cs: " << cses[arr_idx] <<
      std::endl;
    arr_idx +=1;
  }
  std::cout << "num points: " << arr_idx << std::endl;

  TGraph* graph = new TGraph(k_num_points, angles, cses);
  graph->SetLineColor(kBlue);
  graph->SetLineWidth(2);
  graph->SetMarkerColor(kBlack);
  graph->SetMarkerStyle(19);
  graph->SetTitle("CS p-4He@200MeV (Moss)");
  graph->GetXaxis()->SetTitle("Lab angle [deg]");
  graph->GetYaxis()->SetTitle("CS");

  return graph;
};

/*
  Function builds a TGraph with Moss CS columns for p-4He@200Mev.

  Params:
  filename: Name of the CSV file with Moss columns.
*/
TGraph* s13::ana::BuildMossCmsCsGraph(const std::string &filename) {
  misc::MessageLogger log("BuildMossCmsCsGraph()");
  log.debug("Building Moss CM CS graph...");
  auto moss_cs = s13::io::load_csv(filename, 2, 2);
  assert(moss_cs.size() == 2);

  TGraph* graph = new TGraph(moss_cs[0].size(), &moss_cs[0][0], &moss_cs[1][0]);
  graph->SetLineColor(kBlue);
  graph->SetLineWidth(2);
  graph->SetMarkerColor(47);
  graph->SetMarkerStyle(28);
  graph->SetTitle("CS p-4He@200MeV (Moss)");
  graph->GetXaxis()->SetTitle("CM angle [deg]");
  graph->GetYaxis()->SetTitle("CS [mb/sr]");

  log.debug("returning Moss CS.");
  return graph;
};

/*
  Function builds a TGraph with Moss CS columns for p-4He@200Mev.

  Params:
  filename: Name of the CSV file with Moss columns.
*/
TGraph* s13::ana::BuildMossAyGraph(const std::string& filename) {
  auto moss_ay = s13::io::load_csv(filename, 2, 2);
  assert(moss_ay.size() == 2);

  TGraph* graph = new TGraph(moss_ay[0].size(), &moss_ay[0][0], &moss_ay[1][0]);
  graph->SetLineColor(2);
  graph->SetLineWidth(4);
  graph->SetMarkerColor(4);
  graph->SetMarkerStyle(21);
  graph->SetTitle("Moss p-4He@200MeV Ay");
  graph->GetXaxis()->SetTitle("Lab angle [deg]");
  graph->GetYaxis()->SetTitle("Ay");
  // graph->SetMinimum(-1.2);
  // graph->SetMaximum(1.2);

  return graph;
};

TGraph* s13::ana::BuildHe6TheorCs(const std::string& filename) {
  s13::misc::MessageLogger log("BuildHe6TheorCs()");
  s13::io::CsvColumnsVector columns = io::load_csv(filename, 2);
  double xes[columns.at(0).size()], ys[columns.at(0).size()];
  log.debug("Loading CSV columns...");
  for (size_t i = 0; i < columns.at(0).size(); ++i) {
    xes[i] = columns.at(0).at(i);
    ys[i] = columns.at(1).at(i);
    log.debug("\tangle: %.1f; CS: %.1f mbarns", xes[i], ys[i]);
  }

  TGraph* graph = new TGraph(columns.at(0).size(), xes, ys);
  graph->SetLineColor(2);
  graph->SetLineWidth(4);
  graph->SetMarkerColor(4);
  graph->SetMarkerStyle(21);
  graph->SetTitle("S13 proposal p-6He@200MeV");
  graph->GetXaxis()->SetTitle("Lab angle [deg]");
  graph->GetYaxis()->SetTitle("Ay");

  return graph;
}


/*
  Function returns a graph of 6He-p kinematics
*/
TGraph* s13::ana::BuildHe6ProtonKinematicsGraph() {
  auto he6_kin_interp = s13::ana::LoadKinematicsData("data/he6_kin.csv");
  std::vector<double> p_thetas, he_thetas;
  for (int th = 50; th < 75; th += 1) {
    p_thetas.push_back(th);
    he_thetas.push_back(he6_kin_interp->Eval(th));
  }

  TGraph* graph = new TGraph(p_thetas.size(), &he_thetas[0], &p_thetas[0]);
  graph->SetLineColor(kRed);
  graph->SetLineWidth(3);
  graph->SetMarkerColor(5);
  graph->SetMarkerStyle(22);
  graph->SetTitle("6he-p kinematics");
  graph->GetXaxis()->SetTitle("He #theta [deg]");
  graph->GetYaxis()->SetTitle("Proton #theta [deg]");
  graph->SetMinimum(50);
  graph->SetMaximum(75);

  return graph;
};

TGraph* s13::ana::BuildHe4ProtonKinematicsGraph() {
  auto he4_kin_interp = s13::ana::LoadKinematicsData("data/he4_kin.csv");
  std::vector<double> p_thetas, he_thetas;
  for (int th = 50; th < 75; th += 1) {
    p_thetas.push_back(th);
    he_thetas.push_back(he4_kin_interp->Eval(th));
  }

  TGraph* graph = new TGraph(p_thetas.size(), &he_thetas[0], &p_thetas[0]);
  graph->SetLineColor(kRed);
  graph->SetLineWidth(3);
  graph->SetMarkerColor(5);
  graph->SetMarkerStyle(22);
  graph->SetTitle("4he-p kinematics");
  graph->GetXaxis()->SetTitle("He #theta [deg]");
  graph->GetYaxis()->SetTitle("Proton #theta [deg]");
  graph->SetMinimum(50);
  graph->SetMaximum(75);

  return graph;
};

void
s13::ana::SetRootStyle(s13::ana::RootStyle style) {
  switch (style) {
  case s13::ana::RootStyle::publication:
    gROOT->SetStyle("Pub");
    gStyle->SetLabelFont(132, "xyz");
    gStyle->SetTitleFont(132, "xyz");
    gStyle->SetLegendFont(132);
    gStyle->SetLabelSize(0.055, "xyz");
    gStyle->SetTitleSize(0.055, "xyz");
    gStyle->SetPadLeftMargin(0.27);
    gStyle->SetPadBottomMargin(0.20);
    gStyle->SetTitleOffset(1.9, "y");
    // gStyle->SetLegendTextSize(0.04);
    break;
  default:
    throw std::invalid_argument("unknown style passed!");
  }
}

s13::ana::TInterpPtr
s13::ana::LoadKinematicsData(std::string filename) {
  auto csv_data = io::load_csv(filename, 2, 2);
  auto interp =
    new TInterp(csv_data[0].size(), ROOT::Math::Interpolation::kLINEAR);
  interp->SetData(csv_data[0], csv_data[1]);
  return TInterpPtr{interp};
}
