/**
 *   \file drawing.hpp
 *   \brief Functions to draw histograms & graphs
 *
 */

#pragma once

#include <TChain.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TLine.h>
#include <TEllipse.h>
#include <TGraphErrors.h>
#include <TCut.h>
#include <TPad.h>

#include "common.hpp"
#include "scattyield.hpp"


namespace s13 {
  namespace ana {


    ////////////////////////////////////////////////////////////
    // histogramming common functions
    ////////////////////////////////////////////////////////////


    void draw_marker_line(Double_t x1, Double_t y1, Double_t x2, Double_t y2,
                          Double_t width = 1);

    void draw_marker_line_ysym(Double_t abs_x, Double_t y1, Double_t y2,
                               Double_t width = 1);

    void draw_marker_circle(Double_t center_x, Double_t center_y, Double_t radius,
                            Double_t width);

    TH1* hdraw(TTree& tree, TString name, TString draw_cmd,
               TString binning, TCut cuts = "", TString title = "",
               TString xaxis_title = "", TString yaxis_title = "",
               TString draw_opts = "colz");

    TH1* draw_theta_corr(TChain& chain, const TCut& cuts,
                         TString title = "", TString draw_opts = "col", bool draw2d = true,
                         int bin_num = 100, TString name = "");

    TH1* draw_theta_corr_2d(TChain& chain, const TCut& cuts,
                            TString title = "", TString draw_opts = "col",
                            int bin_num = 100, TString name = "");

    TH1* draw_theta_corr_1d(TChain& chain, const TCut& cuts,
                            TString title = "", TString draw_opts = "col",
                            int bin_num = 100, TString name = "");

    TGraph*
    BuildAyGraph(ScatteringYield yield, TString title);

    TGraph*
    BuildAsymmetryGraphErrors(ScatteringYield yield,
                              TString title);

    TGraph*
    BuildTGraphErrors(ScatteringYield yield,
                      TString title,
                      TString xaxis_title = "Lab angle [deg]",
                      TString yaxis_title = "Value");

    TGraph*
    BuildCsGraphFromData(const s13::io::CsvColumnsVector& csv_data,
                         std::string title);

    TGraph* Build4HeAyGraph();

    TMultiGraph* PrepareAyMultiGraph(TString title,
                                     std::vector<TGraph*> graphs = {});

    void AddHe6TheorCses(TMultiGraph* mg);

    TH1* init_ecut_bin_hist(double bin_min_theta, double bin_max_theta,
                            double xlim = 1);

    TH1* init_phi_bin_hist(double bin_min_theta, double bin_max_theta,
                           double xlim = 90);

    TH1* init_vertR_bin_hist(double bin_min_theta, double bin_max_theta,
                             double rlim = 20);

  }
}
