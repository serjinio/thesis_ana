
#pragma once

#include <assert.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>

#include <TTree.h>
#include <TChain.h>
#include <TGraph.h>
#include <TCut.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TPaveText.h>
#include "Math/Interpolator.h"

#include "msgoutput.hpp"


namespace s13 {

  namespace ana {

    using TInterp = ROOT::Math::Interpolator;
    using TInterpPtr = std::shared_ptr<TInterp>;

    /*
      Defines possible selections of ESPRI detectors
    */
    enum class Espri {left, right, both};

    std::string espri_selector_to_str(Espri espri_sel);

    /*
      Defines DSDCs which can be used as a source of angular data for He
    */
    enum class Dsdcs {s1dc, fdc0};

    /*
      Defines possible dataset types
    */
    enum class DatasetType {he4, he6, carbon};

    std::string ds_type_to_str(DatasetType ds_type);

    /*
      Defines possible polarization choces
    */
    enum class PolarizationDirection {up, down, both};

    std::string pol_dir_to_str(PolarizationDirection pol_dir);

    /*
      Convenience function to format a TString
    */
    template <typename... Ts>
    TString tstrfmt(const char* fmtstr, Ts&&... args) {
      return TString::Format(fmtstr, std::forward<Ts>(args)...);
    }

    TH1* make_th1(int bins_no, double min_x, double max_x,
                  TString title = "",
                  TString xlabel = "X", TString ylabel = "Counts");

    TH2* make_th2(int nbins_x, double min_x, double max_x,
                  int nbins_y, double min_y, double max_y,
                  TString title = "",
                  TString xlabel = "X", TString ylabel = "Y");

    THStack* make_thstack(std::vector<TH1*> hists, TString title = "Histograms stack");

    THStack* draw_thstack(std::vector<TH1*> hists, TString title = "Histograms stack",
                          TString draw_opt = "");

    TPaveText* draw_text(TString text, double x1 = 0.1, double y1 = 0.8,
                         double x2 = 0.4, double y2 = 0.9,
                         TString draw_opts = "");

    /*
      A line read from CSV file and separated into tokens.
    */
    class CSVRow {

    public:
      std::string const& operator[](std::size_t index) const {
        return m_data[index];
      }

      std::size_t size() const {
        return m_data.size();
      }

      void clear() {
        m_data.clear();
      }

      void readNextRow(std::istream& str) {
        std::string line;
        std::getline(str, line);

        std::stringstream lineStream(line);
        std::string cell;

        m_data.clear();
        while(std::getline(lineStream, cell, ',')) {
          m_data.push_back(cell);
        }
      }

    private:
      std::vector<std::string> m_data;
    };

    std::istream& operator>>(std::istream& str, CSVRow& columns);


    /*
      Iterator to read columns from CSV file
    */
    class CSVIterator {

    public:
      typedef std::input_iterator_tag     iterator_category;
      typedef CSVRow                      value_type;
      typedef std::size_t                 difference_type;
      typedef CSVRow*                     pointer;
      typedef CSVRow&                     reference;

      CSVIterator(std::istream& str, int skip_nrows = 0) :
        m_str(str.good() ? &str:NULL) {
        for (int i = 0; i < skip_nrows; i++) {
          ++(*this);
        }
      }

      CSVIterator() : m_str(NULL) {}

      // Pre Increment
      CSVIterator& operator++() {
        if (m_str) {
          if (!((*m_str) >> m_row)) {
            m_str = NULL;
          }
        }
        return *this;
      }

      // Post increment
      CSVIterator operator++(int) {
        CSVIterator tmp(*this);
        ++(*this);
        return tmp;
      }

      CSVRow const& operator*() const {
        return m_row;
      }

      CSVRow const* operator->() const {
        return &m_row;
      }

      bool operator==(CSVIterator const& rhs) {
        return ((this == &rhs) ||
                ((this->m_str == NULL) && (rhs.m_str == NULL)));
      }

      bool operator!=(CSVIterator const& rhs) {
        return !((*this) == rhs);
      }

    private:
      std::istream* m_str;
      CSVRow m_row;
    };

    /*
      Class to make phi cut dependent on proton theta as a linear function
    */
    class LinearInterp {

    public:

      /*
        a & b params - coefficients for linear fn: y = ax + b
      */
      LinearInterp(Double_t a, Double_t b) :
        a_{a}, b_{b} {}

      /*
        Return phi distribution sigma for a given proton angle
      */
      Double_t eval(Double_t x) const {
        return a_ * x + b_;
      }

    private:
      Double_t a_;
      Double_t b_;
    };

    template<typename T, typename TParam>
    class NamedType {
    private:
      T v_;

    public:
      explicit NamedType(const T& value) : v_{value} {};

      template<typename T_ = T>
      explicit NamedType(T&& value,
                         typename std::enable_if<
                         !std::is_reference<T_>{}>::type = 0) :
        v_{std::move(value)} {};
      T& get() { return v_; }
      const T& get() const { return v_; }
    };

  }

    namespace io {

      using CsvColumnsVector = std::vector< std::vector<double> >;

      /*
        Loads CSV data into CsvColumnsVector data structure.
      */
      CsvColumnsVector
      load_csv(std::string filename, size_t columns_number, size_t skip_rows = 1);

    }

    namespace ana {

      static constexpr Double_t D2R = 0.0174533;

      /*
        Class serves as a container for fit params.
        Knows how to read & write them from file.
      */
      class FitParams {

      public:

        FitParams() :
          logger_{"FitParams"}, m_params_num{0} {};
        FitParams(size_t params_num) :
          logger_{"FitParams"}, m_params_num{params_num} {};

        FitParams(const FitParams& fit_params, size_t from, size_t to) :
          logger_{"FitParams"}, m_params_num{from - to} {
            if (to <= from) {
              throw std::invalid_argument("'to' param must be larger than 'from'");
            }

            for (int i = from; i < 3; ++i) {
              params.push_back(fit_params.params[i]);
              errors.push_back(fit_params.errors[i]);
            }
        }

        void SetParamsNum(size_t params_num) {
          m_params_num = params_num;
        }

        void SetErrors(const Double_t *errors_array) {
          errors.clear();
          for (size_t i = 0; i < m_params_num; i++) {
            errors.push_back(errors_array[i]);
          }
        }

        void SetParams(const Double_t *params_array) {
          params.clear();
          for (size_t i = 0; i < m_params_num; i++) {
            params.push_back(params_array[i]);
          }
        }

      public:
        std::vector<Double_t> params;
        std::vector<Double_t> errors;

      private:
        misc::MessageLogger logger_;
        size_t m_params_num;
      };

      std::ostream& operator<<(std::ostream& os, const FitParams& obj);
      std::istream& operator>>(std::istream& is, FitParams& obj);

      /*
        Fits passed 1D hist with a gaussian.
        Returns: FitParams
      */
      FitParams fit_hist(TH1* hist);

      /*
        Loads Ay Moss columns for p-4He@200MeV.
      */
      std::unique_ptr<ROOT::Math::Interpolator>
      LoadHe4MossAy(std::string filename);

      /*
        Loads Ay Moss columns for p-4He@200MeV.
      */
      std::unique_ptr<ROOT::Math::Interpolator>
      LoadHe4MossLabCs(std::string filename);

      /*
      Loads Ay Moss columns for p-4He@200MeV.
      */
      std::unique_ptr<ROOT::Math::Interpolator>
      LoadHe4MossCmsCs(std::string filename);

      std::unique_ptr<ROOT::Math::Interpolator>
      LoadRdcEtoThetaConv(std::string filename);

      class HeAbangFdc2YConverter {

      public:
        /*
          Constructs converter. filename - is a file with
          coefficients to convert He a&b angles into Y position
          on FDC2
        */
        HeAbangFdc2YConverter(std::string filename) :
          logger_{"HeAbangFdc2YConverter"}, coeffs_file_{filename} {
            load_coeffs();
        }

        double fdc2_ypos(double he_aang, double he_bang) {
          assert(std::abs(he_aang) < 10);
          double ypos = interp_coeffs_->Eval(he_aang) * he_bang
            + interp_consts_->Eval(he_aang);
          return ypos;
        }

      private:

        void load_coeffs() {
          logger_.debug("Loading conversion coefficients from file: %s...",
                        coeffs_file_.c_str());
          std::ifstream csv_fs(coeffs_file_);
          std::vector<double> he_aangs;
          std::vector<double> fit_coeffs;
          std::vector<double> fit_constants;

          while(csv_fs) {
            std::string line;
            if (!getline(csv_fs, line)) {
              break;
            }

            if (line[0] == '#') {
              logger_.trace("skipping comment line");
              continue;
            }

            logger_.trace("got a line from input file: %s", line.c_str());

            std::stringstream line_ss(line);
            std::vector<double> record;
            while (line_ss) {
              std::string field = "";
              if (!getline(line_ss, field, ',' )) break;
              if (field == "") break;
              logger_.trace("got a non-empty field: %s", field.c_str());
              record.push_back(std::stod(field));
            }
            assert(record.size() == 5);

            he_aangs.push_back(record.at(0));
            fit_constants.push_back(record.at(1));
            fit_coeffs.push_back(record.at(2));
            logger_.trace("Loaded data row (he_aang, fit_const, fit_coeff): "
                          "(%.1f, %.1f, %.1f)",
                          record.at(0), record.at(1), record.at(2));
          }

          interp_coeffs_.reset(new ROOT::Math::Interpolator());
          interp_coeffs_->SetData(he_aangs, fit_coeffs);
          interp_consts_.reset(new ROOT::Math::Interpolator());
          interp_consts_->SetData(he_aangs, fit_constants);
        }

        misc::MessageLogger logger_;

        std::string coeffs_file_;
        std::unique_ptr<ROOT::Math::Interpolator> interp_coeffs_;
        std::unique_ptr<ROOT::Math::Interpolator> interp_consts_;

      };

      /*
        Function builds a TGraph with Moss CS columns for p-4He@200Mev.

        Params:
          filename: Name of the CSV file with Moss columns.
      */
      TGraph* BuildMossLabCsGraph(const std::string &filename);

      /*
      Function builds a TGraph with Moss CS columns for p-4He@200Mev.

      Params:
        filename: Name of the CSV file with Moss columns.
      */
      TGraph* BuildMossCmsCsGraph(const std::string &filename);

      /*
        Function builds a TGraph with Moss CS columns for p-4He@200Mev.

        Params:
          filename: Name of the CSV file with Moss columns.
      */
      TGraph* BuildMossAyGraph(const std::string& filename);

      /*
        Builds a TGraph with theoretical He6 CS from S13 proposal.
      */
      TGraph* BuildHe6TheorCs(const std::string& filename);

      /*
        Function returns a graph of 6He-p kinematics
      */
      TGraph* BuildHe6ProtonKinematicsGraph();

      /*
        Function returns a graph of 6He-p kinematics
      */
      TGraph* BuildHe4ProtonKinematicsGraph();

      /*
        Different style presets for ROOT hists.
      */
      enum class RootStyle { publication };

      /*
        Sets up style of hist elements
      */
      void SetRootStyle(s13::ana::RootStyle style);

      /*
        Loads kinematics data
      */
      TInterpPtr LoadKinematicsData(std::string filename);
    }
}
