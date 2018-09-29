/**
 *   \file scattevt.hpp
 *   \brief Definition of ScatteringEvent - used in TTree traversal
 *
 */


#pragma once

#include <assert.h>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>

#include <TFile.h>
#include <TTree.h>

#include <common.hpp>


namespace s13 {
  namespace ana {

    class ScatteringYield {

    public:

      ScatteringYield(int nbins, Double_t range_start_, Double_t range_end_,
                      Double_t initial_yield = 0) :
        bins_num_{nbins}, range_start_{range_start_}, range_end_{range_end_},
        overflow_bin_{0}, underflow_bin_{0},
        log_{"ScatteringYield"} {
          assert(range_start_ < range_end_);
          Reset(initial_yield);
        };

      ScatteringYield(TH1* hist) {
        assert(hist != nullptr);
        range_start_ = hist->GetXaxis()->GetXmin();
        range_end_ = hist->GetXaxis()->GetXmax();
        bins_num_ = hist->GetXaxis()->GetNbins();

        Reset(0);
        for (int i = 1; i <= bins_num_; i++) {
          SetBinValue(i, hist->GetBinContent(i));
        }
        overflow_bin_ = hist->GetBinContent(bins_num_ + 1);
        underflow_bin_ = hist->GetBinContent(0);
      }


      bool IsInRange(Double_t argument) {
        if (argument < range_start_ || argument > range_end_) {
          return false;
        } else {
          return true;
        }
      }

      int GetBinNo(Double_t argument) {
        double bin_no = (argument - range_start_) / bin_width_;
        if (bin_no < 0) {
          return -1;
        } else {
          return static_cast<int>(bin_no);
        }
      }

      int GetBinNoChecked(Double_t argument) {
        CheckArg(argument);
        return GetBinNo(argument);
      }

      void AddEvent(Double_t argument) {
        // log_.debug("Adding event with arg: %.2f", argument);
        int bin_no = GetBinNo(argument);
        //cerr << "bin no: " << bin_no << std::endl;
        if (bin_no < 0) {
          log_.trace("Got underflow event with argument: %.1f", argument);
          underflow_bin_ += 1;
          return;
        } else if (bin_no >= GetBinsNumber()){
          log_.trace("Got overflow event with argument: %.1f", argument);
          overflow_bin_ += 1;
          return;
        }

        if (ys_[bin_no] == 0) {
          ys_[bin_no] = 1;
          yserr_[bin_no] = 1;
        } else {
          ys_[bin_no] += 1;
          yserr_[bin_no] = std::sqrt(std::pow(yserr_[bin_no], 2) + 1);
        }
      }

      void CheckArg(Double_t argument) const {
        if (argument < range_start_ || argument > range_end_) {
          std::string msg("Event out of assigned range.");
          throw std::invalid_argument(msg);
        }
      }

      bool CheckArgNoThrow(Double_t argument) const {
        if (argument < range_start_ || argument > range_end_) {
          return false;
        } else {
          return true;
        }
      }

      void CheckBinIdx(int bin_idx) const {
        if (bin_idx >= bins_num_ || bin_idx < 0) {
          TString msg = TString::Format("Bin #%i is out of range", bin_idx);
          throw std::out_of_range(msg);
        }
      }

      Double_t GetBinWidth() const {
        return bin_width_;
      }

      int GetBinsNumber() const {
        return bins_num_;
      }

      double RangeStart() const {
        return range_start_;
      }

      double RangeEnd() const {
        return range_end_;
      }

      Double_t GetBinValue(Double_t argument) const {
        CheckArg(argument);
        int bin_no = static_cast<int>((argument - range_start_) / bin_width_);
        return ys_[bin_no];
      }

      Double_t GetBinValue(int bin_idx) const {
        CheckBinIdx(bin_idx);
        return ys_[bin_idx];
      }

      Double_t GetBinError(Double_t argument) const {
        CheckArg(argument);
        int bin_no = static_cast<int>((argument - range_start_) / bin_width_);
        return yserr_[bin_no];
      }

      Double_t GetBinError(int bin_idx) const {
        CheckBinIdx(bin_idx);
        return yserr_[bin_idx];
      }

      Double_t GetBinArg(int bin_idx) const {
        CheckBinIdx(bin_idx);
        return xs_[bin_idx];
      }

      Double_t GetOverflow() const {
        return overflow_bin_;
      }

      Double_t GetUnderflow() const {
        return underflow_bin_;
      }

      void SetBinValue(Double_t argument, Double_t value) {
        CheckArg(argument);
        int bin_no = static_cast<int>((argument - range_start_) / bin_width_);
        ys_[bin_no] = value;
        yserr_[bin_no] = std::sqrt(std::abs(value));
      }

      void SetBinError(Double_t argument, Double_t error) {
        CheckArg(argument);
        int bin_no = static_cast<int>((argument - range_start_) / bin_width_);
        yserr_[bin_no] = error;
      }

      void ScaleBin(int bin_idx, Double_t factor) {
        assert(bin_idx < bins_num_);
        ys_[bin_idx] *= factor;
        // TODO: verify, not sure if this is correct way to compute error
        yserr_[bin_idx] *= factor;
      }

      std::size_t Size() {
        assert(ys_.size() == xs_.size());
        return ys_.size();
      }

      Double_t TotalEventCount() const {
        Double_t total = 0;
        for_each(ys_.begin(), ys_.end(), [&total] (Double_t yield) {
          total += yield;
        });

        return total;
      }

      std::vector<double>::iterator
      ArgsRangeStartIter() {
        return xs_.begin();
      }

      std::vector<double>::iterator
      ArgsRangeEndIter() {
        return xs_.end();
      }

      std::vector<double>::iterator
      ValsRangeStartIter() {
        return ys_.begin();
      }

      std::vector<double>::iterator
      ValsRangeEndIter() {
        return ys_.end();
      }

      std::vector<double>::iterator
      ErrsRangeStartIter() {
        return yserr_.begin();
      }

      std::vector<double>::iterator
      ErrsRangeEndIter() {
        return yserr_.end();
      }

      ScatteringYield& operator=(const ScatteringYield& other) {
        bins_num_ = other.bins_num_;
        range_start_ = other.range_start_;
        range_end_ = other.range_end_;
        bin_width_ = other.bin_width_;
        xs_ = other.xs_;
        ys_ = other.ys_;
        yserr_ = other.yserr_;
        overflow_bin_ = other.overflow_bin_;
        underflow_bin_ = other.underflow_bin_;
        return *this;
      }

      ScatteringYield invert() const {
        ScatteringYield tmp = *this;
        for (int i = 0; i < bins_num_; i++) {
          tmp.yserr_[i] = tmp.yserr_[i] / std::pow(tmp.ys_[i], 2);
          tmp.ys_[i] = 1 / tmp.ys_[i];
        }

        tmp.overflow_bin_ = 1 / overflow_bin_;
        tmp.underflow_bin_ = 1 / underflow_bin_;

        return tmp;
      }

      ScatteringYield operator-() const {
        ScatteringYield tmp = *this;
        for (int i = 0; i < bins_num_; i++) {
          tmp.ys_[i] = -tmp.ys_[i];
        }

        tmp.overflow_bin_ = -overflow_bin_;
        tmp.underflow_bin_ = -underflow_bin_;

        return tmp;
      }

      ScatteringYield multiply(Double_t factor) const {
        ScatteringYield tmp = *this;
        for (int i = 0; i < bins_num_; i++) {
          tmp.ys_[i] = tmp.ys_[i] * factor;
          tmp.yserr_[i] = tmp.yserr_[i] * std::abs(factor);
        }

        tmp.overflow_bin_ *= factor;
        tmp.underflow_bin_ *= factor;

        return tmp;
      }

      ScatteringYield pow(Double_t power) const {
        ScatteringYield tmp = *this;
        for (int i = 0; i < bins_num_; i++) {
          tmp.ys_[i] = std::pow(tmp.ys_[i], power);
          tmp.yserr_[i] = std::abs(power * std::pow(tmp.ys_[i], power - 1)) *
            tmp.yserr_[i];
        }

        tmp.overflow_bin_ = std::pow(tmp.overflow_bin_, power);
        tmp.underflow_bin_ = std::pow(tmp.underflow_bin_, power);

        return tmp;
      }

      ScatteringYield abs() const {
        ScatteringYield tmp = *this;
        for (int i = 0; i < bins_num_; i++) {
          tmp.ys_[i] = std::abs(tmp.ys_[i]);
        }

        tmp.overflow_bin_ = std::abs(overflow_bin_);
        tmp.underflow_bin_ = std::abs(underflow_bin_);

        return tmp;
      }

      ScatteringYield sqrt() const {
        ScatteringYield tmp = *this;
        for (int i = 0; i < bins_num_; i++) {
          tmp.yserr_[i] = tmp.yserr_[i]
            / ( 2 * std::sqrt( std::abs(tmp.ys_[i]) ) );
          tmp.ys_[i] = std::sqrt(tmp.ys_[i]);
        }

        tmp.overflow_bin_ = std::sqrt(overflow_bin_);
        tmp.underflow_bin_ = std::sqrt(underflow_bin_);

        return tmp;
      }

      ScatteringYield& operator+=(const ScatteringYield& other) {
        assert(other.bins_num_ == bins_num_ &&
               other.range_start_ == range_start_ &&
               other.range_end_ == range_end_);
        for (int i = 0; i < bins_num_; i++) {
          ys_[i] += other.ys_[i];
          yserr_[i] = std::sqrt(std::pow(yserr_[i], 2) +
                               std::pow(other.yserr_[i], 2));
        }

        overflow_bin_ += other.overflow_bin_;
        underflow_bin_ += other.underflow_bin_;

        return *this;
      }

      ScatteringYield& operator-=(const ScatteringYield& other) {
        *this += -other;
        return *this;
      }

      ScatteringYield& operator*=(const ScatteringYield& other) {
        assert(other.bins_num_ == bins_num_ &&
               other.range_start_ == range_start_ &&
               other.range_end_ == range_end_);
        for (int i = 0; i < bins_num_; i++) {
          Double_t new_y = ys_[i] * other.ys_[i];
          // Double_t new_yerr =
          //   std::sqrt(std::pow(yserr_[i]/ys_[i], 2) +
          //             std::pow(other.yserr_[i]/other.ys_[i], 2))
          //   * std::abs(new_y);
          Double_t new_yerr = std::pow(yserr_[i] * other.ys_[i], 2) +
            std::pow(other.yserr_[i] * ys_[i], 2);
          new_yerr = std::sqrt(new_yerr);
          ys_[i] = new_y;
          yserr_[i] = new_yerr;
        }

        overflow_bin_ *= other.overflow_bin_;
        underflow_bin_ *= other.underflow_bin_;

        return *this;
      }

      ScatteringYield& operator/=(const ScatteringYield& other) {
        *this *= other.invert();
        return *this;
      }

      std::ostream& Serialize(std::ostream& os) {
        os << bins_num_ << "," << range_start_ << "," << range_end_ << std::endl;
        for (int i = 0; i < bins_num_; i++) {
          os << xs_[i] << "," << ys_[i] << "," << yserr_[i] << std::endl;
        }

        os << underflow_bin_ << "," << overflow_bin_ << std::endl;

        return os;
      }

      std::istream& Deserialize(std::istream& is) {
        xs_.clear();
        ys_.clear();
        yserr_.clear();
        overflow_bin_ = 0;
        underflow_bin_ = 0;

        int current_line = -1;
        while (is) {
          std::string line;
          std::getline(is, line);
          if (line == "") {
            continue;
          }

          std::stringstream lineStream(line);
          std::string cell;
          if (current_line == -1) {
            // first input line - range information
            std::getline(lineStream, cell, ',');
            bins_num_ = std::stod(cell);
            std::getline(lineStream, cell, ',');
            range_start_ = std::stod(cell);
            std::getline(lineStream, cell);
            range_end_ = std::stod(cell);
            bin_width_ = (range_end_ - range_start_) / bins_num_;
          } else if (current_line == bins_num_) {
            // last line - overflow & underflow bins information
            std::getline(lineStream, cell, ',');
            underflow_bin_ = std::stod(cell);
            std::getline(lineStream, cell, ',');
            overflow_bin_ = std::stod(cell);
          } else {
            // bins information
            std::getline(lineStream, cell, ',');
            double arg = std::stod(cell);
            assert(arg > range_start_ && arg < range_end_ &&
                   "Argument value out of range.");
            xs_.push_back(arg);
            std::getline(lineStream, cell, ',');
            ys_.push_back(std::stod(cell));
            std::getline(lineStream, cell);
            yserr_.push_back(std::stod(cell));
          }

          if (current_line == bins_num_) {
            break;
          }
          ++current_line;
        }


        return is;
      }

      void Reset(Double_t initial_yield = 0) {
        xs_.clear();
        ys_.clear();
        yserr_.clear();
        underflow_bin_ = 0;
        overflow_bin_ = 0;
        bin_width_ = (range_end_ - range_start_) / bins_num_;
        for (int i = 0; i < bins_num_; i++) {
          xs_.push_back(range_start_ + bin_width_ * i + bin_width_ / 2);
          ys_.push_back(initial_yield);
          yserr_.push_back(std::sqrt(initial_yield));
        }
      }

      ScatteringYield Empty() const {
        ScatteringYield tmp(*this);
        tmp.Reset();
        return tmp;
      }

      TH1* AsHist(TString title = "", TString xlabel = "X",
                  TString ylabel = "Counts") const {
        TH1* hist =
          s13::ana::make_th1(bins_num_, range_start_, range_end_,
                             title, xlabel, ylabel);
        for (int i = 1; i <= bins_num_; i++) {
          hist->SetBinContent(i, ys_[i - 1]);
          hist->SetBinError(i, yserr_[i - 1]);
        }

        hist->SetBinContent(0, underflow_bin_);
        hist->SetBinContent(bins_num_ + 1, overflow_bin_);

        return hist;
      }

      TH1* AsHist(double xmin, double xmax,
                  TString title = "",
                  TString xlabel = "X",
                  TString ylabel = "Counts") const {
        int hist_bins_num = 0;
        for (int i = 1; i <= bins_num_; i++) {
          double arg = GetBinArg(i - 1);
          if (arg < xmin || arg > xmax) continue;
          hist_bins_num += 1;
        }
        double hist_bin_width = (xmax - xmin) / hist_bins_num;

        TH1* hist =
          s13::ana::make_th1(hist_bins_num, xmin, xmax,
                             title, xlabel, ylabel);

        for (int i = 1; i <= bins_num_; i++) {
          double arg = GetBinArg(i - 1);
          int cur_bin = (arg - xmin) / hist_bin_width + 1;
          if (cur_bin < 1 || cur_bin > hist_bins_num) continue;

          hist->SetBinContent(cur_bin, ys_[i - 1]);
          hist->SetBinError(cur_bin, yserr_[i - 1]);
        }

        hist->SetBinContent(0, underflow_bin_);
        hist->SetBinContent(bins_num_ + 1, overflow_bin_);

        return hist;
      }

    private:

      int bins_num_;
      Double_t range_start_;
      Double_t range_end_;
      std::vector<Double_t> xs_;
      std::vector<Double_t> ys_;
      std::vector<Double_t> yserr_;
      double overflow_bin_;
      double underflow_bin_;

      static int m_hist_counter;

      Double_t bin_width_;

      s13::misc::MessageLogger log_;
    };

    ScatteringYield operator+(ScatteringYield lhs, const ScatteringYield& rhs);
    ScatteringYield operator-(ScatteringYield lhs, const ScatteringYield& rhs);
    ScatteringYield operator/(ScatteringYield lhs, const ScatteringYield& rhs);
    ScatteringYield operator*(s13::ana::ScatteringYield lhs,
                              const s13::ana::ScatteringYield& rhs);
    std::ostream& operator<<(std::ostream& os,
                             const s13::ana::ScatteringYield& obj);

    ScatteringYield
    weighted_average(const ScatteringYield& lhs, const ScatteringYield& rhs);

  }
}
