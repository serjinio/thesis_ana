
#pragma once

#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

#include <cassert>
#include <string>
#include <array>
#include <vector>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>

#include "Math/Interpolator.h"
#include <TMath.h>
#include <TVector3.h>


namespace s13 {
  namespace common {

    constexpr Double_t D2R = 0.0174533;

    std::string GetUserHomeDir();
    std::string GetAnaRootDir();

    /*
      A line read from CSV file and separated into tokens.
    */
    class CSVRow {

    public:
      inline std::string const& operator[](std::size_t index) const {
        return m_data[index];
      }

      inline std::size_t size() const {
        return m_data.size();
      }

      inline void clear() {
        m_data.clear();
      }

      void readNextRow(std::istream& str);

    private:
      std::vector<std::string> m_data;
    };

    std::istream& operator>>(std::istream& str, CSVRow& data);

    /*
      Iterator to read data from CSV file
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
        ++(*this);

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
      Class presents scattering kinematics based on tabulated data.
    */
    class ScatteringKinematics {

    public:
      using KinSeries = std::vector<Double_t>;

      ScatteringKinematics() :
        m_is_data_dirty{true},
        m_interpolator{1000, ROOT::Math::Interpolation::kLINEAR} {}

      ScatteringKinematics(KinSeries fragment_angles, KinSeries recoil_angles) :
        m_fragment_angles{fragment_angles}, m_recoil_angles{recoil_angles},
        m_is_data_dirty{true},
        m_interpolator{1000, ROOT::Math::Interpolation::kLINEAR} {}

      ScatteringKinematics(const ScatteringKinematics& other) :
        m_is_data_dirty{true},
        m_interpolator{1000, ROOT::Math::Interpolation::kLINEAR} {
          SetData(other.m_fragment_angles, other.m_recoil_angles);
        }

      ScatteringKinematics& operator=(const ScatteringKinematics& rhs) {
        SetData(rhs.m_fragment_angles, rhs.m_recoil_angles);
        return *this;
      }

      void SetData(KinSeries fragment_angles, KinSeries recoil_angles) {
        m_fragment_angles = fragment_angles;
        m_recoil_angles = recoil_angles;
        m_is_data_dirty = true;
      }

      void AddPoint(Double_t fragment_angle, Double_t recoil_angle) {
        m_fragment_angles.push_back(fragment_angle);
        m_recoil_angles.push_back(recoil_angle);
        m_is_data_dirty = true;
      }

      // Double_t GetRecoilAngle(Double_t fragment_angle) {
      //   return 0;
      // }

      Double_t GetFragmentAngle(Double_t recoil_angle) {
        if (m_is_data_dirty) {
          // std::cout << "setting data for interpolator, dataset length: "
          //           << m_recoil_angles.size() << std::endl;
          assert(m_recoil_angles.size() == m_fragment_angles.size() &&
                 "Invalid input dataset for interpolator!");
          m_interpolator.SetData(m_recoil_angles, m_fragment_angles);
          m_is_data_dirty = false;
        }

        return m_interpolator.Eval(recoil_angle);
      }

    private:
      KinSeries m_fragment_angles;
      KinSeries m_recoil_angles;
      bool m_is_data_dirty;
      ROOT::Math::Interpolator m_interpolator;
    };


    std::unique_ptr<ROOT::Math::Interpolator>
    LoadProtonKinE(std::string filename);

    std::unique_ptr<ROOT::Math::Interpolator>
    LoadStoppingPower(std::string filename);

  }
}
