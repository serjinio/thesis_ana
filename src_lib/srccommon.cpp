
#include <assert.h>

#include <memory>
#include <vector>

#include "srccommon.hpp"


using namespace s13::common;


std::string s13::common::GetUserHomeDir() {
  auto *pw = getpwuid(getuid());
  return std::string{pw->pw_dir};
}

std::string s13::common::GetAnaRootDir() {
  return GetUserHomeDir() + "/s13_ana";
}

void CSVRow::readNextRow(std::istream& str) {
   std::string line;
   std::getline(str, line);

   std::stringstream lineStream(line);
   std::string cell;

   m_data.clear();
   while(std::getline(lineStream, cell, ',')) {
     m_data.push_back(cell);
   }
 }

std::istream& s13::common::operator>>(std::istream& str, CSVRow& data) {
  data.readNextRow(str);
  return str;
}

std::unique_ptr<ROOT::Math::Interpolator>
s13::common::LoadProtonKinE(std::string filename) {
  std::ifstream ifile(filename);
  std::vector<Double_t> p_angles, p_energies;

  for (CSVIterator iter(ifile, 2); iter != CSVIterator(); iter++) {
    auto csv_row = *iter;
    assert(csv_row.size() == 4 && "Invalid input file format");

    Double_t p_a = std::stod(csv_row[2]) * D2R;
    Double_t p_e = std::stod(csv_row[3]);

    p_angles.push_back(p_a);
    p_energies.push_back(p_e);
  }

  auto interp = new ROOT::Math::Interpolator{
    static_cast<unsigned int>(p_angles.size()),
    ROOT::Math::Interpolation::kLINEAR};
  interp->SetData(p_angles, p_energies);

  using interp_ptr_t = std::unique_ptr<ROOT::Math::Interpolator>;
  return interp_ptr_t(interp);
}

std::unique_ptr<ROOT::Math::Interpolator>
s13::common::LoadStoppingPower(std::string filename) {
  std::cout << "Loading stopping power file: " <<
    filename << "..." << std::endl;
  std::ifstream ifile(filename);
  if (!ifile) {
    throw std::invalid_argument("Invalid input filename specified!");
  }
  std::vector<Double_t> p_energies, stopping_powers;

  for (CSVIterator iter(ifile, 2); iter != CSVIterator(); iter++) {
    auto csv_row = *iter;
    assert(csv_row.size() == 4 && "Invalid input file format");

    Double_t p_e = std::stod(csv_row[0]);
    Double_t stp = std::stod(csv_row[3]);

    p_energies.push_back(p_e);
    stopping_powers.push_back(stp);
  }

  auto interp = new ROOT::Math::Interpolator{
    static_cast<unsigned int>(p_energies.size()),
    ROOT::Math::Interpolation::kLINEAR};
  interp->SetData(p_energies, stopping_powers);

  using interp_ptr_t = std::unique_ptr<ROOT::Math::Interpolator>;
  return interp_ptr_t(interp);
}
