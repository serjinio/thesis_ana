

#include "scattyield.hpp"


int s13::ana::ScatteringYield::m_hist_counter = 0;


s13::ana::ScatteringYield
s13::ana::operator+(s13::ana::ScatteringYield lhs,
                    const s13::ana::ScatteringYield& rhs) {
  return lhs += rhs;
}

s13::ana::ScatteringYield
s13::ana::operator-(s13::ana::ScatteringYield lhs,
                    const s13::ana::ScatteringYield& rhs) {
  return lhs -= rhs;
}

s13::ana::ScatteringYield
s13::ana::operator/(s13::ana::ScatteringYield lhs,
                    const s13::ana::ScatteringYield& rhs) {
  return lhs /= rhs;
}

s13::ana::ScatteringYield
s13::ana::operator*(s13::ana::ScatteringYield lhs,
                    const s13::ana::ScatteringYield& rhs) {
  return lhs *= rhs;
}

std::ostream&
s13::ana::operator<<(std::ostream& os,
                     const s13::ana::ScatteringYield& obj) {
  os << "Yields in range: ("
     << obj.RangeStart() << ", " << obj.RangeEnd()
     << "). Total events: " << obj.TotalEventCount()
     << ". Bin counts:" << std::endl;
  os << "Underflow events: " << obj.GetUnderflow()
     << "; overflow events: " << obj.GetOverflow() << "." << std::endl;

  for (int i = 0; i < obj.GetBinsNumber(); i++) {
    os << "\tAngle: " << obj.GetBinArg(i)
       << "; yield: " << obj.GetBinValue(i)
       << "; stat. error: " << obj.GetBinError(i)
       << std::endl;
  }

  return os;
}

s13::ana::ScatteringYield
s13::ana::weighted_average(const ScatteringYield& lhs,
                           const ScatteringYield& rhs) {
  assert(rhs.GetBinsNumber() == lhs.GetBinsNumber());
  ScatteringYield res{lhs};
  res.Reset();

  for (int i = 0; i < lhs.GetBinsNumber(); i++) {
    double val = 1/lhs.GetBinError(i) * lhs.GetBinValue(i) +
      1/rhs.GetBinError(i) * rhs.GetBinValue(i);
    double norfact = 1/lhs.GetBinError(i) + 1/rhs.GetBinError(i);
    val /= norfact;
    double err = std::sqrt(2.) / norfact;
    res.SetBinValue(res.GetBinArg(i), val);
    res.SetBinError(res.GetBinArg(i), err);
  }

  return res;
}
