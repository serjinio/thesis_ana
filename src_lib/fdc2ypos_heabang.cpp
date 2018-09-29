

#include "fdc2ypos_heabang.hpp"


using Fdc2YposAlgPtr = std::shared_ptr<s13::ana::Fdc2YposHeAbangAlg>;


Fdc2YposAlgPtr
s13::ana::make_fdc2_heabang_alg(s13::ana::ElasticScatteringCuts& cuts) {
  Fdc2YposAlgPtr alg_ptr(new s13::ana::Fdc2YposHeAbangAlg(cuts));
  return alg_ptr;
}
