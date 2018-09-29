

#include <memory>

#include "common.hpp"


namespace s13 {
  namespace ana{

    std::shared_ptr<ROOT::Math::Interpolator>
    make_theta_sigma_interp(s13::ana::DatasetType ds_type,
                            Dsdcs dsdc_selector, Espri espri_selector);

    std::shared_ptr<ROOT::Math::Interpolator>
    make_theta_mean_interp(s13::ana::DatasetType ds_type,
                           Dsdcs dsdc_selector, Espri espri_selector);

  }
}
