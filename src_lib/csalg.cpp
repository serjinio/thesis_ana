

#include "csalg.hpp"


// definition of static members
int s13::ana::SignalBgMonitorAlg::hist_counter = 0;

/*
  Factory for SignalBgMonitorAlg
*/
std::shared_ptr<s13::ana::SignalBgMonitorAlg>
s13::ana::make_signal_bg_monitor_alg(int bin_num,
                                     double theta_range_start,
                                     double theta_range_end,
                                     s13::ana::ElasticScatteringCuts& cuts) {
  std::shared_ptr<s13::ana::SignalBgMonitorAlg>
    p_mon(new SignalBgMonitorAlg(bin_num, theta_range_start, theta_range_end, cuts));
  return p_mon;
}

std::shared_ptr<s13::ana::SignalBgMonitorAlg>
s13::ana::make_signal_bg_monitor_alg_from_file(std::istream& is) {
  std::shared_ptr<s13::ana::SignalBgMonitorAlg>
    p_mon(new SignalBgMonitorAlg());
  p_mon->deserialize(is);
  return p_mon;
}

std::shared_ptr<s13::ana::ElasticScatteringYieldsAlg>
s13::ana::make_elastic_scattering_yield_alg_from_file(std::istream& is) {
  std::shared_ptr<s13::ana::ElasticScatteringYieldsAlg>
    p_alg(new ElasticScatteringYieldsAlg());
  p_alg->deserialize(is);
  return p_alg;
}

std::shared_ptr<s13::ana::ElasticScatteringYieldsAlg>
s13::ana::make_elastic_scattering_yield_alg(int bin_num,
                                            double theta_range_start,
                                            double theta_range_end,
                                            s13::ana::ElasticScatteringCuts& cuts) {
  std::shared_ptr<s13::ana::ElasticScatteringYieldsAlg>
    alg(new s13::ana::ElasticScatteringYieldsAlg(bin_num, theta_range_start,
                                                 theta_range_end, cuts));
  return alg;
}
