

#include "cuts_conf.hpp"


int main() {
  std::string sample_input = "#siski\n\n\nprop1\t : val1\nprop2: val2";
  std::string input_fname = "sample_cuts.conf";
  std::fstream fstr(input_fname);
  std::stringstream is(sample_input);
  s13::io::parse_cuts_config(fstr);
}
