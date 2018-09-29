/**
 *   \file cli.hpp
 *   \brief Declarations for command line options handling
 *
 */


#pragma once

#include <thread>
#include <boost/program_options.hpp>

#include "common.hpp"
#include "msgoutput.hpp"


namespace po = boost::program_options;


namespace s13 {
  namespace misc {

    /*
      Simple structure to hold parsed command line options.
    */
    class CliOptions
    {
    public:
      CliOptions() :
        require_cuts_conf_{true}, log_{"CliOptions"} {
      }

      po::options_description& opts_description() {
        return opts_desc_;
      }

      po::positional_options_description pos_opts_description() {
        return pos_opts_desc_;
      }

      s13::ana::DatasetType dataset_type() const {
        return ds_type_;
      }

      s13::ana::PolarizationDirection pol_direction() const {
        return pol_dir_;
      }

      std::string pol_direction_str() const {
        if (pol_dir_ == s13::ana::PolarizationDirection::up) {
          return "pol_up";
        } else if (pol_dir_ == s13::ana::PolarizationDirection::down) {
          return "pol_down";
        } else {
          return "pol_both";
        }
      }

      int dataset_size() const {
        return ds_size_;
      }

      std::string dataset_type_str() const {
        if (vm_.count("ds")) {
            return vm_["ds"].as<std::string>();
        } else {
          return "";
        }
      }

      s13::ana::Espri espri_selector() const {
        return espri_selector_;
      }

      std::string espri_selector_str() const {
        if (vm_.count("espri")) {
          return vm_["espri"].as<std::string>();
        } else {
          return "";
        }
      }

      int threads_num() const {
        return threads_num_;
      }

      bool use_mt() const {
        return use_mt_;
      }

      std::string cuts_config() const {
        return cuts_config_;
      }

      std::string input_filename() {
        return input_filename_;
      }

      std::string output_file_basename(std::string prefix = "", std::string suffix = "") {
        if (output_filename_ != "") {
          return prefix + output_filename_ + suffix;
        } else {
          return prefix + dataset_type_str() + suffix;
        }
      }

      bool require_cuts_conf_opt() {
        return require_cuts_conf_;
      }

      void require_cuts_conf_opt(bool value) {
        require_cuts_conf_ = value;
      }

      void parse_options(int argc, char** argv) {
        make_common_options_description();
        po::store(po::command_line_parser(argc, argv).options(opts_desc_)
                  .positional(pos_opts_desc_).run(), vm_);
        po::notify(vm_);

        parse_common_options();
        parse_user_options();
      }

    private:

      void make_common_options_description() {
        auto cuts_conf_val = po::value<std::string>();
        if (require_cuts_conf_) {
          cuts_conf_val->required();
        }
        opts_desc_.add_options()
          ("help", "produce help message")
          ("ds", po::value<std::string>(), "type of dataset [he4|he6]")
          ("pol", po::value<std::string>(), "polarization direction [up|down|both]")
          ("ds-size", po::value<int>(), "number of runs to use from the dataset")
          ("espri", po::value<std::string>(), "which ESPRI to use [both|left|right]")
          ("mt", "indicates to use multi-threading, number of threads to use "
           "will be determined basing on hardware")
          ("threads,j", po::value<int>(), "number of threads to use, "
           "use if you want to force usage of a certain number of threads")
          ("cuts-conf", cuts_conf_val, "filename for cuts configuration")
          ("infile", po::value<std::string>(), "input filename")
          ("out,f", po::value<std::string>(), "output filename");

        pos_opts_desc_.add("ds", 1);
      }

      void parse_common_options() {
        if (vm_.count("help")) {
          std::cout << "List of allowed options:\n";
          std::cout << opts_desc_ << "\n";
          exit(0);
        }

        if (vm_.count("ds")) {
          std::string dataset = vm_["ds"].as<std::string>();
          parse_dataset_type_option(dataset);
        }

        if (vm_.count("pol")) {
          std::string pol_dir_str = vm_["pol"].as<std::string>();
          parse_pol_direction_option(pol_dir_str);
        }

        if (vm_.count("ds-size")) {
          ds_size_ = vm_["ds-size"].as<int>();
        }

        if (vm_.count("espri")) {
          std::string espri_sel = vm_["espri"].as<std::string>();
          parse_espri_selector_option(espri_sel);
        }

        if (vm_.count("hodf-he6")) {
          if (dataset_type() == s13::ana::DatasetType::he4) {
            throw std::invalid_argument("'hodf-he6' option cannot be used "
                                        "with he4 dataset type");
          }
          hodf_he6_cut_ = vm_["hodf-he6"].as<bool>();
        }

        if (vm_.count("mt")) {
          use_mt_ = true;
          threads_num_ = std::thread::hardware_concurrency();
        }

        if (vm_.count("threads")) {
          threads_num_ = vm_["threads"].as<int>();
          use_mt_ = true;
        }

        if (vm_.count("cuts-conf")) {
          cuts_config_ = vm_["cuts-conf"].as<std::string>();
        }

        if (vm_.count("out")) {
          output_filename_ = vm_["out"].as<std::string>();
        }

        if (vm_.count("infile")) {
          input_filename_ = vm_["infile"].as<std::string>();
        }
      }

      void parse_dataset_type_option(std::string ds_arg) {
        std::string ds_type_he4 = "he4";
        std::string ds_type_he6 = "he6";
        std::string ds_type_carbon = "carbon";
        if (ds_arg == ds_type_he4) {
          ds_type_= s13::ana::DatasetType::he4;
        } else if (ds_arg == ds_type_he6) {
          ds_type_ = s13::ana::DatasetType::he6;
        } else if (ds_arg == ds_type_carbon) {
          ds_type_ = s13::ana::DatasetType::carbon;
        } else {
          throw std::invalid_argument("Command line parameter only accepts "
                                      "'he4' and 'he6' as parameter values");
        }
      }

      void parse_pol_direction_option(std::string ds_arg) {
        std::string pol_dir_up = "up";
        std::string pol_dir_down = "down";
        std::string pol_dir_both = "both";
        if (ds_arg == pol_dir_up) {
          pol_dir_= s13::ana::PolarizationDirection::up;
        } else if (ds_arg == pol_dir_down) {
          pol_dir_ = s13::ana::PolarizationDirection::down;
        } else if (ds_arg == pol_dir_both) {
          pol_dir_ = s13::ana::PolarizationDirection::both;
        } else {
          throw std::invalid_argument("Command line parameter only accepts "
                                      "'up', 'down' and 'both' as parameter values");
        }
      }
      void parse_espri_selector_option(std::string espri_arg) {
        std::string espri_left = "left";
        std::string espri_right = "right";
        std::string espri_both = "both";
        if (espri_arg == espri_left) {
          espri_selector_ = s13::ana::Espri::left;
        } else if (espri_arg == espri_right) {
          espri_selector_ = s13::ana::Espri::right;
        } else if (espri_arg == espri_both) {
          espri_selector_ = s13::ana::Espri::both;
        } else {
          throw std::invalid_argument("Command line paramenter for ESPRI "
                                      "selection only accepts 'left', 'right' "
                                      "and 'both' as parameter values.");
        }
      }

      /*
        Should be customized in derived classes to add custom options
      */
      virtual void parse_user_options() {};

      s13::ana::DatasetType ds_type_ = s13::ana::DatasetType::he4;
      s13::ana::PolarizationDirection pol_dir_ = s13::ana::PolarizationDirection::both;
      int ds_size_ = -1;
      s13::ana::Espri espri_selector_ = s13::ana::Espri::both;
      bool hodf_he6_cut_ = true;
      int threads_num_ = 1;
      bool use_mt_ = false;
      std::string cuts_config_ = "";
      std::string output_filename_ = "";
      std::string input_filename_ = "";

      po::options_description opts_desc_;
      po::positional_options_description pos_opts_desc_;
      po::variables_map vm_;

      // specifies if --cuts-conf should be a required argument
      bool require_cuts_conf_;

      MessageLogger log_;
    };

  }
}
