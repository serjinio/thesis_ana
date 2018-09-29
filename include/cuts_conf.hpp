/**
 *   \file cuts_conf.hpp
 *   \brief Code for parsing cuts configuration files.
 *
 */


#pragma once

#include <assert.h>

#include <cctype>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>

#include <msgoutput.hpp>
#include <elacuts.hpp>


namespace s13 {
  namespace io {

    enum class CutsConfigTokenType {eof, na, newline, word, colon, hashtag};
    enum class CutsConfigAstNodeType {root, name, value, bin_op_assign};

    class CutsConfigError : public std::logic_error {
    public:
      explicit CutsConfigError(const std::string& what) :
        std::logic_error(what) {}
    };

    class CutsConfigAstNode {

    public:
      CutsConfigAstNode(CutsConfigAstNodeType type, std::string token) :
        token_{token}, type_{type} {}

      CutsConfigAstNode& add_child(CutsConfigAstNode child) {
        children_.push_back(child);
        return children_.at(children_.size() - 1);
      }

      int num_children() const {
        return children_.size();
      }

      bool is_leaf() const {
        return children_.empty();
      }

      std::string value() const {
        return token_;
      }

      CutsConfigAstNodeType type() const {
        return type_;
      }

      CutsConfigAstNode& operator[] (int idx) {
        return children_.at(idx);
      }

    private:

      std::string token_;
      CutsConfigAstNodeType type_;
      std::vector<CutsConfigAstNode> children_;
    };

    class CutsConfigAst {

    public:
      CutsConfigAst(std::string tree_name) :
        root_{CutsConfigAstNodeType::root, tree_name},
        logger_("CutsConfigAst") {}

      s13::ana::ElasticScatteringCuts eval() {
        s13::ana::ElasticScatteringCuts cuts;
        for (int i = 0; i < root().num_children(); i++) {
          auto& node = root()[i];
          if (node.type() != CutsConfigAstNodeType::bin_op_assign) {
            throw CutsConfigError("Malformed AST!");
          }
          eval_assign_op(node, cuts);
        }
        return cuts;
      }

      size_t
      count_node_values(const CutsConfigAstNode& node) {
        std::istringstream iss(node.value());
        std::string item;
        size_t count = 0;
        while(std::getline(iss, item, ',')) {
          count++;
        }
        logger_.trace("Getting AST node values count. Node str value: %s. "
                      "Counted values number: %d",
                      node.value().data(), count);
        return count;
      }

      std::pair<std::string, size_t>
      find_range_params() {
        for (int i = 0; i < root().num_children(); i++) {
          auto& node = root()[i];
          if (node.type() != CutsConfigAstNodeType::bin_op_assign &&
              node.num_children() != 2) {
            throw CutsConfigError("Malformed AST - binary assign AST node"
                                  " must have exactly 2 child nodes!");
          }
          size_t count = count_node_values(node[1]);
          if (count > 1) {
            logger_.info("Found range variable '%s' with %d values.",
                         node.value().c_str(), count);
            return std::pair<std::string, size_t>(node.value(), count);
          }
        }
        return std::pair<std::string, size_t>("", 0);
      }

      std::vector<s13::ana::ElasticScatteringCuts>
      eval_range() {
        std::vector<s13::ana::ElasticScatteringCuts> cuts_range;
        std::pair<std::string, size_t> range_params = find_range_params();
        if (std::get<1>(range_params) == 0) {
          // in case we cannot identify range variable, then just return one cut
          return {eval()};

          // throw CutsConfigError("Invalid AST for range evaluation - no "
          //                       "range variable was found!");
        }

        for (size_t k = 0; k < std::get<1>(range_params); k++) {
          s13::ana::ElasticScatteringCuts cuts;
          for (int i = 0; i < root().num_children(); i++) {
            auto& node = root()[i];
            if (node.type() != CutsConfigAstNodeType::bin_op_assign) {
              throw CutsConfigError("Malformed AST!");
            }
            eval_assign_op_range(node, cuts,
                                 std::get<0>(range_params), k);
          }
          cuts_range.push_back(cuts);
        }

        logger_.info("Evaluated cuts config with range variable, "
                     "returning %d cut objects.", cuts_range.size());
        return cuts_range;
      }

      CutsConfigAstNode& root() {
        return root_;
      }

      std::string to_str() {
        std::string retval;
        for (int i = 0; i < root().num_children(); i++) {
          auto& node = root()[i];
          if (node.type() != CutsConfigAstNodeType::bin_op_assign) {
            throw CutsConfigError("Malformed AST!");
          }
          retval += assign_op_to_str(node);
        }
        return retval;
      }

    private:

      void eval_assign_op(CutsConfigAstNode& node,
                          s13::ana::ElasticScatteringCuts& cuts) {
        assert(node.num_children() == 2 && "Malformed assign node.");
        assign_cuts_prop(node[0], node[1], cuts);
      }

      void eval_assign_op_range(CutsConfigAstNode& node,
                                s13::ana::ElasticScatteringCuts& cuts,
                                std::string range_var_name, size_t range_idx) {
        assert(node.num_children() == 2 && "Malformed assign node.");
        if (node[0].value() == range_var_name) {
          assign_cuts_prop(node[0], node[1], cuts, range_idx);
        } else {
          assign_cuts_prop(node[0], node[1], cuts);
        }
      }

      void assign_cuts_prop(CutsConfigAstNode& prop_name_node,
                            CutsConfigAstNode& prop_value_node,
                            s13::ana::ElasticScatteringCuts& cuts,
                            size_t value_idx = 0) {
        std::string pname = prop_name_node.value();
        std::string val = prop_value_node.value();
        auto bval = [this, &val] () -> bool {
          return get_bool_prop_value(val); };
        auto fval = [this, &val, value_idx] () -> float {
          return get_float_prop_value(val).at(value_idx); };


        if (pname == "reaction_trigger_cut") {
          cuts.apply_reaction_trigger_cut = bval();
        } else if (pname == "beam_trigger_cut") {
          cuts.apply_beam_trigger_cut = bval();
        } else if (pname == "beam_abang_cut") {
          cuts.apply_beam_abang_cut = bval();
        } else if (pname == "beam_aang_width") {
          cuts.beam_aang_width = fval();
        } else if (pname == "beam_bang_width") {
          cuts.beam_bang_width = fval();
        } else if (pname == "vertexXY_cut") {
          cuts.apply_vertexXY_cut = bval();
        } else if (pname == "vertexXY_elliptical_cut") {
          cuts.apply_vertexXY_elliptical_cut = bval();
        } else if (pname == "vertexXY_Xdisp") {
          cuts.vertexXY_Xdisp = fval();
        } else if (pname == "vertexXY_Ydisp") {
          cuts.vertexXY_Ydisp = fval();
        } else if (pname == "vertexXY_ellipse_xR") {
          cuts.vertexXY_ellipse_xR = fval();
        } else if (pname == "vertexXY_ellipse_yR") {
          cuts.vertexXY_ellipse_yR = fval();
        } else if (pname == "phi_cut") {
          cuts.apply_phi_cut = bval();
        } else if (pname == "theta_cut") {
          cuts.apply_theta_cut = bval();
        } else if (pname == "hodf_he6_cut") {
          cuts.apply_hodf_he6_cut = bval();
        } else if (pname == "hodf_he4_cut") {
          cuts.apply_hodf_he4_cut = bval();
        } else if (pname == "subtract_he6_carbon_bg") {
          cuts.subtract_he6_carbon_bg = bval();
        } else if (pname == "espri_aang_cut") {
          cuts.apply_espri_aang_cut = bval();
        } else if (pname == "use_espri_aang_corr") {
          cuts.use_espri_aang_corr = bval();
        } else if (pname == "espri_vertz_cut") {
          cuts.apply_espri_vertz_cut = bval();
        } else if (pname == "espri_min_e_de_cut") {
          cuts.apply_espri_min_e_de_cut = bval();
        } else if (pname == "espri_e_cut") {
          cuts.apply_espri_e_cut = bval();
        } else if (pname == "espri_e_cut_use_he6_rel_e_corr") {
          cuts.espri_e_cut_use_he6_rel_e_corr = bval();
        } else if (pname == "espri_e_cut_use_he4_rel_e_corr") {
          cuts.espri_e_cut_use_he4_rel_e_corr = bval();
        } else if (pname == "espri_y_cut") {
          cuts.apply_espri_y_cut = bval();
        } else if (pname == "do_user_cut") {
          cuts.apply_user_cut = bval();
        } else if (pname == "vertexXY_radius") {
          cuts.vertexXY_radius = fval();
        } else if (pname == "phi_width") {
          if (val == "variable") {
            cuts.phi_width = 3;
            cuts.use_var_phi_cut = true;
          } else {
            cuts.phi_width = fval();
            cuts.use_var_phi_cut = true;
          }
        } else if (pname == "theta_width") {
          if (val == "variable") {
            cuts.theta_width = -1;
            cuts.use_var_theta_cut = true;
          } else {
            cuts.theta_width = fval();
            cuts.use_var_theta_cut = false;
          }
        } else if (pname == "espri_e_width_stdev") {
          cuts.espri_e_width_stdev = fval();
        } else if (pname == "espri_aang_width") {
          cuts.espri_aang_width = fval();
        } else if (pname == "espri_vertz_width") {
          cuts.espri_vertz_width = fval();
        } else if (pname == "espri_y_width") {
          cuts.espri_y_width = fval();
        } else if (pname == "espri_selector") {
          cuts.espri_selector = get_espri_selector_val(prop_value_node.value());
        } else if (pname == "dsdc_selector") {
          cuts.dsdc_selector = get_dsdc_selector_val(prop_value_node.value());
          if (cuts.dsdc_selector == s13::ana::Dsdcs::fdc0) {
            logger_.debug("setting DSDC selector to FDC0");
          } else if (cuts.dsdc_selector == s13::ana::Dsdcs::s1dc) {
            logger_.debug("setting DSDC selector to S1DC");
          }
        } else if (pname == "user_cut") {
          cuts.user_cut = prop_value_node.value().c_str();
        } else if (pname == "rdc_xnhitplanes_cut") {
          cuts.apply_rdc_xnhitplanes_cut = bval();
        } else if (pname == "rdc_xnhitplanes") {
          cuts.rdc_xnhitplanes = fval();
        } else if (pname == "rdc_ynhitplanes_cut") {
          cuts.apply_rdc_ynhitplanes_cut = bval();
        } else if (pname == "rdc_ynhitplanes") {
          cuts.rdc_ynhitplanes = fval();
        } else {
          std::string msg = "Unknown property name: " + pname;
          throw CutsConfigError(msg);
        }
      }

      bool get_bool_prop_value(std::string value) {
        if (value == "true") {
          return true;
        } else if (value == "false") {
          return false;
        } else {
          std::string msg = "Invalid boolean property value: " + value;
          throw CutsConfigError(msg);
        }
      }

      std::vector<double>
      get_float_prop_value(std::string values) {
        std::istringstream iss(values);
        std::string item;
        std::vector<std::string> items;
        while(std::getline(iss, item, ',')) {
          items.push_back(item);
        }
        if (items.size() == 0) {
          std::string msg = "Invalid float property value: " + values;
          throw CutsConfigError(msg);
        }

        logger_.trace("Getting float values from a node (str value: '%s'). "
                      "Items vector size: %d.", values.data(), items.size());
        std::vector<double> fvalues;
        for (size_t k = 0; k < items.size(); k++) {
          auto item_value = items[k];
          logger_.trace("\tGetting float value from string item: %s...",
                        item_value.data());
          for (size_t i = 0; i < item_value.size(); i++) {
            if (isalpha(item_value[i])) {
              std::string msg = "Invalid float property value: " + item_value;
              throw CutsConfigError(msg);
            }
          }
          double v = std::stod(item_value);
          logger_.trace("\tGot float value: %.2f.", v);
          fvalues.push_back(v);
        }
        return fvalues;
      }

      s13::ana::Espri get_espri_selector_val(std::string str_val) {
        if (str_val == "left") {
          return s13::ana::Espri::left;
        } else if (str_val == "right") {
          return s13::ana::Espri::right;
        } else if (str_val == "both") {
          return s13::ana::Espri::both;
        } else {
          std::string msg = "Invalid value for ESPRI selector property: "
            + str_val;
          throw CutsConfigError(msg);
        }
      }

      s13::ana::Dsdcs get_dsdc_selector_val(std::string str_val) {
        if (str_val == "s1dc") {
          return s13::ana::Dsdcs::s1dc;
        } else if (str_val == "fdc0") {
          return s13::ana::Dsdcs::fdc0;
        } else {
          std::string msg = "Invalid value for DSDC selector property: "
            + str_val;
          throw CutsConfigError(msg);
        }
      }

      std::string assign_op_to_str(CutsConfigAstNode& node) {
        if (node.num_children() != 2) {
          throw CutsConfigError("Malformed assgin node in AST!");
        }
        return node[0].value() + " : " + node[1].value() + "\n";
      }

      CutsConfigAstNode root_;
      s13::misc::MessageLogger logger_;
  };

    class CutsConfigLexer {

    public:
      CutsConfigLexer(std::istream& is) :
        is_{is}, line_is_empty_{true},
        token_type_{CutsConfigTokenType::na}, token_value_{""},
        is_eof_{false}, logger_{"CutsConfigLexer"} {
          if (!is.eof()) {
            consume();
          }
        }

      CutsConfigTokenType next_token() {
        if (is_eof_) {
          return token_type();
        }
        token_type_ = get_token_type();
        token_value_ = get_token_value();
        return token_type();
      }

      CutsConfigTokenType token_type() {
        return token_type_;
      }

      std::string token_name() {
        if (token_type_ == CutsConfigTokenType::eof) {
          return "eof";
        } else if (token_type_ == CutsConfigTokenType::colon) {
          return "colon";
        } else if (token_type_ == CutsConfigTokenType::hashtag) {
          return "hashtag";
        } else if (token_type_ == CutsConfigTokenType::newline) {
          return "newline";
        } else if (token_type_ == CutsConfigTokenType::word) {
          return "word";
        }
        return "n/a";
      }

      std::string token_value() {
        return token_value_;
      }

      bool is_eof() {
        return is_eof_;
      }

    private:

      void consume() {
        is_.get(next_char_);
      }

      void consume_line() {
        while (next_char_ != '\n') {
          consume();
        }
        consume(); // consume the newline char itself
      }

      CutsConfigTokenType get_token_type() {
        CutsConfigTokenType type;

        if (is_.eof()) {
          type = CutsConfigTokenType::eof;
          is_eof_ = true;
        } else if (next_char_ == ':') {
          type = CutsConfigTokenType::colon;
        } else if (next_char_ == '#') {
          consume_line();
          type = get_token_type();
        } else if (next_char_ == ' ' || next_char_ == '\t') {
          consume();
          type = get_token_type();
        } else if (next_char_ == '\n') {
          if (line_is_empty_) {
            consume();
            type = get_token_type();
          } else {
            type = CutsConfigTokenType::newline;
          }
        } else {
          type = CutsConfigTokenType::word;
        }

        if (type == CutsConfigTokenType::newline) {
          line_is_empty_ = true;
        } else {
          line_is_empty_ = false;
        }

        return type;
      }

      std::string get_token_value() {
        std::string retval = "";

        if (token_type() == CutsConfigTokenType::eof) {
          return "<eof>";
        } else if (token_type() != CutsConfigTokenType::word) {
          retval += next_char_;
          consume();
        } else {
          retval = get_word();
        }
        return retval;
      }

      std::string get_word() {
        std::string word = "";
        // logger_.debug("Getting word, next_char_ is: '%c'", next_char_);
        while ((isalnum(next_char_) ||
                next_char_ == '_' || next_char_ == '.' || next_char_ == ',' ||
                next_char_ == '-') && !is_.eof()) {
          word += next_char_;
          consume();
        }
        return word;
      }

      std::istream& is_;
      char next_char_;
      bool line_is_empty_;
      CutsConfigTokenType token_type_;
      std::string token_value_;
      bool is_eof_;

      s13::misc::MessageLogger logger_;
    };

    class CutsConfigParser {

    public:
      CutsConfigParser(CutsConfigLexer& lexer) :
        lexer_{lexer}, logger_{"CutsConfigParser"} {}

      CutsConfigAst parse() {
        CutsConfigAst ast("cuts");
        current_node_ = &ast.root();
        consume();
        while(!lexer_.is_eof()) {
          assign_r();
        }
        return ast;
      }

    private:

      void assign_r() {
        CutsConfigAstNode& rule_node =
          current_node_->add_child(CutsConfigAstNode(CutsConfigAstNodeType::bin_op_assign,
                                                    lexer_.token_value()));
        CutsConfigAstNode* prev_node = current_node_;
        current_node_ = &rule_node;
        prop_name_r(); colon_r(); value_r(); newline_r();
        current_node_ = prev_node;
      }

      void prop_name_r() {
        match(CutsConfigTokenType::word);
        current_node_->add_child(CutsConfigAstNode(CutsConfigAstNodeType::name,
                                                   lexer_.token_value()));
        consume();
      }

      void colon_r() {
        match(CutsConfigTokenType::colon);
        consume();
      }

      void value_r() {
        match(CutsConfigTokenType::word);
        current_node_->add_child(CutsConfigAstNode(CutsConfigAstNodeType::value,
                                                   lexer_.token_value()));
        consume();
      }

      void newline_r() {
        if (lexer_.token_type() == CutsConfigTokenType::newline
            || lexer_.token_type() == CutsConfigTokenType::eof) {
          consume();
        }
      }

      void consume() {
        lexer_.next_token();
        logger_.trace("Consuming token. Next token type: %s. Token value: '%s'",
                      lexer_.token_name().c_str(), lexer_.token_value().c_str());
      }

      void match(CutsConfigTokenType type) {
        if (lexer_.token_type() == type) {
          return;
        } else {
          throw CutsConfigError("Invalid syntax!");
        }
      }

      CutsConfigLexer& lexer_;

      CutsConfigAstNode* current_node_;

      s13::misc::MessageLogger logger_;
    };

    s13::ana::ElasticScatteringCuts
    parse_cuts_config(std::istream& is);

    s13::ana::ElasticScatteringCuts
    parse_cuts_config(std::istream&& is);

    std::vector<s13::ana::ElasticScatteringCuts>
    parse_cuts_range_config(std::istream& is);

    std::vector<s13::ana::ElasticScatteringCuts>
    parse_cuts_range_config(std::istream&& is);

  }
}
