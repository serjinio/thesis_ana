
#include <utility>

#include "cuts_conf.hpp"


using namespace s13::io;


s13::ana::ElasticScatteringCuts
s13::io::parse_cuts_config(std::istream&& is) {
  return parse_cuts_config(is);
}

s13::ana::ElasticScatteringCuts
s13::io::parse_cuts_config(std::istream& is) {
  s13::misc::MessageLogger logger("parse_cuts_config");
  if (!is) {
    logger.error("Specified input stream is not valid!");
    throw std::invalid_argument("Specified input stream is not valid!");
  }

  s13::ana::ElasticScatteringCuts cuts;
  try {
    CutsConfigLexer lexer(is);
    CutsConfigParser parser(lexer);
    CutsConfigAst ast = parser.parse();
    logger.debug("parsed AST of cuts config file:\n%s", ast.to_str().c_str());
    cuts = ast.eval();
  } catch (s13::io::CutsConfigError& err) {
    logger.error("Error during parsing of cuts configs: %s", err.what());
    throw;
  }
  return cuts;
}

std::vector<s13::ana::ElasticScatteringCuts>
s13::io::parse_cuts_range_config(std::istream&& is) {
  return parse_cuts_range_config(is);
}

std::vector<s13::ana::ElasticScatteringCuts>
s13::io::parse_cuts_range_config(std::istream& is) {
  s13::misc::MessageLogger logger("parse_cuts_range_config");
  if (!is) {
    logger.error("Specified input stream is not valid!");
    throw std::invalid_argument("Specified input stream is not valid!");
  }

  std::vector<s13::ana::ElasticScatteringCuts> cuts;
  try {
    CutsConfigLexer lexer(is);
    CutsConfigParser parser(lexer);
    CutsConfigAst ast = parser.parse();
    logger.debug("parsed AST of cuts config file:\n%s", ast.to_str().c_str());
    cuts = ast.eval_range();
  } catch (s13::io::CutsConfigError& err) {
    logger.error("Error during parsing of cuts configs: %s", err.what());
    throw;
  }
  return cuts;
}
