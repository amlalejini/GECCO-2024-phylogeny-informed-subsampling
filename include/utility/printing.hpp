#pragma once

#include <iostream>
#include <sstream>
#include "emp/base/vector.hpp"
#include "emp/tools/string_utils.hpp"

namespace utils {

template<typename T>
void PrintVector(std::ostream& os, const emp::vector<T>& vec, bool quotes=false) {
  if (quotes) os << "\"";
  os << "[";
  for (size_t i = 0; i < vec.size(); ++i) {
    if (i) os << ",";
    os << emp::to_string(vec[i]);
  }
  os << "]";
  if (quotes) os << "\"";
}


template<typename T>
void PrintMapping(
  std::ostream& os,
  const std::unordered_map<std::string, T>& mapping,
  const std::string& pair_sep = ":",
  const std::string& mapping_sep = " ",
  bool quotes = false
) {
  if (quotes) os << "\"";
  os << "{";
  size_t cnt = 0;
  for (const auto& pair : mapping) {
    os << pair.first << pair_sep << pair.second;
    ++cnt;
    if (cnt != mapping.size()) os << mapping_sep;
  }
  os << "}";
  if (quotes) os << "\"";
}

} // End utils namespace