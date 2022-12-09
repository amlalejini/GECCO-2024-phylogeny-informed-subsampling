#pragma once

#include <iostream>
#include <sstream>
#include "emp/base/vector.hpp"
#include "emp/tools/string_utils.hpp"

namespace utils {

template<typename T>
void PrintVector(std::ostream& os, const emp::vector<T>& vec) {
  os << "[";
  for (size_t i = 0; i < vec.size(); ++i) {
    if (i) os << ",";
    os << emp::to_string(vec[i]);
  }
  os << "]";
}

}