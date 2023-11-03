#pragma once

namespace utils {

bool IsClose(double value, double target, double eps) {
  return (value <= (target + eps)) && (value >= (target - eps));
}

}