#pragma once

#include <functional>
#include <utility>

#include "emp/base/Ptr.hpp"

#include "BaseProblemHardware.hpp"

namespace psynth {

template<typename TAG_T>
class ProgSynthHardwareComponent {
public:
  using tag_t = TAG_T;

protected:
  // Have a pointer<base-problem-specific-component> that is cast as necessary by problem
  emp::Ptr<BaseProblemHardware> prob_hw = nullptr;
  bool prob_hw_init = false;

  /// @brief Function used to cleanup the problem hardware. Defaults to naive delete.
  std::function<void()> cleanup_prob_hw = [this]() {
    if (prob_hw != nullptr) {
      prob_hw.Delete();
    }
  };

  std::function<void()> reset_prob_hw = [this]() { ; };

  tag_t input_tag;

public:
  ~ProgSynthHardwareComponent() {
    cleanup_prob_hw();
  }

  template<typename PROB_HW_T>
  void CreateProblemHardware() {
    // Cleanup existing problem hardware (if any)
    cleanup_prob_hw();
    // (Re-)Configure cleanup
    cleanup_prob_hw = [this]() {
      if (prob_hw == nullptr) return;
      auto hw_ptr = prob_hw.Cast<PROB_HW_T>();
      hw_ptr.Delete();
      prob_hw_init = false;
    };
    // (Re-)Configure reset
    reset_prob_hw = [this]() {
      auto hw_ptr = prob_hw.Cast<PROB_HW_T>();
      hw_ptr->Reset();
    };

    prob_hw = emp::NewPtr<PROB_HW_T>();
    prob_hw_init = true;

  }

  void Reset() {
    emp_assert(prob_hw_init);
    // Reset problem hardware
    reset_prob_hw();
  }

  template<typename PROB_HW_T>
  PROB_HW_T& GetProbHW() {
    emp_assert(prob_hw_init);
    return *(prob_hw.Cast<PROB_HW_T>());
  }

  void SetInputTag(const tag_t& tag) {
    input_tag = tag;
  }

  const tag_t& GetInputTag() {
    return input_tag;
  }

};

}