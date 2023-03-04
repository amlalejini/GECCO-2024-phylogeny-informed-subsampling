#pragma once

#include <string>
#include <functional>
#include <map>

#include "emp/base/Ptr.hpp"
#include "emp/datastructs/map_utils.hpp"

#include "psb/TestCaseSet.hpp"
#include "psb/readers/BaseProblemReader.hpp"

#include "problems/problems.hpp"

namespace psynth {

class ProblemManager {
public:

  template<typename PROBLEM_T>
  static std::function<void(ProblemManager& manager)> BuildProblemSetupFunction() {
    return [](ProblemManager& manager) {
      manager.ConfigureProblem<PROBLEM_T>();
    };
  }

protected:

  emp::Ptr<psb::BaseTestCaseSet> training_set;
  emp::Ptr<psb::BaseTestCaseSet> testing_set;
  // emp::Ptr<psb::BaseProblemUtility> problem_util;

  std::function<void()> cleanup = [](){ ; };
  std::function<void(const std::string&)> load_training_set;
  std::function<void(const std::string&)> load_testing_set;

  bool configured = false;

  /// Directory of available problems.
  std::map<
    std::string,
    std::function<void(ProblemManager& manager)>
  > problem_dir = {
    {"small-or-large", BuildProblemSetupFunction<problems::SmallOrLarge>()},
    {"median", BuildProblemSetupFunction<problems::Median>()}
  };

  template<typename PROBLEM_T>
  void ConfigureProblem() {
    using READER_T = typename PROBLEM_T::reader_t;
    // Cleanup any previous version of the problem.
    cleanup();

    // Configure cleanup for this version.
    cleanup = [this]() {
      if (training_set != nullptr) {
        auto training_set_ptr = training_set.Cast<psb::TestCaseSet<READER_T>>();
        training_set_ptr.Delete();
        training_set = nullptr;
      }
      if (testing_set != nullptr) {
        auto testing_set_ptr = testing_set.Cast<psb::TestCaseSet<READER_T>>();
        testing_set_ptr.Delete();
        testing_set = nullptr;
      }
    };

    // Configure training set load function
    load_training_set = [this](const std::string& filename) {
      auto training_set_ptr = training_set.Cast<psb::TestCaseSet<READER_T>>();
      training_set_ptr->LoadTests(filename);
    };

    // Configure testing set load function
    load_testing_set = [this](const std::string& filename) {
      auto testing_set_ptr = testing_set.Cast<psb::TestCaseSet<READER_T>>();
      testing_set_ptr->LoadTests(filename);
    };

    // Create new training/testing sets
    training_set = emp::NewPtr<psb::TestCaseSet<READER_T>>();
    testing_set = emp::NewPtr<psb::TestCaseSet<READER_T>>();

    configured = true;
  }

public:

  ~ProblemManager() {
    cleanup();
  }

  void ConfigureProblem(const std::string& problem_name) {
    emp_assert(emp::Has(problem_dir, problem_name), "Unknown problem", problem_name);
    problem_dir[problem_name](*this);
  }

  void LoadTestingSet(const std::string& filename) {
    emp_assert(configured);
    load_training_set(filename);
  }

  void LoadTrainingSet(const std::string& filename) {
    emp_assert(configured);
    load_testing_set(filename);
  }

  size_t GetTestingSetSize() const {
    return testing_set->GetSize();
  }

  size_t GetTrainingSetSize() const {
    return training_set->GetSize();
  }

  // template<typename INST_LIB_T>
  // void AddProblemInstructions(INST_LIB_T& inst_lib) {

  // }
  // TODO - load input into SGP hardware
  // TODO - check output
  // TODO - configure instruction set

  // TODO - name to functionality map

  bool IsValidProblem(const std::string& problem_name) {
    return emp::Has(problem_dir, problem_name);
  }

};




}