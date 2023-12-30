#pragma once

#include <string>
#include <functional>
#include <map>

#include "emp/base/Ptr.hpp"
#include "emp/datastructs/map_utils.hpp"

#include "psb/TestCaseSet.hpp"
#include "psb/readers/BaseProblemReader.hpp"

#include "ProgSynthOrg.hpp"
#include "TestResult.hpp"

#include "problems/problems.hpp"
#include "problems/BaseProblem.hpp"

namespace psynth {

template<typename HARDWARE_T>
class ProblemManager {
public:

  using sgp_hardware_t = HARDWARE_T;
  using inst_lib_t = typename sgp_hardware_t::inst_lib_t;
  using sgp_program_t = typename sgp_hardware_t::program_t;
  using org_t = ProgSynthOrg<sgp_program_t>;
  using event_lib_t = typename sgp_hardware_t::event_lib_t;

  template<typename PROBLEM_T>
  static std::function<void(ProblemManager& manager)> BuildProblemSetupFunction() {
    return [](ProblemManager& manager) {
      manager.ConfigureProblem<PROBLEM_T>();
    };
  }

protected:

  emp::Ptr<psb::BaseTestCaseSet> training_set = nullptr;
  emp::Ptr<psb::BaseTestCaseSet> testing_set = nullptr;
  emp::Ptr<problems::BaseProblem> problem_util = nullptr;

  std::function<void()> cleanup = [](){ ; };
  std::function<void(const std::string&)> load_training_set;
  std::function<void(const std::string&)> load_testing_set;
  std::function<void(inst_lib_t&)> add_problem_instructions;
  std::function<void(event_lib_t&)> add_problem_events;
  std::function<void(sgp_hardware_t&)> add_problem_hardware;
  std::function<void(sgp_hardware_t&, org_t&, size_t)> init_training_case;
  std::function<void(sgp_hardware_t&, org_t&, size_t)> init_testing_case;
  std::function<TestResult(sgp_hardware_t&, org_t&, size_t)> eval_output_training;
  std::function<TestResult(sgp_hardware_t&, org_t&, size_t)> eval_output_testing;
  std::function<double(void)> get_max_test_score;

  // std::string problem_name;
  bool configured = false;

  /// Directory of available problems.
  std::map<
    std::string,
    std::function<void(ProblemManager& manager)>
  > problem_dir = {
    {"small-or-large", BuildProblemSetupFunction<problems::SmallOrLarge>()},
    {"median", BuildProblemSetupFunction<problems::Median>()},
    {"grade", BuildProblemSetupFunction<problems::Grade>()},
    {"fizz-buzz", BuildProblemSetupFunction<problems::FizzBuzz>()},
    {"snow-day", BuildProblemSetupFunction<problems::SnowDay>()},
    {"smallest", BuildProblemSetupFunction<problems::Smallest>()},
    {"bouncing-balls", BuildProblemSetupFunction<problems::BouncingBalls>()},
    {"for-loop-index", BuildProblemSetupFunction<problems::ForLoopIndex>()},
    {"gcd", BuildProblemSetupFunction<problems::GCD>()},
    {"dice-game", BuildProblemSetupFunction<problems::DiceGame>()}
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
      if (problem_util != nullptr) {
        auto problem_util_ptr = problem_util.Cast<PROBLEM_T>();
        problem_util_ptr.Delete();
        problem_util_ptr = nullptr;
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

    // Configure add instructions function
    add_problem_instructions = [this](inst_lib_t& inst_lib) {
      auto problem_util_ptr = problem_util.Cast<PROBLEM_T>();
      problem_util_ptr->AddInstructions(inst_lib);
    };

    // Configure add events function
    add_problem_events = [this](event_lib_t& event_lib) {
      auto problem_util_ptr = problem_util.Cast<PROBLEM_T>();
      problem_util_ptr->AddEvents(event_lib);
    };

    problem_util = emp::NewPtr<PROBLEM_T>();

    // Configure function to add problem-specific hardware
    add_problem_hardware = [this](sgp_hardware_t& hw) {
      auto problem_util_ptr = problem_util.Cast<PROBLEM_T>();
      problem_util_ptr->ConfigureHardware(hw);
    };

    // Configure function to initialize a training example on a hardware unit
    init_training_case = [this](sgp_hardware_t& hw, org_t& org, size_t test_id) {
      auto training_set_ptr = training_set.Cast<psb::TestCaseSet<READER_T>>();
      auto problem_util_ptr = problem_util.Cast<PROBLEM_T>();
      problem_util_ptr->InitTest(hw, org, training_set_ptr->GetTest(test_id));
    };

    init_testing_case = [this](sgp_hardware_t& hw, org_t& org, size_t test_id) {
      auto testing_set_ptr = testing_set.Cast<psb::TestCaseSet<READER_T>>();
      auto problem_util_ptr = problem_util.Cast<PROBLEM_T>();
      problem_util_ptr->InitTest(hw, org, testing_set_ptr->GetTest(test_id));
    };

    eval_output_training = [this](
      sgp_hardware_t& hw,
      org_t& org,
      size_t test_id
    ) -> TestResult {
      auto training_set_ptr = training_set.Cast<psb::TestCaseSet<READER_T>>();
      auto problem_util_ptr = problem_util.Cast<PROBLEM_T>();
      return problem_util_ptr->EvaluateOutput(hw, org, training_set_ptr->GetTest(test_id));
    };

    eval_output_testing = [this](
      sgp_hardware_t& hw,
      org_t& org,
      size_t test_id
    ) -> TestResult {
      auto testing_set_ptr = testing_set.Cast<psb::TestCaseSet<READER_T>>();
      auto problem_util_ptr = problem_util.Cast<PROBLEM_T>();
      return problem_util_ptr->EvaluateOutput(hw, org, testing_set_ptr->GetTest(test_id));
    };

    get_max_test_score = [this]() -> double {
      auto problem_util_ptr = problem_util.Cast<PROBLEM_T>();
      return problem_util_ptr->max_test_score;
    };

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
    load_testing_set(filename);
  }

  void LoadTrainingSet(const std::string& filename) {
    emp_assert(configured);
    load_training_set(filename);
  }

  size_t GetTestingSetSize() const {
    emp_assert(configured);
    return testing_set->GetSize();
  }

  size_t GetTrainingSetSize() const {
    emp_assert(configured);
    return training_set->GetSize();
  }

  // template<typename INST_LIB_T>
  void AddProblemInstructions(inst_lib_t& inst_lib) {
    emp_assert(configured);
    add_problem_instructions(inst_lib);
  }

  void AddProblemEvents(event_lib_t& event_lib) {
    emp_assert(configured);
    add_problem_events(event_lib);
  }

  void AddProblemHardware(sgp_hardware_t& hw) {
    emp_assert(configured);
    add_problem_hardware(hw);
  }

  void InitTrainingCase(sgp_hardware_t& hw, org_t& org, size_t test_id) {
    emp_assert(configured);
    init_training_case(hw, org, test_id);
  }

  void InitTestingCase(sgp_hardware_t& hw, org_t& org, size_t test_id) {
    emp_assert(configured);
    init_testing_case(hw, org, test_id);
  }

  void InitCase(sgp_hardware_t& hw, org_t& org, size_t test_id, bool training) {
    emp_assert(configured);
    if (training) {
      init_training_case(hw, org, test_id);
    } else {
      init_testing_case(hw, org, test_id);
    }
  }

  TestResult EvaluateOutput(sgp_hardware_t& hw, org_t& org, size_t test_id, bool training) {
    emp_assert(configured);
    return (training) ? eval_output_training(hw, org, test_id) : eval_output_testing(hw, org, test_id);
  }

  double GetMaxTestScore() {
    emp_assert(configured);
    return get_max_test_score();
  }

  bool IsValidProblem(const std::string& problem_name) {
    return emp::Has(problem_dir, problem_name);
  }

};




}