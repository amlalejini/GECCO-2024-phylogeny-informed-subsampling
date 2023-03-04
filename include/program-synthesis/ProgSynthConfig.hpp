#pragma once

#include "emp/config/config.hpp"

namespace psynth {

EMP_BUILD_CONFIG(ProgSynthConfig,
  GROUP(WORLD, "How should the world be setup?"),
  VALUE(SEED, int, 0, "Random number seed"),
  VALUE(POP_SIZE, size_t, 512, "Population size."),

  GROUP(PSB_PROBLEM, "Problem-related settings"),
  VALUE(PROBLEM, std::string, "small-or-large", "Problem to solve"),
  VALUE(TESTING_SET_PATH, std::string, "testing.json", "Path to testing set (json)"),
  VALUE(TRAINING_SET_PATH, std::string, "training.json", "Path to training set (json)"),

  GROUP(SGP_CPU, "SignalGP Virtual CPU"),
  VALUE(MAX_ACTIVE_THREAD_CNT, size_t, 8, "Maximum number of active threads that can run simultaneously on a SGP virtual CPU."),
  VALUE(MAX_THREAD_CAPACITY, size_t, 16, "Maximum thread capacity.")

)

}