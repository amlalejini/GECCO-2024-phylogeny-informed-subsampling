// This is the main function for the NATIVE version of this project.

#include <iostream>
#include <limits>

#include "emp/base/vector.hpp"
#include "emp/config/command_line.hpp"
#include "emp/config/ArgManager.hpp"

#include "program-synthesis/ProgSynthConfig.hpp"
#include "program-synthesis/ProgSynthWorld.hpp"

#include "psb/readers/TestCaseReaders.hpp"


int main(int argc, char* argv[])
{
  psynth::ProgSynthConfig config;
  config.Read("prog_synth.cfg", false);
  auto args = emp::cl::ArgManager(argc, argv);
  if (args.ProcessConfigOptions(config, std::cout, "prog_synth.cfg", "diagnostics-macros.h") == false) exit(0);
  if (args.TestUnknown() == false) exit(0);  // If there are leftover args, throw an error.

  std::cout << "==============================" << std::endl;
  std::cout << "|    How am I configured?    |" << std::endl;
  std::cout << "==============================" << std::endl;
  config.Write(std::cout);
  std::cout << "==============================\n"
            << std::endl;

  psynth::ProgSynthWorld world(config);
  // world.Run();

  // pm.ConfigureProblem("")
  // pm.ConfigureProblem<psb::readers::SmallOrLarge>();
  // pm.ConfigureProblem<psb::readers::Bowling>();
  // pm.ConfigureProblem<psb::readers::BouncingBalls>();

}