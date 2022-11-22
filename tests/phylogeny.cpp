#define CATCH_CONFIG_MAIN

#include "Catch2/single_include/catch2/catch.hpp"
#include "Phylogeny.hpp"

TEST_CASE("Test phenotype_info", "[Phylogeny]") {

  std::cout << "Hello!" << std::endl;

  emp::vector<double> phenotype{0.0, 2.0, 0.0};
  emp::vector<bool> traits_evaluated{false, true, false};

  phenotype_info info;
  info.RecordPhenotype(
    phenotype,
    traits_evaluated
  );

  std::cout << "Phenotype: " << info.GetPhenotype() << std::endl;
  std::cout << "Evaluated: " << info.GetTraitsEvaluated() << std::endl;

  phenotype[0] = 1.0;
  traits_evaluated[0] = true;
  std::cout << "---" << std::endl;
  std::cout << "Phenotype: " << info.GetPhenotype() << std::endl;
  std::cout << "Evaluated: " << info.GetTraitsEvaluated() << std::endl;

};