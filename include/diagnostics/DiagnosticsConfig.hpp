/**
 * Code adapted from https://github.com/jgh9094/ECJ-2022-suite-of-diagnostics-for-selection-schemes/ by Jose Hernandez
*/
#pragma once

#include "emp/config/config.hpp"

namespace diag {

EMP_BUILD_CONFIG(DiagnosticsConfig,
  GROUP(WORLD, "How should the world be setup?"),
  VALUE(POP_SIZE, size_t, 512, "Population size."),
  VALUE(MAX_GENS, size_t, 50000, "Maximum number of generations."),
  VALUE(SEED, int, 1, "Random number seed."),
  VALUE(INIT_POP_RAND, bool, true, "Do we start randomly (true) or from the lowest point (false)"),

  GROUP(DIAGNOSTICS, "How are the diagnostics setup?"),
  VALUE(DIAGNOSTIC, std::string, "exploitation", "Which diagnostic should be run? Options: \nexploitation \nstruct-exploitation \nstrong-ecology \nweak-ecology \nexploration"),
  VALUE(VALLEY_CROSSING, bool, false, "Do we add valley-crossing layer to diagnostics?"),
  VALUE(GENE_LOWER_BND, double, 0.0, "Lower bound for random starts."),
  VALUE(GENE_UPPER_BND, double, 1.0, "Upper bound for random starts."),
  VALUE(ACCURACY, double, 0.99, "Accuracy percentage needed to be considered an optimal trait"),
  VALUE(TARGET, double, 100.0, "Upper bound for random starts."),
  VALUE(CREDIT, double, 0.00, "Maximum credit a solution can get on an objective if applicable"),
  VALUE(DIAGNOSTIC_DIMENSIONALITY, size_t, 100, "Number of traits an organism has (i.e., genome/phenotype dimensionality"),

  GROUP(MUTATIONS, "Mutation rates for organisms."),
  VALUE(MUTATE_PER_SITE_RATE, double, 0.007, "Probability of instructions being mutated"),
  VALUE(MUTATE_MEAN, double, 0.0, "Mean of Gaussian Distribution for mutations"),
  VALUE(MUTATE_STD, double, 1.0, "Standard Deviation of Gaussian Distribution for mutations"),

  GROUP(EVALUTION, "How are organisms evaluated?"),
  VALUE(EVAL_MODE, std::string, "full", "Evaluation mode. Options:\nfull\ncohort\ndown-sample"),

  GROUP(EVALUATION_COHORT, "Cohort evaluation settings"),
  VALUE(NUM_COHORTS, size_t, 2, "How many cohorts should we divide the tests and organisms into?"),

  GROUP(EVALUATION_DOWNSAMPLE, "Down-sample evaluation settings"),
  VALUE(TEST_DOWNSAMPLE_RATE, double, 0.5, "Proportion of tests to use in test down-sample"),

  GROUP(SELECTION, "Selection scheme"),
  VALUE(SELECTION, std::string, "truncation", "Which selection are we doing? Options: \ntruncation \ntournament \nfitness-sharing \nlexicase \nlexicase-eps lexicase-even-lead \nnondominated-sorting \nnovelty"),

  GROUP(TRUNCATION, "Parameters for truncation."),
  VALUE(TRUNC_SIZE, size_t, 8, "Parameter estimate for truncation size t."),

  GROUP(TOURNAMENT, "Parameters for tournament."),
  VALUE(TOURNAMENT_SIZE, size_t, 8, "Parameter estimate for tournament size."),

  GROUP(FITSHARING, "Parameters for fitness sharing."),
  VALUE(FITNESS_SHARING_SIGMA, double, 1.0, "Parameter estimate for proportion of similarity threshold sigma (based on maximum distance between solutions)."),
  VALUE(FITNESS_SHARING_ALPHA, double, 1.0, "Parameter estimate for penalty function shape alpha."),
  VALUE(FITNESS_SHARING_APPLI, bool, false, "Fitness sharing applied: 0->Genome, 1->Phenotype"),

  GROUP(NOVELTY, "Parameters for novelty search."),
  VALUE(NOVELTY_K, size_t, 15, "Parameter estimate k-nearest neighbors."),
  VALUE(NOVELTY_PMIN, double, 10.0, "Minimum novelty score needed to enter archive."),
  VALUE(NOVELTY_UP, double, 0.25, "Increase pmin by this much."),
  VALUE(NOVELTY_DOWN, double, 0.05, "Decrease pmin by this much."),
  VALUE(NOVELTY_RI, double, 0.00001, "Probability of random solution inserted in archive."),
  VALUE(NOVELTY_GEN, size_t, 500, "Number of generations to lower pmin."),
  VALUE(NOVELTY_CAP, size_t, 512, "Cap on number of solutions allowed in the archive"),
  VALUE(NOVELTY_CQS, bool, false, "Do we cap solutions?"),

  GROUP(LEXICASE, "Parameters for lexicase."),
  VALUE(LEX_EPS, double, 0.0, "Parameter estimate for lexicase epsilon."),
  VALUE(LEX_DS_RATE, double, 0.5, "Percent of test cases to sample for down-sampled lexicase?"),
  VALUE(LEX_DS_MODE, std::string, "random", "How should test cases be down sampled? Options: \nnone \nrandom \nmaxmin-full \nmaxmin-pop-sample"),
  VALUE(LEX_DS_POP_RATE, double, 0.5, "Percent of population to use when computing maxmin test case sample"),

  GROUP(NONDOMINATED_SORTING, "Parameters for nondominated sorting algorithm."),
  VALUE(NDS_MAX, double, 9000000000000000.0, "Large dummy number ."),
  VALUE(NDS_RED, double, 0.99, "Default large value for different Pareto groups."),
  VALUE(NDS_ALP, double, 2.0, "Alpha value for Pareto fitness sharing."),
  VALUE(NDS_SIG, double, 0.1, "Sigma value for Pareto fitness sharing."),

  GROUP(OUTPUT, "Output rates for OpenWorld"),
  VALUE(OUTPUT_DIR, std::string, "./output/", "What directory are we dumping all this data"),
  VALUE(OUTPUT_PHYLO_DATA_INTERVAL, size_t, 10, "How often should we output data?"),
  VALUE(OUTPUT_SUMMARY_DATA_INTERVAL, size_t, 10, "How often should we output summary data?"),
  VALUE(PRINT_INTERVAL, size_t, 1, "How often do we print run status information?")
)

}