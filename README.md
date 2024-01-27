# Phylogeny-informed subsampling

<!-- TODO - Badges! -->
<!-- [![supplemental](https://img.shields.io/badge/go%20to-supplemental%20material-ff69b4)](https://lalejini.com/GPTP-2023-phylogeny-informed-evaluation/bookdown/book/) -->
<!-- [![DOI](https://zenodo.org/badge/569026105.svg)](https://zenodo.org/badge/latestdoi/569026105) -->
<!-- [![OSF](https://img.shields.io/badge/data%20%40%20OSF-10.17605%2FOSF.IO%2FWXCKN-blue)](https://osf.io/wxckn/) -->

## Overview

### Abstract

> Phylogenies (ancestry trees) tell the evolutionary history of an evolving population.
  In evolutionary computing, phylogenies reveal how evolutionary algorithms steer populations through a search space by illuminating the step-by-step evolution of solutions.
  To date, phylogenetic analyses have almost exclusively been applied in post-hoc analyses of evolutionary algorithms for performance tuning and research.
  Here, we apply phylogenetic information at runtime to augment parent selection procedures that use training sets to assess candidate solution quality.
  We propose phylogeny-informed fitness estimation, thinning a fraction of costly training case evaluations by substituting the fitness profiles of near relatives as a heuristic estimate.
  We evaluate phylogeny-informed fitness estimation in the context of the down-sampled lexicase and cohort lexicase selection algorithms on two diagnostic analyses and four genetic programming (GP) problems.
  Our results indicate that phylogeny-informed fitness estimation can mitigate the drawbacks of down-sampled lexicase, improving diversity maintenance and search space exploration.
  However, the extent to which phylogeny-informed fitness estimation improves problem-solving success for GP varies by problem, subsampling method, and subsampling level.
  This work serves as an initial step toward improving evolutionary algorithms by exploiting runtime phylogenetic analysis.

## Repository guide

- `docs/` contains supplemental documentation for our methods.
- `experiments/` contains HPC job submission scripts, configuration files, and data analyses for all experiments.
- `include/` contains C++ implementations of experiment software (header only).
- `scripts/` contains generically useful scripts used in this work.
- `source/` contains .cpp files that can be compiled to run our experiments.

