# Phylogeny-informed fitness estimation

[![supplemental](https://img.shields.io/badge/go%20to-supplemental%20material-ff69b4)](https://lalejini.com/phylogeny-informed-evaluation/bookdown/book/)
[![DOI](https://zenodo.org/badge/569026105.svg)](https://zenodo.org/badge/latestdoi/569026105)
[![OSF](https://img.shields.io/badge/data%20%40%20OSF-10.17605%2FOSF.IO%2FWXCKN-blue)](https://osf.io/wxckn/)

## Overview

### Abstract

> Phylogenies (ancestry trees) depict the evolutionary history of an evolving population.
In evolutionary computing, a phylogeny can reveal how an evolutionary algorithm steers a population through a search space, illuminating the step-by-step process by which any solutions evolve.
Thus far, phylogenetic analyses have primarily been applied as post-hoc analyses used to deepen our understanding of existing evolutionary algorithms.
Here, we investigate whether phylogenetic analyses can be used at runtime to augment parent selection procedures during an evolutionary search.
Specifically, we propose phylogeny-informed fitness estimation, which exploits a population's phylogeny to estimate fitness evaluations.
We demonstrate phylogeny-informed fitness estimation in the context of the down-sampled lexicase and cohort lexicase selection algorithms on two diagnostic problems and four genetic programming (GP) problems.
Our results indicate that phylogeny-informed fitness estimation can mitigate the drawbacks of down-sampled lexicase, improving diversity maintenance and search space exploration.
However, the extent to which phylogeny-informed fitness estimation improves problem-solving success for GP varies by problem, subsampling method, and subsampling level.
This work serves as an initial step toward improving evolutionary algorithms by exploiting runtime phylogenetic analysis.

## Repository guide

- `docs/` contains supplemental documentation for our methods.
- `experiments/` contains HPC job submission scripts, configuration files, and data analyses for all experiments.
- `include/` contains C++ implementations of experiment software (header only).
- `scripts/` contains generically useful scripts used in this work.
- `source/` contains .cpp files that can be compiled to run our experiments.

