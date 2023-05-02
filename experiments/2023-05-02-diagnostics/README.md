# 2023-05-02 - Diagnostics experiments for 2023 GPTP

Parameters

- Diagnostics: exploitation-rate, exploration, contradictory-objectives
- Subsampling rates: 0.01, 0.05, 0.1, 0.5, 1
- Evaluation modes: down-sample, cohort, full
- _Generations_ fixed across runs (measuring loss relative to full)
  - Why not evaluations?
    - Want to run full for enough generations to max out (see original diagnostics paper)
    - To hold evaluations constant, we'd end up running 1% for 500k - 1000k generations, which would run extremely slowly for the contradictory objectives diagnostic (result in huge phylogenies because of deep branches)
    - We're asking about mitigating drawback. Not comparing across
- 10 replicates each treatment