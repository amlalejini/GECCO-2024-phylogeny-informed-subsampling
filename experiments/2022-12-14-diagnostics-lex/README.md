# 2022-12-14 Testing phylogeny-informed evaluation with selection scheme diagnostics (lexicase selection)

Exploratory experiment conditions:

- Replicates: 20
- Generations: 50,000
  - Evaluations: 2560000000
- Traits: 100
- Population size: 512
- Selection methods:
  - Lexicase
- Sampling/partitioning methods:
  - Random down-sample
    - No estimation
    - Ancestor estimation
    - Relative estimation
  - Cohort partitioning
    - No estimation
    - Ancestor estimation
    - Relative estimation
  - None
- Diagnostics
  - Exploitation
  - Multi-path exploration
  - Contradictory objectives