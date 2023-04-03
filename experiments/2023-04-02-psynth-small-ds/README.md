# 2023-04-02 - Program synthesis small down-sample sizes

Push sampling to more extreme levels. E.g., where not all output categories can be represented in down-sample.

Problems:

- Grade
- FizzBuzz

Sample sizes

- 0.01
- 0.02
- 0.04

Depth of search:

- 5
- 10

Generations:

- 300 (for full)

Other changes:

- Change to single function
- Get rid of tag mutations (simplifies genotypes)
- Expand range on instruction arguments

```
set PRG_MIN_FUNC_CNT 1        # Minimum number of functions per program.
set PRG_MAX_FUNC_CNT 1        # Maximum number of functions per program.
set PRG_MIN_FUNC_INST_CNT 1   # Minimum number of instructions per function.
set PRG_MAX_FUNC_INST_CNT 128  # Maximum number of instructions per function.
set PRG_INST_MIN_ARG_VAL -8   # Minimum instruction-argment value
set PRG_INST_MAX_ARG_VAL 8   # Maximum instruction-argument value

set MUT_RATE_INST_ARG_SUB 0.001    # InstArgSub rate
set MUT_RATE_INST_SUB 0.001        # InstSub rate
set MUT_RATE_INST_INS 0.001        # InstIns rate
set MUT_RATE_INST_DEL 0.001        # InstDel rate
set MUT_RATE_SEQ_SLIP 0.05         # SeqSlip rate
set MUT_RATE_FUNC_DUP 0.0          # FuncDup rate
set MUT_RATE_FUNC_DEL 0.0          # FuncDel rate
set MUT_RATE_INST_TAG_BF 0.0       # InstArgTagBF rate
set MUT_RATE_FUNC_TAG_BF 0.0       # FuncTagBF rate
set MUT_RATE_INST_TAG_SINGLE_BF 0  # Per-tag single bit flip rate
set MUT_RATE_FUNC_TAG_SINGLE_BF 0  # Per-tag single bit flip rate
set MUT_RATE_INST_TAG_SEQ_RAND 0   # Per-tag sequence randomization rate
set MUT_RATE_FUNC_TAG_SEQ_RAND 0   # Per-tag sequence randomization rate
```