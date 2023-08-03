# Phylogeny-based Estimation Lineage Distance Analysis

This section compares the phylogenetic distances of retrieved references in relative and ancestor-based modes.
We test the hypothesis that relative-based lookup retrieves closer relatives particularly under cohort (as opposed to down-sample) lexicase, due to the availability of a candidate solution's unevaluated traits among contemporary peers.
We find the hypothesis to be supported.
These analyses explore downsample rate 0.5 on the diagnostic problem set.
Other downsample rates and problem domains remain to be explored.

The more computationally-expensive relative-based lookup does retrieve significantly more closely related trait references under cohort lexicase.
The relative improvement of mean reference phylogenetic distance, compared to ancestor-based lookup, ranges between 17 and 24%.
Figure PELDA-1, below, shows distributions of phylogenetic lookup distances under cohort selection.

In comparison, down-sample lexicase hosts essentially negligible differences in relatedness between phylogeny lookup methods.
We observe effect sizes of at most 3% on mean lookup distance.
Only one surveyed difference in lookup distance is significant: under the contradictory objectives diagnostic, where ancestor-based lookup outperforms relative-based lookup by 1%.

Table PELDA-1 enumerates statistical results.
The notebook source for figures and statistics reported in this section is available [here](https://github.com/mmore500/phylogeny-informed-evaluation/blob/03f77e07e2d714ee97943aab1e8f79a7a470705d/trait-estimation-error.ipynb).

## Figure PELDA-1

Phylogenetic lookup distance distributions under cohort lexicase.

![](https://raw.githubusercontent.com/mmore500/phylogeny-informed-evaluation/03f77e07e2d714ee97943aab1e8f79a7a470705d/teeplots/EVAL_MODE%3Dcohort%2Bcol%3Ddiagnostic%2Brow%3Deval-mode%2Bviz%3Dfacet-barplot%2Bx%3Deval-fit-est-mode%2By%3Dtraits-estimation-dist%2Bext%3D.png)
*a. barplot, with bootstrapped 95% CI*

![](https://raw.githubusercontent.com/mmore500/phylogeny-informed-evaluation/03f77e07e2d714ee97943aab1e8f79a7a470705d/teeplots/EVAL_MODE%3Dcohort%2Bcol%3Ddiagnostic%2Brow%3Deval-mode%2Bshowfliers%3DFalse%2Bviz%3Dfacet-boxplot%2Bx%3Deval-fit-est-mode%2By%3Dtraits-estimation-dist%2Bext%3D.png)
*b. boxplot, with rendering of outlying values disabled*

![](https://raw.githubusercontent.com/mmore500/phylogeny-informed-evaluation/03f77e07e2d714ee97943aab1e8f79a7a470705d/teeplots/EVAL_MODE%3Dcohort%2Bcol%3Ddiagnostic%2Brow%3Deval-mode%2Bshowfliers%3DTrue%2Bviz%3Dfacet-boxplot%2Bx%3Deval-fit-est-mode%2By%3Dtraits-estimation-dist%2Bext%3D.png)
*c. boxplot, with rendering of outlying values enabled*


*Three visualizations of phylogenetic lookup distances under cohort lexicase, with downsampling rate 0.5.*
*Lookup distance zero indicates that the trait was cached within the taxon of the genome performing lookup.*

## Figure PELDA-2

Phylogenetic lookup distance distributions under down-sample lexicase.

![](https://raw.githubusercontent.com/mmore500/phylogeny-informed-evaluation/03f77e07e2d714ee97943aab1e8f79a7a470705d/teeplots/EVAL_MODE%3Ddown-sample%2Bcol%3Ddiagnostic%2Brow%3Deval-mode%2Bviz%3Dfacet-barplot%2Bx%3Deval-fit-est-mode%2By%3Dtraits-estimation-dist%2Bext%3D.png)
*a. barplot, with bootstrapped 95% CI*


![](https://raw.githubusercontent.com/mmore500/phylogeny-informed-evaluation/03f77e07e2d714ee97943aab1e8f79a7a470705d/teeplots/EVAL_MODE%3Ddown-sample%2Bcol%3Ddiagnostic%2Brow%3Deval-mode%2Bshowfliers%3DFalse%2Bviz%3Dfacet-boxplot%2Bx%3Deval-fit-est-mode%2By%3Dtraits-estimation-dist%2Bext%3D.png)
*b. boxplot, with rendering of outlying values disabled*


![](https://raw.githubusercontent.com/mmore500/phylogeny-informed-evaluation/03f77e07e2d714ee97943aab1e8f79a7a470705d/teeplots/EVAL_MODE%3Ddown-sample%2Bcol%3Ddiagnostic%2Brow%3Deval-mode%2Bshowfliers%3DTrue%2Bviz%3Dfacet-boxplot%2Bx%3Deval-fit-est-mode%2By%3Dtraits-estimation-dist%2Bext%3D.png)
*c. boxplot, with rendering of outlying values enabled*

*Three visualizations of phylogenetic lookup distances under down-sample lexicase, with downsampling rate 0.5.*
*Lookup distance zero indicates that the trait was cached within the taxon of the genome performing lookup.*

# Table PELDA-1

| lexicase    | diagnostic               | lowest-distance estimator | mean lookup distance improvement | Mann-Whitney signif. |
|-------------|--------------------------|------------------------|-------------------------------------------|-----------------|
| cohort      | contradictory objectives | relative               | 17%                           | __*p < 3e-10*__ |
| cohort      | exploitation rate        | relative               | 19%                           | __*p < 0.002*__ |
| cohort      | multipath exploration    | relative               | 24%                           | __*p < 1e-22*__ |
| down-sample | contradictory objectives | relative               | 3%                            | *p = 0.89*      |
| down-sample | exploitation rate        | relative               | 3%                            | *p = 0.09*      |
| down-sample | multipath exploration    | ancestor               | 1%                            | __*p < 0.007*__ |

Statistical comparison of phylogenetic lookup distance under relative-based versus ancestor-based queries.
Bolding according to significance threshold *alpha* = 0.05.
Reported *p* values are Bonferonni-corrected.

<!--
('cohort', 'contradictory-objectives')
MannwhitneyuResult(statistic=1276252337.0, pvalue=3.7206180707822e-11)
p value conservatively corrected for 6 comparisons 2.2323708424693202e-10
mean_ancestor=2.045695766659915 mean_relative=1.7556809184481394
mean_ancestor > mean_relative
effect size 0.1651865354144886
---

('cohort', 'exploitation-rate')
MannwhitneyuResult(statistic=1117502586.0, pvalue=0.00023161688988127547)
p value conservatively corrected for 6 comparisons 0.0013897013392876528
mean_ancestor=1.948423761526914 mean_relative=1.6431501057082452
mean_ancestor > mean_relative
effect size 0.18578561676024544
---

('cohort', 'multipath-exploration')
MannwhitneyuResult(statistic=1043169007.5, pvalue=1.4926101471916312e-23)
p value conservatively corrected for 6 comparisons 8.955660883149787e-23
mean_ancestor=1.815751896474788 mean_relative=1.469510022271715
mean_ancestor > mean_relative
effect size 0.23561722543941407
---

('down-sample', 'contradictory-objectives')
MannwhitneyuResult(statistic=60951497.0, pvalue=0.14845954910462886)
p value conservatively corrected for 6 comparisons 0.8907572946277731
mean_ancestor=1.0861538461538462 mean_relative=1.055948275862069
mean_ancestor > mean_relative
effect size 0.028605160860855187
---

('down-sample', 'exploitation-rate')
MannwhitneyuResult(statistic=50826474.0, pvalue=0.014886657281090248)
p value conservatively corrected for 6 comparisons 0.08931994368654149
mean_ancestor=1.202870813397129 mean_relative=1.1722513089005235
mean_ancestor > mean_relative
effect size 0.026120256180668426
---

('down-sample', 'multipath-exploration')
MannwhitneyuResult(statistic=44594474.0, pvalue=0.0010354724688072807)
p value conservatively corrected for 6 comparisons 0.006212834812843684
mean_ancestor=1.1633879781420764 mean_relative=1.1768
mean_relative > mean_ancestor
effect size 0.011528417097228898
---
 -->
