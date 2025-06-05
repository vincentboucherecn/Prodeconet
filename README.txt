REPLICATION PACKAGE FOR: THE INTERACTIONS BETWEEN PRODUCTION AND ECOLOGICAL NETWORKS - V. Boucher

---
ATN.zip : This file contains the code and the output of the simulated data effectively used in the paper
The file is too big to be stored on GitHub and is available here:
---
To recplcate the analysis of Section 6:
 - Step 1: run gen_data.R - produces baseline ecosystems (without fishing)
 - Step 2: run comparison0.R - simulates the economy with isolated sectors
 - Step 3: run comparison_cal.R - calibrates TFP for production network economies
 - Step 4: run comparison.R - simulates the economy with production network
 - Step 5: run spagplot.R - computes volatility indices
 - Step 6: run analysis.R - Produces the analysis of Section 6.

 - function.R - (NOT RUN) functions definitions
 - parameters.R - (NOT RUN) parameters definitions
 - checkreasons.R (OTHER) - For information only, produces summary of rejected ecosystems
---

---
example non-existence.R : Produces Figure C.1
---

---
Stylized_example.R: Produces Figure 7
---
