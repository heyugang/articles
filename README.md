# Code Explanation

---

This section contains the core code and code submission scripts used in the article.

## Core Code

The core code of the article, including SNP selection and breed classification:

- **R Script Submission Files for Background Execution**
    - `sbatch_Rscript_fst.sh`: Submission file for the FST method
    - `sbatch_Rscript_relf.sh`: Submission file for the Relief-F method
    - `sbatch_Rscript_mrmr.sh`: Submission file for the mrmr method

- **R Code Files to Run in Conda Environment on the Server**
    - `fst_loop_script.R`: R code for running the FST method on 5 training sets
    - `relf_loop_script.R`: R code for running the Relief-F method on the first training set
    - `mrmr_loop_script.R`: R code for running the mrmr method on the first training set
