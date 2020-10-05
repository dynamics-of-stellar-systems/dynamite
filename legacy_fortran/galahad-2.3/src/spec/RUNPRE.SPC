BEGIN RUNPRE SPECIFICATIONS (default)
!  write-initial-sif
!  initial-sif-file-name                        initial.sif
!  initial-sif-file-device                      51
!  write-presolved-sif
!  presolved-sif-file-name                      presolved.sif
!  presolved-sif-file-device                    52
END RUNPRE SPECIFICATIONS

BEGIN PRESOLVE SPECIFICATIONS (DEFAULT)
!  printout-device                              6
!  error-printout-device                        6
   print-level                                  TRACE
!  presolve-termination-strategy                REDUCED_SIZE
!  maximum-number-of-transformations            1000000
!  maximum-number-of-passes                     25
!  constraints-accuracy                         1.0D-6
!  dual-variables-accuracy                      1.0D-6
!  primal-constraints-analysis-frequency        1
!  dual-constraints-analysis-frequency          1
!  singleton-columns-analysis-frequency         1
!  doubleton-columns-analysis-frequency         1
!  unconstrained-variables-analysis-frequency   1
!  dependent-variables-analysis-frequency       1
!  row-sparsification-frequency                 1
!  maximum-percentage-row-fill                  -1
!  transformations-buffer-size                  50000
!  transformations-file-device                  52
!  transformations-file-status                  KEEP
!  transformations-file-name                    transf.sav
!  primal-feasibility-check                     NONE
!  dual-feasibility-check                       NONE
!  active-multipliers-sign                      POSITIVE
!  inactive-multipliers-value                   LEAVE_AS_IS
!  active-dual-variables-sign                   POSITIVE
!  inactive-dual-variables-value                LEAVE_AS_IS
!  primal-variables-bound-status                TIGHTEST
!  dual-variables-bound-status                  TIGHTEST
!  constraints-bound-status                     TIGHTEST
!  multipliers-bound-status                     TIGHTEST
!  infinity-value                               1.0D19
!  pivoting-threshold                           1.10D-10
!  minimum-relative-bound-improvement           1.0D-10
!  compute-quadratic-value
!  compute-objective-constant
!  compute-gradient
!  compute-Hessian
!  compute-constraints-matrix
!  compute-primal-variables-values
!  compute-primal-variables-bounds
!  compute-dual-variables-values
!  compute-dual-variables-bounds
!  compute-constraints-values
!  compute-constraints-bounds
!  compute-multipliers-values
!  compute-multipliers-bounds
END PRESOLVE SPECIFICATIONS
