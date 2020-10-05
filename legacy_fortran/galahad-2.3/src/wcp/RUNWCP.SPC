BEGIN RUNWCP SPECIFICATIONS
!  write-problem-data                                F
!  problem-data-file-name                            WCP.data
!  problem-data-file-device                          26
!  write-initial-sif                                 F
!  initial-sif-file-name                             INITIAL.SIF
!  initial-sif-file-device                           51
!  no-primal-objective                               F
!  no-dual-objective                                 F
!  scale-problem                                     0
!  pre-solve-problem                                 F
!  write-presolved-sif                               F
!  presolved-sif-file-name                           PRESOLVE.SIF
!  presolved-sif-file-device                         53
!  write-scaled-sif                                  NO
!  scaled-sif-file-name                              SCALED.SIF
!  scaled-sif-file-device                            58
!  solve-problem                                     NO
!  print-full-solution                               F
!  write-solution  
!  solution-file-name                                WCPSOL.d
!  solution-file-device                              62
!  write-result-summary                              
!  result-summary-file-name                          WCPRES.d
!  result-summary-file-device                        47
END RUNWCP SPECIFICATIONS  
  
BEGIN WCP SPECIFICATIONS  
!  error-printout-device                             6
!  printout-device                                   6
!  print-level                                       1
!  maximum-number-of-iterations                      1000
!  start-print                                       -1
!  stop-print                                        -1
!  initial-point-used                                0
!  factorization-used                                0
!  maximum-column-nonzeros-in-schur-complement       35
!  initial-integer-workspace                         10000
!  initial-real-workspace                            10000
!  maximum-refinements                               1
!  maximum-poor-iterations-before-infeasible         200
!  perturbation-strategy                             2
!  restore-problem-on-output                         2
!  infinity-value                                    1.0D+10
!  primal-accuracy-required                          1.0D-5
!  dual-accuracy-required                            1.0D-5
!  complementary-slackness-accuracy-required         1.0D-5
!  initial-bound-perturbation                        -1.0
!  perturbation-small                                -1.0
!  reduce-perturbation-factor                        0.25
!  reduce-perturbation-multiplier                    0.01
!  implicit-variable-tolerance                       1.0D-5
!  mininum-initial-primal-feasibility                1.0
!  mininum-initial-dual-feasibility                  1.0
!  target-barrier-parameter                          -1.0
!  target-barrier-accept-fraction                    1.0
!  increase-barrier-parameter-by                     2.0
!  required-infeasibility-reduction                  0.99
!  pivot-tolerance-used                              1.0D-12
!  pivot-tolerance-used-for-dependencies             0.5
!  zero-pivot-tolerance                              1.0D-12
!  alpha-scaling-tolerance                           1.0D-2
!  identical-bounds-tolerance                        1.0D-15
!  insufficiently-feasible-tolerance                 0.01
!  maximum-cpu-time-limit                            -1.0
!  remove-linear-dependencies                        T
!  treat-zero-bounds-as-general                      F
!  just-find-feasible-point                          F
!  use-corrector-step                                F
!  balance-initial-complementarity                   F
!  record-x-status                                   T
!  record-c-status                                   T
END WCP SPECIFICATIONS  
  
BEGIN PRESOLVE SPECIFICATIONS  
!  printout-device                                   6
!  error-printout-device                             6
!  print-level                                       DEBUG
!  presolve-termination-strategy                     REDUCED_SIZE
!  maximum-number-of-transformations                 1000000
!  maximum-number-of-passes                          25
!  constraints-accuracy                              1.0D-6
!  dual-variables-accuracy                           1.0D-6
!  allow-dual-transformations                        F
!  remove-redundant-variables-constraints            F
!  primal-constraints-analysis-frequency             0
!  dual-constraints-analysis-frequency               0
!  singleton-columns-analysis-frequency              0
!  doubleton-columns-analysis-frequency              0
!  unconstrained-variables-analysis-frequency        0
!  dependent-variables-analysis-frequency            0
!  row-sparsification-frequency                      0
!  maximum-percentage-row-fill                       -1
!  transformations-buffer-size                       50000
!  transformations-file-device                       57
!  transformations-file-status                       KEEP
!  transformations-file-name                         transf.sav
!  primal-feasibility-check                          NONE
!  dual-feasibility-check                            NONE
!  active-multipliers-sign                           POSITIVE
!  inactive-multipliers-value                        LEAVE_AS_IS
!  active-dual-variables-sign                        POSITIVE
!  inactive-dual-variables-value                     LEAVE_AS_IS
!  primal-variables-bound-status                     TIGHTEST
!  dual-variables-bound-status                       TIGHTEST
!  constraints-bound-status                          TIGHTEST
!  multipliers-bound-status                          TIGHTEST
!  infinity-value                                    1.0D19
!  pivoting-threshold                                1.10D-10
!  minimum-relative-bound-improvement                1.0D-10
!  maximum-growth-factor                             1.0D+8
   compute-quadratic-value
   compute-objective-constant
   compute-gradient
   compute-Hessian
   compute-constraints-matrix
   compute-primal-variables-values
   compute-primal-variables-bounds
   compute-dual-variables-values
   compute-dual-variables-bounds
   compute-constraints-values
   compute-constraints-bounds
   compute-multipliers-values
   compute-multipliers-bounds
END PRESOLVE SPECIFICATIONS
