BEGIN RUNQPA SPECIFICATIONS
!  write-problem-data                                F
!  problem-data-file-name                            QPA.data
!  problem-data-file-device                          26
!  write-initial-sif                                 T
!  initial-sif-file-name                             INITIAL.SIF
!  initial-sif-file-device                           51
!  initial-rho-g                                     1.0D+6
!  initial-rho-b                                     1.0D+6
!  scale-problem                                     0
!  pre-solve-problem                                 F
!  write-presolved-sif                               F
!  presolved-sif-file-name                           PRESOLVE.SIF
!  presolved-sif-file-device                         52
!  solve-problem                                    
!  print-full-solution                               F
!  write-solution
!  solution-file-name                                QPASOL.d
!  solution-file-device                              62
!  write-result-summary                              
!  result-summary-file-name                          QPARES.d
!  result-summary-file-device                        47
END RUNQPA SPECIFICATIONS

BEGIN QPA SPECIFICATIONS
!  error-printout-device                             6
!  printout-device                                   6
!  print-level                                       1
!  maximum-number-of-iterations                      100000
!  start-print                                       -1
!  stop-print                                        -1
!  factorization-used                                0
!  maximum-column-nonzeros-in-schur-complement       35
!  maximum-dimension-of-schur-complement             75
!  initial-integer-workspace                         10000
!  initial-real-workspace                            10000
!  maximum-refinements                               0
!  maximum-infeasible-iterations-before-rho-increase 100
!  maximum-number-of-cg-iterations                   10
!  preconditioner-used                               1
!  semi-bandwidth-for-band-preconditioner            0
!  full-max-fill-ratio                               10
!  deletion-strategy                                 0
!  restore-problem-on-output                         2
!  residual-monitor-interval                         1
!  cold-start-strategy                               1
!  infinity-value                                    1.0D+10
!  feasibility-tolerance                             1.0D-6
!  minimum-objective-before-unbounded                -1.0D+32
!  increase-rho-g-factor                             2.0
!  increase-rho-b-factor                             2.0
!  infeasible-g-required-improvement-factor          0.75
!  infeasible-b-required-improvement-factor          0.75
!  pivot-tolerance-used                              1.0D-8
!  pivot-tolerance-used-for-dependencies             0.5
!  zero-pivot-tolerance                              1.0D-12
!  multiplier-tolerance                              1.0D-8
!  inner-iteration-relative-accuracy-required        0.0
!  inner-iteration-absolute-accuracy-required        1.0E-8
!  maximum-cpu-time-limit                            -1.0
!  treat-zero-bounds-as-general                      F
!  solve-qp                                          T
!  solve-within-bounds                               T
!  temporarily-perturb-constraint-bounds             T
!  array-syntax-worse-than-do-loop                   F
END QPA SPECIFICATIONS

BEGIN PRESOLVE SPECIFICATIONS
!  printout-device                                   6
!  error-printout-device                             6
!  print-level                                       TRACE
!  presolve-termination-strategy                     REDUCED_SIZE
!  maximum-number-of-transformations                 1000000
!  maximum-number-of-passes                          25
!  constraints-accuracy                              1.0D-6
!  dual-variables-accuracy                           1.0D-6
!  allow-dual-transformations
!  remove-redundant-variables-constraints
!  primal-constraints-analysis-frequency             1
!  dual-constraints-analysis-frequency               1
!  singleton-columns-analysis-frequency              1
!  doubleton-columns-analysis-frequency              1
!  unconstrained-variables-analysis-frequency        1
!  dependent-variables-analysis-frequency            1
!  row-sparsification-frequency                      1
!  maximum-percentage-row-fill                       -1
!  transformations-buffer-size                       50000
!  transformations-file-device                       52
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
!  maximum-growth-factor                             1.0D8
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
