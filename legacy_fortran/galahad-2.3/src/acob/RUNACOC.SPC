BEGIN RUNACOB SPECIFICATIONS
! write-problem-data                              no
! problem-data-file-name                          ACOB.data
! problem-data-file-device                        26
! write-solution                                  yes
! solution-file-name                              ACOBSOL.d
! solution-file-device                            62
! write-result-summary                            
! result-summary-file-name                        ACOBRES.d
! result-summary-file-device                      47
! check-all-derivatives                           no
! check-derivatives                               
! check-element-derivatives                       
! check-group-derivatives                         
! ignore-derivative-bugs                          no
! ignore-element-derivative-bugs                  no
! ignore-group-derivative-bugs                    no
! get-scaling-factors                             no
! scaling-print-level                             1
! use-scaling-factors                             no
! use-constraint-scaling-factors                  no
! use-variable-scaling-factors                    no
! maximizer-sought                                no
! restart-from-previous-point                     no
! restart-data-file-name                          ACOBSAVE.d
! restart-data-file-device                        59
! save-data-for-restart--every                    0
! one-norm-penalty                                T
END RUNACOB SPECIFICATIONS

BEGIN ACOB SPECIFICATIONS
! error-printout-device                           6
! printout-device                                 6
! alive-device                                    
  print-level                                     1
  maximum-number-of-iterations                    2000
  maximum-number-of-minor-iterations              1
! start-print                                     55
! stop-print                                      200
! iterations-between-printing                     1
! gradient-accuracy-required                      1.0D-5
  initial-regularisation-weight                   1.0D+0
! successful-iteration-tolerence                  0.01
! very-successful-iteration-tolerence             0.9
! regularisation-weight-decrease-factor           0.5
! regularisation-weight-increase-factor           2.0
  print-full-solution                             no
! space-critical                                  no
! deallocate-error-fatal                          no
! alive-filename                                  ALIVE.d
END ACOB SPECIFICATIONS

BEGIN GLRT SPECIFICATIONS
! error-printout-device                           6
! printout-device                                 6
  print-level                                     0
  maximum-number-of-iterations                    25
! tri-diagonal-solve-frequency                    1
  stopping-rule                                   1
! relative-accuracy-required                      1.0E-8 
! absolute-accuracy-required                      0.0
! fraction-optimality-required                    1.0
! constant-term-in-objective                      0.0
! zero-gradient-tolerance                         2.0E-15
! two-norm-regularisation                         T
! space-critical                                  F
! deallocate-error-fatal                          F
! output-line-prefix                              ""
END GLRT SPECIFICATIONS

BEGIN SBLS SPECIFICATIONS
! error-printout-device                           6
! printout-device                                 6
  print-level                                     0
! initial-workspace-for-unsymmetric-solver        0
! initial-integer-workspace                       1000
! initial-real-workspace                          1000
! maximum-refinements                             1
  preconditioner-used                             4
! semi-bandwidth-for-band-preconditioner          5
  factorization-used                              2
! maximum-column-nonzeros-in-schur-complement     35
! ordering-used                                   3
! scaling-used                                    0
! has-a-changed                                   2
! has-h-changed                                   2
! has-c-changed                                   2
! minimum-diagonal                                1.0D-5
! pivot-tolerance-used                            1.0D-12
! pivot-tolerance-used-for-basis                  0.5
! zero-pivot-tolerance                            1.0D-12
! static-pivoting-diagonal-perturbation           0.0D+0
! level-at-which-to-switch-to-static              0.0D+0
! find-basis-by-transpose                         F
! check-for-reliable-basis                        T
! remove-linear-dependencies                      T
! perturb-to-make-+ve-definite                    T
! space-critical                                  F
! deallocate-error-fatal                          F
END SBLS SPECIFICATIONS
