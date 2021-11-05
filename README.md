Algorithm speed tests:


linear_solver_type:

ceres::DENSE_QR = 18.218, 17.143, 17.068 s

ceres::DENSE_NORMAL_CHOLESKY = 16.243, 17.080, 17.070 s

ceres::SPARSE_SCHUR = 17.225, 17.067, 17.074 s

ceres::ITERATIVE_SCHUR = 17.189, 17.082, 17.099 s (some failures)

ceres::CGNR = A looooot of failures

ceres::DENSE_SCHUR = 16.074, 16.053, 17.086 s

ceres::SPARSE_NORMAL_CHOLESKY = 16.062, 17.090, 17.065 s




trust_region_strategy_type:

ceres::DOGLEG = 17.158, 17.061, 18.107 s


I'm pretty sure LM is the default (I didn't change the trust region type for the linear solver type testing)




Profile-guided optimisation (PGO) is giving reasonable speed improvements, 109s vs 119s for one run.


options.minimizer_type = ceres::LINE_SEARCH improves speed hugely (3.5x) but reduced chisq doubles.
