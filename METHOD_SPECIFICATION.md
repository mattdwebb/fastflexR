# Method

Historical state-group pooled FLEX: cluster and time effects, cluster-X and time-X slopes, treated cluster-time cell intercepts and X slopes. Target uses weighted cell masses and historical unweighted cell covariate means. Redundant post-cell regressors allowed only if the target is estimable. Cohort grouping unsupported.

At fixed sample, assignment and design the scalar is linear in Y, with fixed
least-squares and aggregation operators. Tests check additivity and scaling;
designs/contrasts must still be rebuilt after deletion or reassignment.
TWFE uses weighted cyclic FE absorption without giant FE dummy matrices.
Other methods currently use bounded transparent SVD designs (12 million estimated
entries), failing explicitly above the cap. No shared runtime package dependency.
