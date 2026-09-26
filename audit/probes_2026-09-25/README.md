# Robustness probes for system-identification-pa (25 Sep 2026)

Independent code: nothing is imported from the repository except in recheck_repo.py.
Run from this folder with numpy and scipy installed.

| script | what it produces |
|---|---|
| core.py | the local model, dose trains, noise model and Fisher-matrix precision |
| fit.py | MAP estimator and record simulator |
| p1_reproduce.py | the paper's headline numbers, reproduced |
| p2_combination.py | Gap A: combined routes at 5% and 10% of normal capacity |
| p2b_mc.py | Gap A: Monte Carlo check of the Fisher bound |
| p3_shape.py, p3b_free_halflife.py, p3c.py | Gap B: misspecified increment time course |
| p4_drift.py | unmodelled drifts over two years |
| p5_noise.py | Gap C: marker variability sweep |
| p6_pool.py | slow peripheral pool behind plasma (about 1 min) |
| recheck_repo.py | Gap B re-checked through the repository's own likelihood (set the repo path at the top) |
