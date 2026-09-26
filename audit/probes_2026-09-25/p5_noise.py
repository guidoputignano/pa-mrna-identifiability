"""How the nominal design responds if unexplained marker variability exceeds the 5% analytical CV."""
import numpy as np
import core
from core import *
grid = np.arange(0.0, 60 * 30.4, 30.4)
for cv in (0.05, 0.10, 0.15, 0.20, 0.30):
    core.ASSAY_CV = cv
    sig, _ = sigma_for(TH_NOM)
    sd25 = med_sd(TH_NOM, grid[:25], range(12))[0]
    n = next((k for k in range(4, 61) if med_sd(TH_NOM, grid[:k], range(12))[0] <= 0.20), None)
    print(f"marker CV {cv:4.0%}: sigma {sig:5.1f}  sd(ln theta) at 25 monthly draws {sd25:.3f}  monthly draws to 0.20: {n if n else '>60'}")
core.ASSAY_CV = 0.05
