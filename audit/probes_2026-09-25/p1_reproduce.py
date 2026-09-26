import numpy as np
from core import *

monthly = np.arange(0.0, 730.0, 30.4)
fortn = np.arange(0.0, 730.0, 14.0)          # on dosing days
off = np.arange(0.5, 730.0, 14.5)            # 14.5 d, offset half a day
dense = np.arange(0.0, 730.0, 2.0)
print("visits: monthly", monthly.size, "14d", fortn.size, "14.5d", off.size, "2-daily", dense.size)
print("sigma nominal", sigma_for(TH_NOM))
for name, th in (("nominal 3000", 3000.0), ("10% 300", 300.0), ("5% 150", 150.0)):
    print(f"{name:14s} monthly {med_sd(th, monthly)}")
print("14 d on dosing days (nominal)", med_sd(3000.0, fortn))
print("14.5 d offset (nominal)      ", med_sd(3000.0, off))
a = med_sd(3000.0, fortn)[0]; b = med_sd(3000.0, monthly)[0]
print("alias cost vs monthly", a / b)
print("k_e exact (prior 0.001) at 5%, monthly", med_sd(150.0, monthly, ke_prior=0.001))
print("k_e exact at nominal, monthly", med_sd(3000.0, monthly, ke_prior=0.001))
# band edges on the monthly design
from scipy.optimize import brentq
f = lambda lt: med_sd(np.exp(lt), monthly)[0] - 0.20
lo = np.exp(brentq(f, np.log(1000), np.log(3000)))
hi = np.exp(brentq(f, np.log(6000), np.log(60000)))
print(f"band: {lo:.0f} to {hi:.0f} = {100*lo/3000:.0f}% to {100*hi/3000:.0f}% of nominal; disease factor {lo/150:.1f}")
