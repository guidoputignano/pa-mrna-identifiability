import sys, numpy as np
sys.path.insert(0,'src')
import model_C.run as C, design.run as D
from common.models import clearance_C_quasistatic, clearance_A
PC=C.PC; PA=D.P
peak,mp = C.clearance_peak(PC)
print(f"Model C clearance peak = {peak:.3f} at m = {mp:.3f}; S_BASE = {C.S_BASE}")
print(f"  largest clearable steady load = {peak/C.S_BASE:.4f}x baseline")
eq = C.equilibria_nonmonotone(clearance_C_quasistatic, PC, C.S_BASE, np.linspace(1e-3,3000.0,120000))[0]
print(f"  equilibria at baseline load: {[round(float(x),4) for x in eq]}")
mbar=float(eq[0]); sep=float(eq[1]) if len(eq)>1 else None
print(f"  attractor {mbar:.4f}, separatrix {sep:.4f} = {sep/mbar:.2f} m_bar")
print(f"  clearance peak sits at {mp/mbar:.2f} m_bar -> fraction of log-band [1,8] m_bar beyond the peak: "
      f"{max(0,(np.log(8)-np.log(mp/mbar))/np.log(8)):.3f}")
print()
for k in (1,8,14.6,20):
    m=k*D.M_BAR; s=clearance_A(m,PA); s0=clearance_A(D.M_BAR,PA)
    print(f"Model A: steady load to hold {k:5.1f} m_bar (m={m:6.1f}) = {s:9.2f} = {s/s0:6.2f}x baseline"
          f"   [over baseline: {s-s0:8.2f} = {(s-s0)/(D.LOAD_FRAC*s0):6.1f} load-bias prior sd]")
print()
print(f"Ratio of the two models' maximum clearable load: Model A needs {clearance_A(8*D.M_BAR,PA)/clearance_A(D.M_BAR,PA):.2f}x "
      f"to hold 8 m_bar; Model C's ceiling is {peak/C.S_BASE:.3f}x, a factor of "
      f"{(clearance_A(8*D.M_BAR,PA)/clearance_A(D.M_BAR,PA))/(peak/C.S_BASE):.2f}")
