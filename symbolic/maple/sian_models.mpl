# =====================================================================
#  SIAN input for the propionate models
#  Requires SIAN: https://github.com/pogudingleb/SIAN
#  Place IdentifiabilityODE.mpl in the same directory.
# =====================================================================

read "IdentifiabilityODE.mpl":

# ---- A1: baseline, time-varying input --------------------------------
printf("\n===== A1: baseline, input s(t) =====\n");
A1 := [
  diff(m(t), t) = s(t) - theta*m(t)/(Km + m(t)) - ke*m(t),
  y(t) = m(t)
]:
IdentifiabilityODE(A1, GetParameters(A1), p_mod = 0.99);

# ---- A2: constant input (steady-state patient) -----------------------
printf("\n===== A2: constant input =====\n");
A2 := [
  diff(m(t), t)  = s0(t) - theta*m(t)/(Km + m(t)) - ke*m(t),
  diff(s0(t), t) = 0,
  y(t) = m(t)
]:
IdentifiabilityODE(A2, GetParameters(A2), p_mod = 0.99);

# ---- A4: null-genotype anchor ----------------------------------------
printf("\n===== A4: null anchor, theta_2 = 0 =====\n");
A4 := [
  diff(m1(t), t) = s1(t) - th1*m1(t)/(Km + m1(t)) - ke*m1(t),
  diff(m2(t), t) = s2(t) - ke*m2(t),
  diff(s1(t), t) = 0,
  diff(s2(t), t) = 0,
  y1(t) = m1(t),
  y2(t) = m2(t)
]:
IdentifiabilityODE(A4, GetParameters(A4), p_mod = 0.99);

# ---- B: both channels saturating -------------------------------------
printf("\n===== B: both saturating =====\n");
B := [
  diff(m(t), t) = s(t) - theta*m(t)/(Km + m(t)) - ke*m(t)/(Ke + m(t)),
  y(t) = m(t)
]:
IdentifiabilityODE(B, GetParameters(B), p_mod = 0.99);

# ---- C: carnitine depletion, c unobserved ----------------------------
printf("\n===== C: carnitine depletion =====\n");
C := [
  diff(m(t), t) = s(t) - theta*m(t)/(Km + m(t)) - k*c(t)*m(t),
  diff(c(t), t) = R - k*c(t)*m(t),
  y(t) = m(t)
]:
IdentifiabilityODE(C, GetParameters(C), p_mod = 0.99);
