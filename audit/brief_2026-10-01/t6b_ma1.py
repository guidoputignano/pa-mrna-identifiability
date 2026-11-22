"""
Task 6, the MA(1) claim checked where it applies: the EQUATION ERROR of the
finite-difference regression, at the TRUE lumped coefficients, against the
fitted OLS residual. Consecutive differences share a sample, so the equation
error should have lag-1 autocorrelation -0.5; the fitted residual need not,
because the least-squares fit can absorb part of the shared noise through the
coefficient on -psi.
"""
import numpy as np

import t6_estimators as T


def main():
    print("  lag-1 autocorrelation of the FD regression's error, median of 300 records")
    for nm in ("A", "B"):
        b = T.bench(nm)
        lt = T.lumped_true(b)
        rng = np.random.default_rng(3)
        r_true, r_clean, r_fit = [], [], []
        for _ in range(300):
            y = b["m"] * np.exp(rng.normal(0, b["cv"], b["m"].size))
            tau = np.diff(b["t"])
            psi = np.diff(y) / tau - T.u_bar(b["ex"], b["t"])
            for mm, store in ((y[:-1], r_true), (b["m"][:-1], r_clean)):
                X = np.stack([np.ones_like(mm), mm, -mm ** 2, -psi], axis=1)
                e = psi * mm - X @ lt
                store.append(np.corrcoef(e[:-1], e[1:])[0, 1])
            r_fit.append(T.est_fd(y, b["t"], b["ex"])["r1"])
        print(f"    bench {nm}: equation error at the true coefficients {np.median(r_true):+.3f} "
              f"(noiseless m in the regressors {np.median(r_clean):+.3f}); fitted OLS residual "
              f"{np.median(r_fit):+.3f}")


if __name__ == "__main__":
    main()
