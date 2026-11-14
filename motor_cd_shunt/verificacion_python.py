"""
Verificación independiente (Python) del motor CD Shunt con control por campo.
Reproduce la Parte A y simula las EDO con solve_ivp para contrastar
con los resultados de Scilab/Xcos.   Requiere: numpy, scipy
    python3 verificacion_python.py
"""
import numpy as np
from scipy.integrate import solve_ivp

Va, Ra, Rf, IL, La, J, B = 120, 0.5, 120, 20, 0.5, 0.5, 0.075
n1, n2 = 1600, 1800
w1, w2 = n1 * 2 * np.pi / 60, n2 * 2 * np.pi / 60

# ---- Parte A.1
If = Va / Rf; Ia = IL - If; Ea = Va - Ia * Ra
Ke = Ea / w1; Kt = Ke; T = Kt * Ia
print(f"A.1  If={If:.4f} A  Ia={Ia:.4f} A  Ea={Ea:.4f} V  Ke=Kt={Ke:.6f}  T={T:.4f} N.m")

# ---- Parte A.2 : w2*K^2 - Va*K + T*Ra = 0
r = np.roots([w2, -Va, T * Ra])
K2 = r.max(); Ia2 = T / K2; If2 = If * K2 / Ke; Rf2 = Va / If2
print(f"A.2  raices={r}  Ke2={K2:.6f}  Ia2={Ia2:.4f} A  If2={If2:.5f} A  Rf'={Rf2:.4f} Ohm")

# ---- Parte A.3 + B : estado estable y simulación
def ss(K, TL):
    w = (K * Va - Ra * TL) / (K * K + Ra * B)
    return w, (Va - K * w) / Ra

for K, lab in [(Ke, "Rf=120  "), (K2, "Rf=136.6")]:
    for TL in (0, 5):
        w, i = ss(K, TL)
        sol = solve_ivp(lambda t, y: [(Va - Ra * y[0] - K * y[1]) / La,
                                      (K * y[0] - B * y[1] - TL) / J],
                        [0, 30], [0, 0], max_step=1e-3, rtol=1e-9, atol=1e-11)
        k = sol.y[0].argmax()
        print(f"{lab} TL={TL}: teor n={w*30/np.pi:8.2f} rpm ia={i:7.4f} A | "
              f"sim(30 s) n={sol.y[1][-1]*30/np.pi:8.2f} ia={sol.y[0][-1]:7.4f} | "
              f"ia_max={sol.y[0][k]:.2f} A @ {sol.t[k]:.2f} s  n_max={sol.y[1].max()*30/np.pi:.1f} rpm")
