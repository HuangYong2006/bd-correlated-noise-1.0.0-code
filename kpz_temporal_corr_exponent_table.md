# Analytic exponent branches for KPZ with temporally correlated noise (d = 1)

Noise: D(omega,k) = D0 + D_theta * k^(-2 rho) * omega^(-2 theta); here rho = 0 (pure temporal).
Conventions: alpha = chi (global roughness), beta = chi/z (growth), z (dynamic).
Values below are COMPUTED from the formulas quoted in the sources (see the main report for URLs).

| theta | SR-KPZ (theta <= 1/6) | Medina 1-loop: chi=(1+4t)/(3+2t), z=5/(3+2t) | Ma&Ma Flory d=1: chi=(1+2t)/(2+t), z=3/(2+t) | Fedorenko lower bound z >= 5/(3+2t) |
|-------|------------------------|-----------------------------------------------|------------------------------------------------|--------------------------------------|
| 0.1667 (1/6) | chi=0.5, beta=0.3333, z=1.5 | chi=0.5, beta=0.3333, z=1.5 | chi=0.6154, beta=0.4444, z=1.3846 | z >= 1.5 |
| 0.20 | chi=0.5, beta=0.3333, z=1.5 | chi=0.5294, beta=0.36, z=1.4706 | chi=0.6364, beta=0.4667, z=1.3636 | z >= 1.4706 |
| 0.25 | chi=0.5, beta=0.3333, z=1.5 | chi=0.5714, beta=0.40, z=1.4286 | chi=0.6667, beta=0.50, z=1.3333 | z >= 1.4286 |
| 0.30 | chi=0.5, beta=0.3333, z=1.5 | chi=0.6111, beta=0.44, z=1.3889 | chi=0.6957, beta=0.5333, z=1.3043 | z >= 1.3889 |
| 0.35 | chi=0.5, beta=0.3333, z=1.5 | chi=0.6486, beta=0.48, z=1.3514 | chi=0.7234, beta=0.5667, z=1.2766 | z >= 1.3514 |
| 0.40 | chi=0.5, beta=0.3333, z=1.5 | chi=0.6842, beta=0.52, z=1.3158 | chi=0.75,   beta=0.60,   z=1.25   | z >= 1.3158 |
| 0.45 | chi=0.5, beta=0.3333, z=1.5 | chi=0.7179, beta=0.56, z=1.2821 | chi=0.7755, beta=0.6333, z=1.2245 | z >= 1.2821 |
| 0.50 | chi=0.5, beta=0.3333, z=1.5 | chi=0.75,   beta=0.60, z=1.25   | chi=0.80,   beta=0.6667, z=1.2   | z >= 1.25 |

Checks (exact algebra, done numerically):
* LR hyperscaling relation z(1+2 theta) - 2 chi - d = 0 (d=1) is satisfied IDENTICALLY by the Medina
  one-loop branch (residual 0 at every theta), because 2 - chi = 5/(3+2 theta).
* The Ma&Ma (Flory/SA) branch VIOLATES that relation for every theta < 1 (residual e.g. -0.385 at
  theta = 1/6, -0.20 at theta = 1/2) and its z lies BELOW the Fedorenko lower bound 5/(3+2 theta)
  for all 0 < theta < 1 (they cross only at theta = 1).
* Medina 1-loop z IS exactly Fedorenko's lower bound: z = 5/(3+2 theta) = 2 - chi.
* The one-loop chi = (1+4 theta)/(3+2 theta) equals the SR value 1/2 exactly at theta = 1/6.
