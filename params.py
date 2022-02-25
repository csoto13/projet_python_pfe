import numpy as np

# population
N = 1

# paramètres fixe pour les simulations sur le temps et pour la prévalence
p = 0.1           # proportion de plantes primées
c1 = 0.3          # coût de virulence sur les infections primaires et secondaires
c2 = 0.3          # coût de virulence sur la production des formes de survies

# conditions initiales
Sr_0 = 0          # plantes résistantes
Is_0 = 0.01       # plantes sensibles infectées par un avirulent
Sr_star_0 = 0     # plantes primées
Js_0 = 0.01       # plantes sensibles infectées par un virulent
Jr_0 = 0.01       # plantes résistantes infectées par un virulent
xi_avir_0 = 0     # intégrale de Ss (plantes sensibles) dans une épidémie d'avirulents
sigma_0 = 0       # intégrale de Sr* (plantes primées) dans une épidémie d'avirulents
xi_vir_0 = 0      # intégrale de Ss dans une épidémie de virulents

# paramètres variant
len_pc = 20
p_, dp = np.linspace(0,1,len_pc, retstep=True)
c_ = np.linspace(0,1,len_pc)

# prévalence
len_pc_prev = 100
p_prev = np.linspace(0,1,len_pc_prev)
c_prev = np.linspace(0,1,len_pc_prev)

# critère de convergence
eps = 1e-6

# nombre d'années pour la simulation sur le temps
year = 15

# paramètres fixes du modèle
# Lambda : taux d'épuisement de l'inoculum primaire
# Theta : taux de transmission de l'inoculum primaire (infections primaires)
# mu : taux de mortalité de l'inoculum primaire
# beta : taux de transmission des agents pathogènes avirulents et virulents (infections secondaires)
# delta : perte d'infecsiosité
# gamma : perte de l'efficacité de la SAR
# rho : efficacité de la SAR (du priming)
# pi_s : taux de "conversion" des plantes infectées en inoculum primaire
# params_mod = Lambda, Theta, mu, beta, delta, gamma, rho, pi_s
params_mod_1 = 0.52, 0.4875, 0.0072, 0.03, 0, 0.01, 0.2, 1
params_mod_2 = 0.52, 0.4875, 0.0072, 0.03, 0, 0, 0.2, 1
params_mod_4 = 0.52, 0.4875, 0.0072, 0.6, 0.01, 0.2, 0.3, 1

# définition des paramètres du tspan
t0 = 0     # début de l'année
tau = 184  # durée saison estivale (été)
T = 365    # durée d'une année (d'une saison)
pas_t = 1

