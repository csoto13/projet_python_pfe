import numpy as np
from scipy.integrate import odeint
from params import *


# définition du tspan
def tspan(t0, tf):
    return np.arange(t0, tf, pas_t)

# modèle été
def SAR_sum(etat, t, params):
    S_s, S_r, I_s, S_r_star, J_s, J_r = etat
    beta, delta, gamma, rho, c1 = params
    S_s_dot = - beta * I_s * S_s - (1 - c1) * beta * (J_s + J_r) * S_s
    S_r_dot = - beta * I_s * S_r - (1 - c1) * beta * (J_s + J_r) * S_r + gamma * S_r_star
    I_s_dot = beta * I_s * S_s - delta * I_s
    S_r_star_dot = beta * I_s * S_r - (1 - rho) * (1 - c1) * beta * (J_s + J_r) * S_r_star - gamma * S_r_star
    J_s_dot = (1 - c1) * beta * (J_s + J_r) * S_s - delta * J_s
    J_r_dot = (1 - c1) * beta * (J_s + J_r) * S_r + (1 - rho) * (1 - c1) * beta * (J_s + J_r) * S_r_star - delta * J_r
    etatdot = [S_s_dot, S_r_dot, I_s_dot, S_r_star_dot, J_s_dot, J_r_dot]
    return etatdot

# modèle hiver
def SAR_wint(etat, t, params):
    P_a, P_v, S_s, S_r, I_s, S_r_star, J_s, J_r = etat
    mu = params
    P_a_dot = - mu * P_a
    P_v_dot = - mu * P_v
    S_s_dot = 0
    S_r_dot = 0
    I_s_dot = 0
    S_r_star_dot = 0
    J_s_dot = 0
    J_r_dot = 0
    etatdot = [P_a_dot, P_v_dot, S_s_dot, S_r_dot, I_s_dot, S_r_star_dot, J_s_dot, J_r_dot]
    return etatdot


# modèle saisonnier après y année(s)
def seasonal_SAR(y, t0, tau, T, p, Sr_0, Is_0, Js_0, Jr_0, c1, c2, params_mod):
    """Renvoit les solutions du modèle après y année(s)"""
    """y          : nombre d'année à simuler
       t0, tau, T : temps initial, durée été, durée année
       _0         : conditions initiales du modèle
       p          : proportion de plantes primées
       c1, c2     : coût de virulence du pathogène virulent sur les infections primaires et secondaire, sur les productions des formes de survies
       eps        : critère d'arrêt 
       params_mod : Lambda, Theta, mu, beta, delta, gamma, rho, pi_s"""

    # paramètres du modèle
    Lambda, Theta, mu, beta, delta, gamma, rho, pi_s = params_mod

    # paramètres des modèles été et hiver
    params_sum = np.array([beta, delta, gamma, rho, c1])
    params_wint = np.array([mu])

    # replantation 
    Ss_0 = (1 - p) * N
    Sr_star_0 = p * N

    P_a = []
    P_v = []
    etat_saison_SI = np.empty((0, 6))
    tspan_saison_SI = []
    tspan_saison_P = []

    # initialisation des compartiments
    Ss_init = Ss_0
    Sr_init = Sr_0
    Is_init = Is_0
    Sr_star_init = Sr_star_0
    Js_init = Js_0
    Jr_init = Jr_0

    for k in range(y):
        # saison été
        etat0_sum = np.array([Ss_init, Sr_init, Is_init, Sr_star_init, Js_init, Jr_init])
        tspan_sum = tspan(t0 + k * T, k * T + tau)
        int_SAR_sum = odeint(SAR_sum, etat0_sum, tspan_sum, args=(params_sum,), hmax=pas_t)

        # fin été, début hiver
        Pa_saut = pi_s * int_SAR_sum[:, 2][-1]
        Pv_saut = (1 - c2) * pi_s * (int_SAR_sum[:, 4][-1] + int_SAR_sum[:, 5][-1])

        etat0_wint = np.array([Pa_saut, Pv_saut, 0, 0, 0, 0, 0, 0])
        tspan_wint = tspan(k * T + tau, (k + 1) * T)

        # regroupement été + hiver pour l'inoculum primaire (P)
        int_SAR_wint = odeint(SAR_wint, etat0_wint, tspan_wint, args=(params_wint,), hmax=pas_t)
        P_a = np.hstack([P_a, None, int_SAR_wint[:, 0]])
        P_v = np.hstack([P_v, None, int_SAR_wint[:, 1]])
        tspan_saison_P = np.hstack([tspan_saison_P, tspan(k * T + tau, (k + 1) * T + 1)])

        # regroupement été + hiver pour les plantes saines, primées et infectées
        etat_saison_SI = np.vstack([etat_saison_SI, int_SAR_sum[:, 0:6], [None, None, None, None, None, None]])
        tspan_saison_SI = np.hstack([tspan_saison_SI, tspan(t0 + k * T, k * T + tau + 1)])

        tspan_saison = np.array([tspan_saison_P, tspan_saison_SI])

        # conditions initiales pour l'année suivante
        Pa_0 = Pa_saut * np.exp(-mu * (T - tau))
        Pv_0 = Pv_saut * np.exp(-mu * (T - tau))
        Ss_init = Ss_0 * np.exp(-(1 / Lambda) * (Theta * Pa_0 + Theta * (1 - c1) * Pv_0))
        Sr_init = 0
        Is_init = ((Theta * Pa_0 * Ss_0) / (Theta * Pa_0 + Theta * (1 - c1) * Pv_0)) * (
                    1 - np.exp(-(1 / Lambda) * (Theta * Pa_0 + Theta * (1 - c1) * Pv_0)))
        Sr_star_init = Sr_star_0 * np.exp(- (1 / Lambda) * (1 - rho) * Theta * (1 - c1) * Pv_0)
        Js_init = ((Theta * (1 - c1) * Pv_0 * Ss_0) / (Theta * Pa_0 + Theta * (1 - c1) * Pv_0)) * (
                    1 - np.exp(-(1 / Lambda) * (Theta * Pa_0 + Theta * (1 - c1) * Pv_0)))
        Jr_init = Sr_star_0 * (1 - np.exp(- (1 / Lambda) * (1 - rho) * Theta * (1 - c1) * Pv_0))

    return P_a, P_v, etat_saison_SI, tspan_saison