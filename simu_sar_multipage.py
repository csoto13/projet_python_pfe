import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
from multipage import save, MultiPage, start_app, clear_cache

from scipy.integrate import odeint
from params import *
import mes_fcts_SAR_mod1 as mf_mod1
import mes_fcts_SAR_mod2 as mf_mod2


start_app() 

app = MultiPage()
app.start_button = "Commencer"
app.navbar_name = "Navigation"
app.next_page_button = "Page suivante"
app.previous_page_button = "Page précédente"

def startpage():
    st.markdown("Pour utiliser l'application de simulation d'une épidémie dans un cadre saisonnier appuyer sur commencer.")

def simu(mod, y, t0, tau, T, p, Sr_0, Is_0, Sr_star_0, Js_0, Jr_0, c1, c2, params_mod):
    
    Lambda, Theta, mu, beta, delta, gamma, rho, pi_s = params_mod

    # récupère la solution
    if mod == 1:
        sol_SAR = mf_mod1.seasonal_SAR(y, t0, tau, T, p, Is_0, Sr_star_0, Js_0, Jr_0, c1, c2, params_mod)
    if mod == 2:
        sol_SAR = mf_mod2.seasonal_SAR(y, t0, tau, T, p, Sr_0, Is_0, Js_0, Jr_0, c1, c2, params_mod)

    # solutions des différents compartiments
    Pa_plot = sol_SAR[0]
    Pv_plot = sol_SAR[1]
    Ss_plot = sol_SAR[2][:,0]
    Sr_plot = sol_SAR[2][:,1]
    Is_plot = sol_SAR[2][:,2]
    Sr_star_plot = sol_SAR[2][:,3]
    Js_plot = sol_SAR[2][:,4]
    Jr_plot = sol_SAR[2][:,5]


    # définition de la prévalence
    W_plot = np.zeros(len(Is_plot))
    for i in range (len(Is_plot)):
        if Is_plot[i] == None:
            W_plot[i] = None
        else:
            W_plot[i] = Is_plot[i]+Js_plot[i]+Jr_plot[i]

    # pour plotter selon les années
    year_plot_P = (sol_SAR[3][0])/T
    year_plot_SI = (sol_SAR[3][1])/T


    # affichage des solutions
    # création d'une figure et d'un système d'axes
    fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(12,6))

    # titre de la figure
    fig1.suptitle(r"Simulation épidémie ; $p={}, c1={}, c2={}, \Lambda={}, \Theta={}, \mu={}, \beta={}, \delta={}, \gamma={}, \rho={}, \pi_s={}$".format(p, c1, c2, Lambda, Theta, mu, beta, delta, gamma, rho, pi_s), va='top', fontsize='14')


    # Plantes saines
    # tracé des simulations par rapport au temps
    ax1.plot(year_plot_SI, Ss_plot, color = 'green', label = '$S_s$')
    ax1.plot(year_plot_SI, Sr_plot, color = 'chartreuse', label='$S_r$')

    # labellisation des axes
    ax1.set_ylabel('Densité de population', fontsize='10')

    # ajout d'un titre
    ax1.set_title('Plantes saines', fontsize='14')

    # légende
    ax1.legend(loc="upper right", prop={'size':8})

    # ajout d'une grille
    ax1.grid(axis='x')


    # Avirulents
    # tracé des simulations par rapport au temps
    ax2.plot(year_plot_SI, Is_plot, color = 'orange', label = '$I_s$')
    ax2.plot(year_plot_SI, Sr_star_plot, color = 'yellow', label='$S_r^*$')
    ax2.plot(year_plot_P, Pa_plot, color = 'C4', label='$P_a$')

    # ajout d'un titre
    ax2.set_title('Avirulents', fontsize='14')

    # légende
    ax2.legend(loc="upper right", prop={'size':8})

    # ajout d'une grille
    ax2.grid(axis='x')


    # Virulents
    # tracé des simulations par rapport au temps
    ax3.plot(year_plot_SI, Js_plot, color = 'cyan', label='$J_s$')
    ax3.plot(year_plot_SI, Jr_plot, color = 'red', label = '$J_r$')
    ax3.plot(year_plot_P, Pv_plot, color = 'magenta', label = '$P_v$')

    # labellisation des axes
    ax3.set_ylabel('Densité de population', fontsize='10')
    ax3.set_xlabel('Temps (années)', fontsize='10')

    # ajout d'un titre
    ax3.set_title('Virulents', fontsize='14')

    # légende
    ax3.legend(loc="upper right", prop={'size':8})

    # ajout d'une grille
    ax3.grid(axis='x')


    # Prévalence
    # tracé des simulations par rapport au temps
    ax4.plot(year_plot_SI, W_plot, color = 'dodgerblue', label='$W$')

    # labellisation des axes
    ax4.set_xlabel('Temps (années)', fontsize='10')

    # ajout d'un titre
    ax4.set_title('Prévalence', fontsize='14')

    # légende
    ax4.legend(loc="upper right", prop={'size':8})

    # ajout d'une grille
    ax4.grid(axis='x')


    # affiche les différentes saisons (été=blanc/hiver=gris)
    for k in range(year):
        x = np.linspace((k*T + tau)/T,((k+1)*T)/T,year*T)
        y = np.ones(len(x))
        ax1.fill_between(x, y, color='0.8')
        ax2.fill_between(x, y, color='0.8')
        ax3.fill_between(x, y, color='0.8')
        ax4.fill_between(x, y, color='0.8')


    #fig1.savefig('Simulation épidémie ; p={:.2f}, c={:.2f}, Lambda={}, Theta={}, mu={}, beta={}, delta={}, gamma={}, rho={}, pi_s={}.pdf'.format(p,c,Lambda, Theta, mu, beta, delta, gamma, rho, pi_s))
    
    st.pyplot(fig1)

    
    
    
def app1(prev_vars):
    
    st.title("Simulation modèle 1 : replante plantes sensibles et résistantes")    
    
    
    p = st.slider('Proportion de plantes résistantes p :', min_value=0.0, max_value=1.0, value = 0.1, step=0.01)
    c1 = st.slider('Coût de virulence des infections primaires et secondaires c_1 :', min_value=0.0, max_value=1.0, value = 0.1, step=0.01)
    c2 = st.slider('Coût de virulence des productions des formes de survies c_2 :', min_value=0.0, max_value=1.0, value = 0.1, step=0.01)
    
    simu(1, year, t0, tau, T, p, Sr_0, Is_0, Sr_star_0, Js_0, Jr_0, c1, c2, params_mod_1)
    
   
        
def app2(prev_vars):
    
    st.title("Simulation modèle 2 : replante plantes sensibles et primées")    
            
    p = st.slider('Proportion de plantes primées p :', min_value=0.0, max_value=1.0, value = 0.1, step=0.01)
    c1 = st.slider('Coût de virulence des infections primaires et secondaires c_1 :', min_value=0.0, max_value=1.0, value = 0.1, step=0.01)
    c2 = st.slider('Coût de virulence des productions des formes de survies c_2 :', min_value=0.0, max_value=1.0, value = 0.1, step=0.01)
    
    simu(2, year, t0, tau, T, p, Sr_0, Is_0, Sr_star_0, Js_0, Jr_0, c1, c2, params_mod_1)
    
    
app.set_initial_page(startpage)
app.add_app("Modèle 1", app1)
app.add_app("Modèle 2", app2)
app.run()