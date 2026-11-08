"""
READ ME

Dans ce fichier nous resolvons le problème pour ensuite tracer des graphes 2D.

Les modifications possibles sont au niveau des constantes, le nom des figures 
et des fichiers csv, les axes et les titres.
"""
from modules.class_fille.amortisseur import Amortisseur
from modules.class_fille.palier_lisse import PalierLisse
from modules.class_fille.traceur_csv import TraceurCSV

R, c, omega, mu = 0.020, 0.002, 500.0, 0.2
f_tau0 = (mu * omega * R) / c

cas_etudes = [(Amortisseur, 'am', 0.25, 0.0), (Amortisseur, 'am', 0.25, 5.0), 
              (Amortisseur, 'am', 0.75, 0.0), (Amortisseur, 'am', 0.75, 5.0), 
              (PalierLisse, 'pl', 0.5, 0.0), (PalierLisse, 'pl', 0.5, 10.0), 
              (PalierLisse, 'pl', 0.8, 0.0), (PalierLisse, 'pl', 0.8, 10.0)]

for Classe, nom, eps, tau in cas_etudes:
    obj = Classe(R, c, eps, omega, mu, tau * f_tau0)
    nom_complet = "Amortisseur" if nom == 'am' else "Palier Lisse"

    chemin_2d = obj.exporter_champ_vitesse_2D_csv(f'champ_vit_2D_{nom}_eps{eps*100:.0f}_tau{tau*10:.0f}.csv')

    TraceurCSV.plot_2D(nom_fichier = chemin_2d, col_x = 'theta_pi', col_y = 'y_star', col_z = 'u_star',
                       nom_pdf = f'Carte_2D_Vitesse_{nom}_eps{eps*100:.0f}_tau{tau*10:.0f}.pdf', 
                       titre = rf'{nom_complet} : Champ de vitesse ($\epsilon = {eps}$, $\tau_0^* = {tau}$)',
                       xlabel = r'$\theta / \pi$', ylabel = r'$y^*$', zlabel = r'$u^*$')