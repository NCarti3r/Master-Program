"""
READ ME

Dans ce fichier nous affichons des graphes 1D et 3D, les données sont extraites 
des fichiers "main_amortisseur" et "main_palier_lisse".

Il faut éviter de modifier ce qu'il y a dans ce fichier au risque de ne pas faire marcher le code,
les seules modifications autorisées sont le nom des axes et le titre.
"""
from modules.class_fille.traceur_csv import TraceurCSV

import main_amortisseur as am
import main_palier_lisse as pal

print("\nGénération des Graphiques")
labels_art_vit = [rf'Article ($\tau_0^* = {t}$)' for t in am.valeurs_star_vit]

# Profils de vitesse (Figures 3a, 3b, 3c)
cas_vit = [('3a', am.fichiers_vit_3a, 'Figure_3a_Vitesse_Amortisseur', r'Amortisseur : Vitesse ($dp^*/dx^* = -4$)'),
           ('3b', pal.fichiers_vit_3b, 'Figure_3b_Vitesse_Palier_Defavorable', r'Palier Lisse : Vitesse ($dp^*/dx^* = -4$)'),
           ('3c', pal.fichiers_vit_3c, 'Figure_3c_Vitesse_Palier_Favorable', r'Palier Lisse : Vitesse ($dp^*/dx^* = +4$)')]
    
for fig, fich_simu, nom_pdf, titre in cas_vit:
    fich_art = [f'Article/ref_vitesse_fig{fig}_tau{str(t).replace(".", "")}.csv' for t in [0, 0.25, 1]]
    TraceurCSV.comparer_csv(fichiers = fich_art + fich_simu, labels = labels_art_vit + am.labels_vit,
                            nom_pdf = nom_pdf, col_x = 'y_star', col_y = 'u_star', titre = titre,
                            xlabel = r'Coordonnée relative $y^*$', ylabel = r'Vitesse $u^*$')
        
# Pressions 2D et 3D (Figures 4, 5, 7, 8) 
cas_pression = [('am', '025', 0.25, am.valeurs_tau0_star, am.resultats[0.25], 'Amortisseur'),
                ('am', '075', 0.75, am.valeurs_tau0_star, am.resultats[0.75], 'Amortisseur'),
                ('pal', '05', 0.5, pal.valeurs_tau0_star_palier, pal.resultats_palier[0.5], 'Palier'),
                ('pal', '08', 0.8, pal.valeurs_tau0_star_palier, pal.resultats_palier[0.8], 'Palier')]
    
for pref, code, eps, liste_tau, dico_res, nom in cas_pression:
    fich_art = [f'Article/ref_{pref}_eps{code}_tau{int(t)}.csv' for t in liste_tau]
    labs_art = [rf'Article ($\tau_0^* = {t}$)' for t in liste_tau]
    TraceurCSV.comparer_csv(fichiers = fich_art + dico_res["fichiers"], labels = labs_art + dico_res["labels"],
                            nom_pdf = f'Figure_{code}_{nom}_Pression', col_x = 'theta_pi', col_y = 'P_star',
                            titre = rf'{nom} : Pression adimensionnelle ($\epsilon = {eps}$)',
                            xlabel = r'Coordonnée angulaire $\theta / \pi$', ylabel = r'Pression $P^*$')
    
    TraceurCSV.plot_3D(fichiers_csv = dico_res["fichiers"], y_values = liste_tau, 
                       x_col = 'theta_pi', z_col = 'P_star', 
                       titre = rf'Surface 3D de la pression ({nom} $\epsilon = {eps}$)', 
                       nom_pdf = f'Surface_3D_Pression_{nom}_{code}', 
                       xlabel = r'Coordonnée angulaire $\theta / \pi$', 
                       ylabel = r'Contrainte seuil $\tau_0^*$', zlabel = r'Pression $P^*$', 
                       colorbar_label = r"Pression adimensionnelle $P^*$")
        
# Frontières des noyaux (Figures 6 et 9)
TraceurCSV.tracer_frontieres_noyau(am.resultats[0.25]["fichiers"][1], 
                                   r'Amortisseur : Frontières du noyau ($\epsilon = 0.25, \tau_0^* = 1.0$)', 
                                   'Figure_6_Frontieres_Noyau')
TraceurCSV.tracer_frontieres_noyau(pal.resultats_palier[0.8]["fichiers"][1], 
                                   r'Palier Lisse : Frontières du noyau ($\epsilon = 0.8, \tau_0^* = 10$)', 
                                   'Figure_9_Frontieres_Noyau')

print("\nFin de Génération des Graphiques")