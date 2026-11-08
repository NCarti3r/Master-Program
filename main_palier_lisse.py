"""
READ ME

Dans ce fichier nous nous mettons dans le cas du Journal Bearing (Palier lisse en français).
Ce fichier main nous permet de lancer la résolution de ce problème.

Comment l'utiliser :
    
Il faut définir les paramètres initiaux ou on laisse ceux qui sont déterminés dans l'article de Tichy.

Ensuite nous pouvons modifier les listes valeurs_star_vit / valeurs_eps_palier / valeurs_tau0_star_palier
pour tester différentes valeurs de la vitesse, valeur d'epsilon et valeur de tau0.

!!!!! Ces valeurs sont adimentionelles !!!!!

Les modifications possibles par la suite sont les noms des fichiers csv 
ou les gradients de pression adimentionels (args = dp_star).
"""
from modules.class_fille.palier_lisse import PalierLisse

if __name__ != '__main__':
    print = lambda *args, **kwargs: None

# Valeurs tirées de l'article 
R_palier, c_jeu = 0.020, 0.002        # Rayon du palier R = 20 mm / Jeu radial c = 2 mm
omega_rot, mu_visc = 500.0, 0.2       # Vitesse angulaire = 500 s^-1 / Viscosité = 0.2 Pa.s

facteur_conversion_tau0 = (mu_visc * omega_rot * R_palier) / c_jeu
valeurs_star_vit = [0.0, 0.25, 1.0]
valeurs_eps_palier = [0.5, 0.8]
valeurs_tau0_star_palier = [0.0, 10.0]

#PALIER LISSE   
print("\nCalcul des vitesses (Palier Lisse)")
fichiers_vit_3b, fichiers_vit_3c = [], []
for tau0_star in valeurs_star_vit:
    tau0_physique = tau0_star * facteur_conversion_tau0
    pal = PalierLisse(R = R_palier, c = c_jeu, epsilon = 0.5, omega = omega_rot, mu = mu_visc, tau0 = tau0_physique)
    #Gradient défavorable
    chem_3b = pal.exporter_vitesse_locale_csv(dp_star = -4.0, tau0_star = tau0_star, nom_fichier = f"vitesse_pal_fig3b_tau{tau0_star}")
    fichiers_vit_3b.append(chem_3b)
    #Gradient favorable
    chem_3c = pal.exporter_vitesse_locale_csv(dp_star = +4.0, tau0_star = tau0_star, nom_fichier = f"vitesse_pal_fig3c_tau{tau0_star}")
    fichiers_vit_3c.append(chem_3c)

resultats_palier = {}

for eps in valeurs_eps_palier:
    eps_str = f"{int(round(eps * 10)):02d}"
    
    resultats_palier[eps] = {"fichiers": [], "labels": []}
    print(f"\nSimulation du Palier Lisse (eps = {eps})")
    
    # Simulation pour chaque valeur de tau0
    for tau0_star in valeurs_tau0_star_palier:
        tau0_physique = tau0_star * facteur_conversion_tau0
        print(f"Calcul pour tau0* = {tau0_star} (tau0 = {tau0_physique} Pa)...")
        
        palier = PalierLisse(R = R_palier, c = c_jeu, epsilon = eps, omega = omega_rot, 
                             mu = mu_visc, tau0 = tau0_physique)
        
        q_final = palier.resoudre_palier()
        print(f" -> Débit adimensionnel trouvé : q* = {q_final / (omega_rot * R_palier * c_jeu):.4f}")
        
        nom_fichier = f"palier_eps{eps_str}_tau{int(tau0_star)}.csv"
        chemin_pression = palier.exporter_pression_csv(nom_fichier = nom_fichier)
        
        resultats_palier[eps]["fichiers"].append(chemin_pression)
        resultats_palier[eps]["labels"].append(rf'$\tau_0^*$ = {tau0_star}')

    #print(f"\nFichiers pour eps={eps} :", resultats_palier[eps]["fichiers"])