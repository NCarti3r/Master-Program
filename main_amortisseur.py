"""
READ ME

Dans ce fichier nous nous mettons dans le cas du Squeeze Film Damper (Amortisseur en français).
Ce fichier main nous permet de lancer la résolution de ce problème.

Comment l'utiliser :
    
Il faut définir les paramètres initiaux ou on laisse ceux qui sont déterminés dans l'article de Tichy.

Ensuite nous pouvons modifier les listes valeurs_star_vit / eps_amortisseur / valeurs_tau0_star
pour tester différentes valeurs de la vitesse, valeur d'epsilon et valeur de tau0.

!!!!! Ces valeurs sont adimentionelles !!!!!

Les modifications possibles par la suite sont les noms des fichiers csv 
ou les gradients de pression adimentionels (args = dp_star).
"""

from modules.class_fille.amortisseur import Amortisseur

if __name__ != '__main__':
    print = lambda *args, **kwargs: None

# Valeurs tirées de l'article 
R_palier, c_jeu = 0.020, 0.002        # Rayon du palier R = 20 mm / Jeu radial c = 2 mm
omega_rot, mu_visc = 500.0, 0.2       # Vitesse angulaire = 500 s^-1 / Viscosité = 0.2 Pa.s

facteur_conversion_tau0 = (mu_visc * omega_rot * R_palier) / c_jeu           
valeurs_star_vit = [0.0, 0.25, 1.0]
eps_amortisseur = [0.25, 0.75]
valeurs_tau0_star = [0.0, 1.0, 5.0]

print("\nCalcul des vitesses (Amortisseur)")
fichiers_vit_3a = []
labels_vit = []
for tau0_star in valeurs_star_vit:
    tau0_physique = tau0_star * facteur_conversion_tau0
    amort = Amortisseur(R = R_palier, c = c_jeu, epsilon = 0.25, omega = omega_rot, mu = mu_visc, tau0 = tau0_physique)
    # Fig. 3a : gradient imposé dp*/dx* = -4 (profil local, pas de solution du palier)
    chemin = amort.exporter_vitesse_locale_csv(dp_star = -4.0, tau0_star = tau0_star, 
                                               nom_fichier = f"vitesse_am_fig3a_tau{tau0_star}")
    fichiers_vit_3a.append(chemin)
    labels_vit.append(rf'Modèle ($\tau_0^* = {tau0_star}$)')
    
resultats = {}

for eps in eps_amortisseur:
    eps_str = f"{int(round(eps * 100)):03d}"
    resultats[eps] = {"fichiers": [], "labels": []}
    print(f"\nSimulation de l'Amortisseur (eps = {eps})")
    
    # Simulation pour chaque valeur de tau0
    for tau0_star in valeurs_tau0_star:
        tau0_physique = tau0_star * facteur_conversion_tau0
        
        amortisseur = Amortisseur(R = R_palier, c = c_jeu, epsilon = eps, omega = omega_rot, 
                                  mu = mu_visc, tau0 = tau0_physique)
        q_final = amortisseur.resoudre_palier()
        print(f" -> Débit adimensionnel trouvé : q* = {q_final / (omega_rot * R_palier * c_jeu):.4f}")
        
        nom_fichier_pression = f"amortisseur_eps{eps_str}_tau{int(tau0_star)}.csv"
        chemin_pression = amortisseur.exporter_pression_csv(nom_fichier = nom_fichier_pression)
        
        resultats[eps]["fichiers"].append(chemin_pression)
        resultats[eps]["labels"].append(rf'$\tau_0^*$ = {tau0_star}')

    #print(f"\nFichiers pour eps = {eps} :", resultats[eps]["fichiers"])