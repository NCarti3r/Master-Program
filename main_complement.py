"""
READ ME

Ce fichier permet d'observer la variation du flux, l'évolution de la vitesse max 
et les coefficients de frottements aux parois.

Les modifications possibles sont au niveau des constantes, le nom des figures, 
des fichiers csv, les axes et les titres et les méthodes utilisées Amortisseur ou Palier Lisse.
"""
from modules.class_fille.amortisseur import Amortisseur
from modules.class_fille.palier_lisse import PalierLisse
from modules.class_fille.traceur_csv import TraceurCSV

# Valeurs tirées de l'article 
R_palier, c_jeu = 0.020, 0.002        # Rayon du palier R = 20 mm / Jeu radial c = 2 mm
omega_rot, mu_visc = 500.0, 0.2       # Vitesse angulaire = 500 s^-1 / Viscosité = 0.2 Pa.s

# Contraintes pariétales

palier = PalierLisse(R = R_palier, c = c_jeu, epsilon = 0.25, omega = omega_rot, mu = mu_visc, tau0 = 2000)

palier.exporter_pression_csv('simu_complementaire.csv',
                             colonnes = ['theta_pi', 'h', 'dp_dx', 'P', 'P_star', 
                                       'ha_star', 'hb_star', 'tau_1', 'tau_2', 
                                       'tau_1_star', 'tau_2_star'])

visu = TraceurCSV('simu_complementaire.csv')
visu.tracer(col_x = 'theta_pi',col_y = ['tau_1_star', 'tau_2_star', 'P_star'],
            secondary_y = ['P_star'],  nom_pdf = 'contraintes_parietales.pdf',
            titre = 'Contraintes pariétales adimensionnelles', xlabel = r'$\theta / \pi$',
            ylabel = r'Contraintes $\tau^*$', ylabel_droit = r'Pression $P^*$')

# Conservation du Flux
amort = Amortisseur(R = 0.02, c = 0.002, epsilon = 0.5, omega = 500, mu = 0.2, tau0 = 2000)
amort.exporter_flux_csv('flux_conservation.csv')

visu = TraceurCSV('flux_conservation.csv')
visu.tracer(col_x = 'theta_pi', col_y = ['q_integre', 'q_final'],
            nom_pdf = 'conservation_flux.pdf', 
            titre = 'Vérification de la conservation de la masse (Flux q)',
            xlabel = r'Coordonnée angulaire $\theta / \pi$',
            ylabel = r'Débit volumique $q$ (m²/s)')

# Vitesse Max
amort.exporter_vitesse_max_csv('vitesse_max.csv')

visu = TraceurCSV('vitesse_max.csv')
visu.tracer(col_x = 'theta_pi', col_y = 'u_max_star', nom_pdf = 'vitesse_max.pdf',
            titre = 'Vitesse maximale adimensionnelle sur le tour',
            xlabel = r'Coordonnée angulaire $\theta / \pi$', 
            ylabel = r'Vitesse maximale $u_{max}^*$')