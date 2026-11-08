#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 01:43:52 2026

@author: cartier

READ ME

Ce fichier permet d'extraire les données résolues dans les fichiers amortisseur 
et palier lisse pour ensuite les lire par le fichier traceur_csv.
"""
import os
import csv

import numpy as np

from modules.model_bingham import ModelBingham


class ExportCSV(ModelBingham):
    def exporter_pression_csv(self, nom_fichier, colonnes = None, dossier = 'Results'):
        if not hasattr(self, 'q_final'):
            self.resoudre_palier()

        if colonnes is None:
            colonnes = ['theta_pi', 
                        'h', 
                        'dp_dx', 
                        'P', 
                        'P_star', 
                        'ha_star', 
                        'hb_star', 
                        'tau_1', 
                        'tau_2', 
                        'tau_1_star', 
                        'tau_2_star']

        os.makedirs(dossier, exist_ok = True)
        if not nom_fichier.endswith('.csv'):
            nom_fichier += '.csv'
        chemin_csv = os.path.join(dossier, nom_fichier)
        
        facteur_tau_star = self.c / (self.mu * self.omega * self.R)

        with open(chemin_csv, 'w', newline = '', encoding = 'utf-8') as f:
            writer = csv.writer(f, delimiter = ';')
            writer.writerow(colonnes)

            for i, theta in enumerate(self.theta_array):
                h = self.calcule_epaisseur(theta)
                cas, dp_dx = self.classer_local(h, self.q_final)
                ha, hb = self.frontieres(h, dp_dx, cas)
                tau_1, tau_2 = self.calculer_contrainte_parietale(h, dp_dx, cas)

                donnees = {'theta': theta, 'theta_pi': theta / np.pi, 'h': h, 
                           'dp_dx': dp_dx, 'P': self.P_array[i], 
                           'P_star': self.P_star_array[i], 'ha': ha, 'hb': hb, 
                           'ha_star': np.clip(ha / h, 0, 1), 
                           'hb_star': np.clip(hb / h, 0, 1), 
                           'tau_1': tau_1,
                           'tau_2': tau_2,
                           'tau_1_star': tau_1 * facteur_tau_star,
                           'tau_2_star': tau_2 * facteur_tau_star,
                           'q': self.q_final, 'cas': cas}
                writer.writerow([donnees[col] for col in colonnes])

        return chemin_csv
    
    def exporter_vitesse_max_csv(self, nom_fichier, dossier = 'Results'):
        """Exporte la vitesse maximale adimensionnelle sur tout le tour du palier."""
        if not hasattr(self, 'q_final'):
            self.resoudre_palier()

        os.makedirs(dossier, exist_ok = True)
        if not nom_fichier.endswith('.csv'):
            nom_fichier += '.csv'
        chemin_csv = os.path.join(dossier, nom_fichier)

        with open(chemin_csv, 'w', newline = '', encoding = 'utf-8') as f:
            writer = csv.writer(f, delimiter = ';')
            writer.writerow(['theta_pi', 'u_max_star'])

            for theta in self.theta_array:
                h = self.calcule_epaisseur(theta)
                cas, dp_dx = self.classer_local(h, self.q_final)

                y_array = np.linspace(0, h, 100)
                u_max = max(self.calculer_vitesse_locale(h, dp_dx, y, cas) for y in y_array)
                
                u_max_star = u_max / (self.omega * self.R)
                writer.writerow([theta / np.pi, u_max_star])

        return chemin_csv

    def exporter_vitesse_locale_csv(self, dp_star, tau0_star, nom_fichier, dossier = 'Results'):
        """Profil de vitesse local avec gradient de pression imposé.
        Adimensionnement local avec h = c et U = omega*R :
            dp* = h^2/(mu U) dp/dx      tau0* = tau0 h/(mu U)
        Aucun débit q ni bouclage de périodicité : dp/dx est une donnée.
        """
        os.makedirs(dossier, exist_ok = True)
        if not nom_fichier.endswith('.csv'):
            nom_fichier += '.csv'
        chemin_csv = os.path.join(dossier, nom_fichier)

        h = self.c
        U = self.omega * self.R
        dp_dx = dp_star * self.mu * U / h**2

        tau0_sauve = self.tau0
        self.tau0 = tau0_star * self.mu * U / h
        try:
            cas = self.cas_depuis_dp(h, dp_dx)
            y_array = np.linspace(0, h, 200)
            with open(chemin_csv, 'w', newline = '', encoding = 'utf-8') as f:
                writer = csv.writer(f, delimiter = ';')
                writer.writerow(['y_star', 'u_star'])
                for y in y_array:
                    u = self.calculer_vitesse_locale(h, dp_dx, y, cas)
                    writer.writerow([y / h, u / U])
        finally:
            self.tau0 = tau0_sauve
        return chemin_csv
    
    def exporter_flux_csv(self, nom_fichier, dossier = 'Results'):
        if not hasattr(self, 'q_final'):
            self.resoudre_palier()
            
        os.makedirs(dossier, exist_ok = True)
        if not nom_fichier.endswith('.csv'):
            nom_fichier += '.csv'
        chemin_csv = os.path.join(dossier, nom_fichier)

        with open(chemin_csv, 'w', newline = '', encoding = 'utf-8') as f:
            writer = csv.writer(f, delimiter = ';')
            writer.writerow(['theta_pi', 'q_integre', 'q_final'])
            
            for theta in self.theta_array:
                q_int = self.calculer_flux_integral(theta)
                writer.writerow([theta / np.pi, q_int, self.q_final])
                
        return chemin_csv
    
    def exporter_champ_vitesse_2D_csv(self, nom_fichier, nb_theta = 60, nb_y = 40, dossier = 'Results'):
        """Exporte la carte 2D complète de la vitesse (theta, y) sur tout le tour du palier."""
        if not hasattr(self, 'q_final'):
            self.resoudre_palier()

        os.makedirs(dossier, exist_ok = True)
        if not nom_fichier.endswith('.csv'):
            nom_fichier += '.csv'
        chemin_csv = os.path.join(dossier, nom_fichier)

        # On fixe les grilles adimensionnelles AVANT la boucle pour éviter tout décalage d'arrondi
        theta_array = np.linspace(0, 2 * np.pi, nb_theta)
        y_star_array = np.linspace(0.0, 1.0, nb_y)

        with open(chemin_csv, 'w', newline = '', encoding = 'utf-8') as f:
            writer = csv.writer(f, delimiter = ';')
            writer.writerow(['theta_pi', 'y_star', 'u_star'])

            for theta in theta_array:
                h = self.calcule_epaisseur(theta)
                cas, dp_dx = self.classer_local(h, self.q_final)
                theta_pi = round(theta / np.pi, 6)
                
                for y_star in y_star_array:
                    y = y_star * h  # Conversion en coordonnée physique pour le calcul
                    u = self.calculer_vitesse_locale(h, dp_dx, y, cas)
                    writer.writerow([theta_pi, round(y_star, 6), u / (self.omega * self.R)])

        return chemin_csv
    