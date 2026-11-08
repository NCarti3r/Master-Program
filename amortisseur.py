#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 00:22:28 2026

@author: cartier


READ ME

Ce fichier permet la résolution dans le cas de l'amortisseur le fichier est appelé
dans le "main_amortisseur". 
"""

import numpy as np

from modules.class_fille.export_csv import ExportCSV

class Amortisseur(ExportCSV):
    """Classe fille pour l'amortisseur à film écrasé (U1 = U2 = -omega * R) """
    
    def __init__(self, R, c, epsilon, omega, mu, tau0):
        super().__init__(R, c, epsilon, omega, mu, tau0)
        self.U1 = -self.omega * self.R
        self.U2 = self.U1
        
    def get_q_init(self):
        denominateur = (2.0 + self.epsilon**2)
        numerateur = -2.0 * self.omega * self.R * self.c * (1.0 - self.epsilon**2)
        return numerateur / denominateur
    
    def dp_dx_cas_b(self, h, q):
        """Résolution de l'équation cubique / (h * dp/dx)^3 + 3 *
    (dp/dx)^2 * [4*mu*(omega*R*h + q) +- h^2 * tau0] - (+- 4 * tau0^3) = 0"""
        if self.tau0 == 0:
            return self.dp_dx_cas_a(h, q)

        flux_relatif = self.omega * self.R * h + q

        if abs(flux_relatif) < 1e-12:
            return 0.0

        if flux_relatif > 0 :
            signe = 1.0 
        else:
            signe = -1.0

        c3 = h**3
        c2 = 3.0 * (4.0 * self.mu * flux_relatif + signe * (h**2) * self.tau0)
        c1 = 0.0
        c0 = -signe * 4.0 * (self.tau0**3)

        racines = np.roots([c3, c2, c1, c0])

        racines_reelles = [r.real for r in racines if abs(r.imag) < 1e-7]
        racines_physiques = [r for r in racines_reelles if (r * signe < 0) and 
                             (abs(h * r) >= 2.0 * self.tau0 - 1e-7)]

        return racines_physiques[0]