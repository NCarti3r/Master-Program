# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 20:55:36 2026

@author: sproc

READ ME

Dans cette Class mère nous déterminons le cas dans lequel nous nous trouvons 
puis résolvons les différentes équations pour calculer le gradient de pression,
les hauteurs ha et hb, la vitesse, les contraintes pariétales et le flux intégrale. 
"""

import numpy as np

from scipy.integrate import cumulative_trapezoid as cumtrapz
from scipy.optimize import root_scalar


class ModelBingham:
    def __init__(self, R, c, epsilon, omega, mu, tau0):
        self.R = R
        self.c = c
        self.epsilon = epsilon

        self.omega = omega
        self.mu = mu
        self.tau0 = tau0

        self.U1 = 0.0
        self.U2 = 0.0

    def calcule_epaisseur(self, theta):
        return self.c * (1 + self.epsilon * np.cos(theta))

    def get_q_init(self):
        raise NotImplementedError("Cette méthode doit être codée dans la classe fille.")

    """Calculs locaux des gradients de pression"""
    def dp_dx_cas_a(self, h, q):
        return ((12 * self.mu) / h**3) * (((self.U1 + self.U2) * h / 2) - q)

    def dp_dx_cas_c(self, h, q):
        den = 9 * (q - h * self.U1)**2
        if den < 1e-300:
            return np.nan
        return (2 * self.mu * (self.U2 - self.U1)**3) / den

    def dp_dx_cas_d(self, h, q):
        den = 9 * (q - h * self.U2)**2
        if den < 1e-300:
            return np.nan
        return (2 * self.mu * (self.U1 - self.U2)**3) / den

    def dp_dx_cas_b(self, h, q):
        raise NotImplementedError("Cette méthode doit être codée dans la classe fille (Amortisseur ou PalierLisse).")

    """Frontières du noyau"""
    def calculer_ha(self, h, dp_dx):
        """Éq. (16c) : valable pour le noyau FLOTTANT (cas b) et pour le test de cas."""
        signe = 1 if dp_dx < 0 else -1
        if dp_dx == 0:
            return 0
        numerateur = ((h * dp_dx)**2 + signe * 4 * h * self.tau0 * dp_dx 
                      + 4 * self.tau0**2
                      + 2 * self.mu * dp_dx * (self.U1 - self.U2))
        denominateur = 2 * dp_dx * (h * dp_dx + signe * 2 * self.tau0)
        return numerateur / denominateur

    def calculer_hb(self, h, dp_dx):
        """Éq. (16d) : valable pour le noyau FLOTTANT (cas b) et pour le test de cas."""
        signe = 1 if dp_dx < 0 else -1
        if dp_dx == 0:
            return h
        numerateur = (h * dp_dx)**2 - 4 * self.tau0**2 + 2 * self.mu * dp_dx * (self.U1 - self.U2)
        denominateur = 2 * dp_dx * (h * dp_dx + signe * 2 * self.tau0)
        return numerateur / denominateur

    def epaisseur_zone_fluide(self, dp_dx):
        """Épaisseur de la zone qui s'écoule dans un cas à noyau attaché (c ou d).
        Profil parabolique de pente nulle au bord du noyau : |dU| = |dp| d^2 / (2 mu)
        (équivalent à Éq. 19b et 22b)."""
        return np.sqrt(2 * self.mu * abs(self.U2 - self.U1) / abs(dp_dx))

    def frontieres(self, h, dp_dx, cas):
        """Retourne (ha, hb) cohérents avec le cas retenu."""
        if cas == 'a':
            return np.nan, np.nan
        if cas == 'c':
            return 0.0, h - self.epaisseur_zone_fluide(dp_dx)
        if cas == 'd':
            return self.epaisseur_zone_fluide(dp_dx), h
        return self.calculer_ha(h, dp_dx), self.calculer_hb(h, dp_dx)

    def classer_local(self, h, q):
        """Détermine le cas (a, b, c ou d) à partir du débit q, et renvoie (cas, dp_dx).

        (a) pas de noyau   : la contrainte newtonienne reste partout >= tau0
                             |mu (U2-U1)/h| - |dp_a| h/2 >= tau0
        (c) noyau attaché en bas, (d) attaché en haut : on calcule le gradient du cas
            (Éq. 20 / 23) puis on vérifie
              - le signe de q - h*U (le profil est monotone dans la zone fluide),
              - 0 < épaisseur de la zone fluide < h,
              - tau0 >= |dp| (h - d) / 2  (généralisation de l'Éq. 30 : le noyau a
                atteint la paroi)
        (b) noyau flottant : sinon.
        """
        dp_a = self.dp_dx_cas_a(h, q)
        if self.tau0 == 0:
            return 'a', dp_a

        dU = self.U2 - self.U1
        if abs(self.mu * dU / h) - abs(dp_a) * h / 2.0 >= self.tau0:
            return 'a', dp_a

        if abs(dU) > 1e-12:
            # Cas (c) : noyau attaché à la paroi inférieure (vitesse U1)
            dp_c = self.dp_dx_cas_c(h, q)
            if np.isfinite(dp_c) and dp_c != 0 and (q - h * self.U1) * dU > 0:
                d = self.epaisseur_zone_fluide(dp_c)
                if 0 < d < h and self.tau0 >= 0.5 * abs(dp_c) * (h - d):
                    return 'c', dp_c

            # Cas (d) : noyau attaché à la paroi supérieure (vitesse U2)
            dp_d = self.dp_dx_cas_d(h, q)
            if np.isfinite(dp_d) and dp_d != 0 and (q - h * self.U2) * (-dU) > 0:
                d = self.epaisseur_zone_fluide(dp_d)
                if 0 < d < h and self.tau0 >= 0.5 * abs(dp_d) * (h - d):
                    return 'd', dp_d

        dp_b = self.dp_dx_cas_b(h, q)
        ha, hb = self.calculer_ha(h, dp_b), self.calculer_hb(h, dp_b)
        if not (0.0 < ha < hb < h):
            return 'a', dp_a
        return 'b', dp_b

    def cas_depuis_dp(self, h, dp_dx):
        """Classement quand dp/dx est imposé (profils locaux, Fig. 3)."""
        if self.tau0 == 0:
            return 'a'
        dU = self.U2 - self.U1
        if abs(self.mu * dU / h) - abs(dp_dx) * h / 2.0 >= self.tau0:
            return 'a'
        ha = self.calculer_ha(h, dp_dx)
        hb = self.calculer_hb(h, dp_dx)
        if 0 < ha < hb < h:
            return 'b'
        if abs(dU) < 1e-12:
            return 'b'
        # noyau attaché : c si dp a le signe de (U2-U1), d sinon
        return 'c' if dp_dx * dU > 0 else 'd'

    def calculer_dp_dx_local(self, h, q, type_palier = 'palier'):
        return self.classer_local(h, q)[1]

    """Calculs locaux de la vitesse"""
    def vitesse_cas_a(self, h, dp_dx, y):
        terme_glissement = self.U1 + ((self.U2 - self.U1) / h) * y
        terme_pression = (1 / (2 * self.mu)) * dp_dx * (y**2 - y * h)
        return terme_glissement + terme_pression

    def vitesse_cas_c(self, h, dp_dx, y, hb):
        # noyau rigide de 0 à hb, puis profil parabolique de pente nulle en hb
        if y <= hb:
            return self.U1
        return self.U2 + (dp_dx / (2 * self.mu)) * ((y - h)**2 - 2 * (hb - h) * (y - h))

    def vitesse_cas_d(self, h, dp_dx, y, ha):
        # noyau rigide de ha à h, profil parabolique de pente nulle en ha en dessous
        if y >= ha:
            return self.U2
        return self.U2 + (dp_dx / (2 * self.mu)) * (y - ha)**2

    def vitesse_cas_b(self, h, dp_dx, y, ha, hb):
        signe = 1 if dp_dx < 0 else -1
        tau0_signe = signe * self.tau0

        num_tau1 = h**2 * dp_dx + signe * 2 * h * self.tau0 + 2 * self.mu * (self.U1 - self.U2)
        den_tau1 = 2 * h * dp_dx + signe * 4 * self.tau0
        tau_1 = - (dp_dx / den_tau1) * num_tau1
        tau_2 = tau_1 + h * dp_dx

        if y <= ha:
            return self.U1 - ((tau0_signe - tau_1) / self.mu) * y + (1 / (2 * self.mu)) * dp_dx * y**2
        elif ha < y < hb:
            return self.U1 - ((tau0_signe - tau_1) / self.mu) * ha + (1 / (2 * self.mu)) * dp_dx * ha**2
        else:
            return self.U2 - ((-tau0_signe - tau_2) / self.mu) * (y - h) + (1 / (2 * self.mu)) * dp_dx * (y - h)**2

    def calculer_vitesse_locale(self, h, dp_dx, y, cas = None):
        """Vitesse en y. Si 'cas' n'est pas fourni, il est déduit de dp_dx."""
        if cas is None:
            cas = self.cas_depuis_dp(h, dp_dx)
        if cas == 'a':
            return self.vitesse_cas_a(h, dp_dx, y)
        ha, hb = self.frontieres(h, dp_dx, cas)
        if cas == 'c':
            return self.vitesse_cas_c(h, dp_dx, y, hb)
        if cas == 'd':
            return self.vitesse_cas_d(h, dp_dx, y, ha)
        return self.vitesse_cas_b(h, dp_dx, y, ha, hb)

    """Résolution globale"""
    def calcule_erreur_periodicite(self, q):
        self.theta_array = np.linspace(0, 2 * np.pi, 200)
        dp_dx_array = np.zeros_like(self.theta_array)
        
        for i, theta in enumerate(self.theta_array):
            h_local = self.calcule_epaisseur(theta)
            dp_dx_array[i] = self.calculer_dp_dx_local(h_local, q)
            
        # dx = R * dtheta, donc dp = R * dp_dx * dtheta
        self.P_array = cumtrapz(dp_dx_array * self.R, self.theta_array, initial=0)
        
        erreur = self.P_array[-1] - self.P_array[0]
        return erreur
    
    def resoudre_palier(self):
        q_init = self.get_q_init()
        try:
            resultat = root_scalar(self.calcule_erreur_periodicite, x0 = q_init,
                                   x1 = q_init+1e-6, method = 'secant')
            self.q_final = resultat.root
            self.calcule_erreur_periodicite(self.q_final)  # P_array cohérent avec q_final
            
        except ValueError:
            print(f"Attention, le solveur n'a pas convergé pour tau0 = {self.tau0}." 
                  "Utilisation de q_init par défaut.")
            self.q_final = q_init
            self.calcule_erreur_periodicite(self.q_final)
            
        # p* = p * c^2 / (mu * omega * R^2)
        facteur_adimensionnement = (self.c**2) / (self.mu * self.omega * self.R**2)
        self.P_star_array = self.P_array * facteur_adimensionnement
        return self.q_final
    
    def calculer_contrainte_parietale(self, h, dp_dx, cas = None):
        """
        Calcule les contraintes de cisaillement pariétales tau_1 (en y = 0) et tau_2 (en y = h)
        selon le régime local ('a', 'b', 'c' ou 'd').
        """
        if cas is None:
            if hasattr(self, 'cas_depuis_dp'):
                cas = self.cas_depuis_dp(h, dp_dx)
            else:
                cas, _, _ = self.identifier_regime(h, dp_dx)

        if cas == 'a':
            du_dy_0 = ((self.U2 - self.U1) / h) - (h / (2.0 * self.mu)) * dp_dx
            signe_cisaillement = 1.0 if du_dy_0 >= 0 else -1.0
            tau_1 = self.mu * du_dy_0 + signe_cisaillement * self.tau0
            tau_2 = tau_1 + h * dp_dx

        elif cas == 'c':
            tau_2 = -self.tau0 + np.sqrt(abs(2.0 * self.mu * dp_dx * (self.U2 - self.U1)))
            tau_1 = tau_2 - h * dp_dx

        elif cas == 'd':
            tau_1 = -self.tau0 - np.sqrt(abs(2.0 * self.mu * dp_dx * (self.U1 - self.U2)))
            tau_2 = tau_1 + h * dp_dx

        else:
            if abs(dp_dx) < 1e-12:
                tau_1 = 0.0
            else:
                signe = 1.0 if dp_dx < 0 else -1.0
                den = 2.0 * h * dp_dx + signe * 4.0 * self.tau0
                if abs(den) < 1e-12:
                    tau_1 = signe * self.tau0
                else:
                    num = (h**2) * dp_dx + signe * 2.0 * h * self.tau0 + 2.0 * self.mu * (self.U1 - self.U2)
                    tau_1 = -(dp_dx / den) * num
            tau_2 = tau_1 + h * dp_dx

        return tau_1, tau_2
    
    def calculer_flux_integral(self, theta, nb_points = 200):
        """
        Calcule le débit volumique local q en intégrant la vitesse
        """
        if not hasattr(self, 'q_final'):
            self.resoudre_palier()
            
        h = self.calcule_epaisseur(theta)
        cas, dp_dx = self.classer_local(h, self.q_final)
        
        y_array = np.linspace(0, h, nb_points)
        u_array = np.zeros_like(y_array)
        
        for i, y in enumerate(y_array):
            u_array[i] = self.calculer_vitesse_locale(h, dp_dx, y, cas)

        q_integre = np.trapezoid(u_array, y_array)
        
        return q_integre
    