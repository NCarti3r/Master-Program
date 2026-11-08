    #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 00:29:54 2026

@author: cartier

READ ME

Ce fichier permet la résolution dans le cas du palier lisse le fichier est appelé
dans le "main_palier_lisse". 
"""
import numpy as np

from modules.class_fille.export_csv import ExportCSV

class PalierLisse(ExportCSV):
  """Classe fille pour le palier lisse (U1 = omega * R / U2 = 0)"""

  def __init__(self, R, c, epsilon, omega, mu, tau0):
    super().__init__(R, c, epsilon, omega, mu, tau0)
    self.U1 = self.omega * self.R
    self.U2 = 0.0

  def get_q_init(self):
      numerateur = self.omega * self.R * self.c * (1.0 - self.epsilon**2)
      denominateur = 2.0 + self.epsilon**2
      return numerateur / denominateur

  def dp_dx_cas_b(self, h, q):
    """Résolution du polynôme de degré 5"""
    if self.tau0 == 0:
      return self.dp_dx_cas_a(h, q)

    dp_dx_newton = self.dp_dx_cas_a(h, q)
    if abs(dp_dx_newton) < 1e-12:
      return 0.0

    if dp_dx_newton < 0:
        signe = 1.0
    else:
        signe = -1.0

    U_somme = self.U1 + self.U2
    U_diff = self.U1 - self.U2
    tau0 = self.tau0
    mu = self.mu

    c5 = -(h**5)
    
    c4 = (-signe * 7.0 * (h**4) * tau0 + 6.0 * mu * (h**3) * U_somme 
          - 12.0 * mu * q * (h**2))
    
    c3 = (-16.0 * (h**3) * (tau0**2) + signe * 24.0 * mu * (h**2) * tau0 * U_somme 
          - signe * 48.0 * mu * q * h * tau0)
    
    c2 = (-signe * 8.0 * (h**2) * (tau0**3) + 24.0 * mu * h * (tau0**2) * U_somme 
          + signe * 12.0 * (mu**2) * tau0 * (U_diff**2) - 48.0 * mu * q * (tau0**2))
    
    c1 = 16.0 * h * (tau0**4)
    
    c0 = signe * 16.0 * (tau0**5)

    racines = np.roots([c5, c4, c3, c2, c1, c0])
    racines_reelles = [r.real for r in racines if abs(r.imag) < 1e-6]

    for r in racines_reelles:
      if r * signe < 0 and abs(h * r + signe * 2.0 * tau0) > 1e-9:
        ha = self.calculer_ha(h, r)
        hb = self.calculer_hb(h, r)
        if 0.0 < ha < hb < h:
          return r

    return dp_dx_newton