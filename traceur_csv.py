#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 02:12:38 2026

@author: cartier

READ ME

Ce fichier est là pour tracer les différentes courbes en essayant de généraliser 
la création de graphes. Pour que les données soient tracées il faut que les fichiers 
respectent l'implémentation demandée.
"""

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

class TraceurCSV:
    def __init__(self, nom_fichier, dossier_source = 'Results', separateur = ';'):
        chemin_csv = nom_fichier if os.path.exists(nom_fichier) else os.path.join(dossier_source, nom_fichier)
        self.df = pd.read_csv(chemin_csv, sep = separateur)

    def tracer(self, col_x, col_y, nom_pdf, titre = "", xlabel = None, ylabel = None, 
               ylabel_droit = None, secondary_y = None, dossier_fig = 'Figures'):
        os.makedirs(dossier_fig, exist_ok = True)
        if not nom_pdf.endswith('.pdf'):
            nom_pdf += '.pdf'
        chemin_pdf = os.path.join(dossier_fig, nom_pdf)

        ax = self.df.plot(x = col_x, y = col_y, secondary_y = secondary_y, mark_right = False, 
                          figsize = (8, 5), grid = True, title = titre)
        
        ax.set_xlabel(xlabel if xlabel else col_x)
        ax.set_ylabel(ylabel if ylabel else str(col_y))
        
        # Condition pour les petites et grandes valeurs sur les axes pour les mettre en puissance de 10

        ymin, ymax = ax.get_ylim()
        if max(abs(ymin), abs(ymax)) < 1e-1 or max(abs(ymin), abs(ymax)) > 1e3:
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

        if secondary_y and hasattr(ax, 'right_ax'):
            ax.right_ax.set_ylabel(ylabel_droit if ylabel_droit else str(secondary_y))
            
            ymin_droit, ymax_droit = ax.right_ax.get_ylim()
            if max(abs(ymin_droit), abs(ymax_droit)) < 1e-1 or max(abs(ymin_droit), abs(ymax_droit)) > 1e3:
                ax.right_ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

        plt.tight_layout()
        plt.savefig(chemin_pdf, format = 'pdf')
        plt.show()

    @staticmethod
    def comparer_csv(fichiers, labels, nom_pdf, col_x = 'theta_pi', col_y = 'P_star',
                     titre = "", xlabel = None, ylabel = None, dossier_source = 'Results',
                     dossier_fig = 'Figures', separateur = ';'):
        """
        :param fichiers: Nom du fichier csv 
        :param labels: Nom des labels pour savoir si un article est lié 
        :param nom_pdf: Nom de la figure enregistrée
        :param col_x: Nom de la colonne CSV à utiliser pour l'axe X
        :param col_y: Nom de la colonne CSV à utiliser pour l'axe Y
        :param titre: Titre du graphe
        :param xlabel, ylabel: kwargs
        """
        
        os.makedirs(dossier_fig, exist_ok = True)
        if not nom_pdf.endswith('.pdf'):
            nom_pdf += '.pdf'
        chemin_pdf = os.path.join(dossier_fig, nom_pdf)

        plt.figure(figsize = (8, 5))
        for fichier, label in zip(fichiers, labels):
            chemin_csv = fichier if os.path.exists(fichier) else os.path.join(dossier_source, fichier)
            df = pd.read_csv(chemin_csv, sep = None, engine = 'python')
            
            if 'Article' in label:
                plt.plot(df[col_x], df[col_y], marker = 'o', linestyle = '', markersize = 4, label = label)
            else:
                plt.plot(df[col_x], df[col_y], linestyle = '-', linewidth = 2, label = label)

        plt.xlabel(xlabel if xlabel else col_x)
        plt.ylabel(ylabel if ylabel else col_y)
        plt.title(titre)
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(chemin_pdf, format = 'pdf')
        plt.show()
        
    @staticmethod
    def tracer_frontieres_noyau(chemin_csv, titre, nom_pdf, dossier_source = 'Results', 
                                dossier_fig = 'Figures'):
        """
        :param chemin_csv: Emplacement du fichier dans un dictionnaire
        :param titre: Titre du graphe
        :param nom_pdf: Nom de la figure enregistrée
        """
        chemin_complet = chemin_csv if os.path.exists(chemin_csv) else os.path.join(dossier_source, chemin_csv)
        os.makedirs(dossier_fig, exist_ok = True)
        if not nom_pdf.endswith('.pdf'):
            nom_pdf += '.pdf'
        chemin_pdf = os.path.join(dossier_fig, nom_pdf)
    
        df = pd.read_csv(chemin_complet, sep = ';', engine = 'python')
        
        plt.figure(figsize = (8, 5))
        plt.plot(df['theta_pi'], df['ha_star'], label = r'$h_a^*$ (Limite inférieure)', 
                 color = 'blue', linewidth = 2)
        plt.plot(df['theta_pi'], df['hb_star'], label = r'$h_b^*$ (Limite supérieure)', 
                 color = 'purple', linewidth = 2)
        plt.fill_between(df['theta_pi'], df['ha_star'], df['hb_star'],
                         color = 'lightgrey', label = 'Noyau rigide')
        
        plt.ylim(0, 1)
        plt.xlabel(r'Coordonnée angulaire $\theta / \pi$')
        plt.ylabel(r'Coordonnée relative $y/h$')
        plt.title(titre)
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(chemin_pdf, format = 'pdf')
        plt.show()
     
    @staticmethod
    def plot_3D(fichiers_csv, y_values, x_col, z_col, titre, nom_pdf, 
                        xlabel, ylabel, zlabel, colorbar_label = "", 
                        cmap = 'plasma', dossier_source = 'Results', dossier_fig = 'Figures'):
        """
        :param fichiers_csv: Liste des noms de fichiers CSV
        :param y_values: Liste ou array des valeurs correspondant à chaque fichier (Axe Y)
        :param x_col: Nom de la colonne CSV à utiliser pour l'axe X
        :param z_col: Nom de la colonne CSV à utiliser pour l'axe Z (hauteur)
        :param xlabel, ylabel, zlabel
        :param colorbar_label: Titre de la barre de couleur
        :param cmap: Nom de la colormap
        """
        os.makedirs(dossier_fig, exist_ok = True)
        
        if not nom_pdf.endswith('.pdf'):
            nom_pdf += '.pdf'
        chemin_pdf = os.path.join(dossier_fig, nom_pdf)
        
        Z_matrix = []
        x_array = None
        
        for fichier in fichiers_csv:
            chemin_complet = fichier if os.path.exists(fichier) else os.path.join(dossier_source, fichier)
            df = pd.read_csv(chemin_complet, sep = None, engine = 'python')
            
            if x_array is None:
                x_array = df[x_col].values
                
            Z_matrix.append(df[z_col].values)
            
        Z_matrix = np.array(Z_matrix)
        X, Y = np.meshgrid(x_array, y_values)
        
        fig = plt.figure(figsize = (10, 7))
        ax = fig.add_subplot(111, projection = '3d')
        
        surf = ax.plot_surface(X, Y, Z_matrix, cmap = cmap, edgecolor = 'none')
        
        if colorbar_label:
            fig.colorbar(surf, label = colorbar_label)
        else:
            fig.colorbar(surf)
            
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)
        ax.set_title(titre)
        
        plt.tight_layout()
        plt.savefig(chemin_pdf, format = 'pdf')
        plt.show()
    
    @staticmethod
    def plot_2D(nom_fichier, col_x, col_y, col_z, nom_pdf,
                           titre = "", xlabel = None, ylabel = None, zlabel = None,
                           dossier_source = 'Results', dossier_fig = 'Figures', separateur = ';'):
        """
        :param nom_fichier: Nom du fichier csv 
        :param col_x: Nom de la colonne CSV à utiliser pour l'axe X
        :param col_y: Nom de la colonne CSV à utiliser pour l'axe Y
        :param col_z: Nom de la colonne CSV à utiliser pour l'axe Z
        :param nom_pdf: Nom de la figure enregistrée
        :param titre: Titre de la figure
        :param xlabel, ylabel, zlabel: kwargs
        """
        os.makedirs(dossier_fig, exist_ok = True)
        if not nom_pdf.endswith('.pdf'):
            nom_pdf += '.pdf'

        chemin = nom_fichier if os.path.exists(nom_fichier) else os.path.join(dossier_source, 
                                                                              nom_fichier)

        df = pd.read_csv(chemin, sep = separateur, usecols = [col_x, col_y, col_z])
        
        grille = df.pivot_table(index = col_y, columns = col_x, values = col_z)

        plt.figure(figsize = (8, 5))
        contour = plt.contourf(grille.columns, grille.index, grille.values, levels = 40, 
                               cmap = 'viridis')
        lignes = plt.contour(grille.columns, grille.index, grille.values, 
                             levels = 10, colors = 'black', linewidths = 0.5)
        plt.clabel(lignes, inline = True, fontsize = 8)
        plt.colorbar(contour, label = zlabel if zlabel else col_z)
        plt.xlabel(xlabel if xlabel else col_x)
        plt.ylabel(ylabel if ylabel else col_y)
        plt.title(titre)
        plt.tight_layout()
        plt.savefig(os.path.join(dossier_fig, nom_pdf), format = 'pdf')
        plt.show()