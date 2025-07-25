import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import math


def PrintSingleSpectra(spectra, numshots, name, fig, axes, axenb):
    if fig == None:
        ax = sns.lineplot(data=spectra, x=0, y=(numshots+1))
        ax.fill_between(spectra.iloc[:,0], spectra.iloc[:,(numshots+2)], spectra.iloc[:,(numshots+3)], alpha=0.2)
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Intensity")
        plt.title(name)
        plt.show()
        plt.close('all')

    else:
        sns.lineplot(data=spectra, x=0, y=(numshots+1), ax=axes[axenb])
        axes[axenb].fill_between(spectra.iloc[:,0], spectra.iloc[:,(numshots+2)], spectra.iloc[:,(numshots+3)], alpha=0.2)
        axes[axenb].set_xlabel("Wavelength (nm)")
        axes[axenb].set_ylabel("Intensity")
        axes[axenb].set_title(name)

def PrintSpectrasSubplots(spectras, numshots, names):
    max_x = 4
    max_y = 4

    # Définir les dimensions du graphique avec les subplots
    nbfigs = len(spectras)
    if nbfigs <= max_x:
        lenx = max_x/2
        leny = math.ceil(nbfigs/lenx)
    else:
        if nbfigs <= max_x*math.ceil(max_y/2):
            lenx = max_x
            leny = math.ceil(max_y/2)

        else:
            if nbfigs <= max_x * max_y:
                lenx = max_x
                leny = max_y

            else:
                print("Too many graphs to plot for now, only the first " + str(max_x*max_y) + " will be plotted.")
                lenx = max_x
                leny = max_y

    # Créer le graphique
    fig, axes = plt.subplots(int(leny), int(lenx))
    for i in range(0, len(spectras), 1):
        PrintSingleSpectra(spectras[i], numshots[i], names[i], fig, axes, i)

    plt.tight_layout()
    plt.show()
    
    plt.close('all')

# Pour avoir un fichier .txt des fichiers libs obtenus
def PrintAvgLIBSTxt(AvgLIBS, title):
    # Ne va imprimer que les deux premières colonnes, sans info, donc juste wavelength et intensité
    AvgLIBS.iloc[:, :2].to_csv(title, index=False, header=False, sep='\t')

