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


# Create fig that will show the intensity as a function of time
def PrintIntensityTime(Peaks, LIBS, Cutoffs, FinalCutoff):

    # Set a Seaborn style
    #sns.set(style="whitegrid")

    # Create the scatter plot
    fig, axs = plt.subplots(1, 1, figsize=(8, 4))
    
    # Add an asymptote for cleaner result
    #plt.axvline(x=(cutoff + 0.5), color='black', linestyle='--', linewidth=1.5)
    
    # Add the lines for each peak
    # PeakIs première colonne est #bin, seconde est #pic dans bin, troisième est la longueur d'onde
    for i in range(len(Peaks)):
        pic = Peaks[i]
        wavelength = pic[2]
        
        # Aller chercher la wavelength dans LIBS
        # 1 dataframe par position dans LIBS, première colonne est lambda, le reste est chaque tir
        intpeak = int(LIBS[0].index[LIBS[0].iloc[:, 0] == wavelength].values[0])
        PeakIntensities = LIBS[0].iloc[intpeak].tolist()[1:-3]
        
        # Add the line to the plot
        PlotIntensityTimeLine(fig, axs, PeakIntensities, Cutoffs[i], wavelength)


    plt.title(f'Vector Plot: First {FinalCutoff} Points removed from LIBS analysis')
    plt.grid()
    plt.legend()
    plt.show()
    
    plt.close('all')
    
# Adds one line to the figure
def PlotIntensityTimeLine(fig, axs, intensities, cutoff, wavelength):
    # Create a DataFrame for Seaborn
    df = pd.DataFrame({
        'Nb_Shots': [a+1 for a in np.arange(len(intensities))],
        'Intensity': intensities,
        'Color': ['red' if i < cutoff else 'blue' for i in range(len(intensities))]
    })
    
    # Add the lineplot
    #sns.lineplot(data=df, x='Nb_Shots', y='Intensity', hue='Color', marker='X', palette={'red', 'blue'}, label=str(wavelength))

    plt.scatter(df['Nb_Shots'], df['Intensity'], c=df['Color'], label=str(wavelength))