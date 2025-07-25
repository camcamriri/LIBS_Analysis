import numpy as np
import pandas as pd

from scipy.signal import find_peaks
from scipy.stats import mannwhitneyu
import statistics
from prettytable import PrettyTable

from Inputs import *
from Outputs import *
from SideFunctions import *


# Calculer le spectre moyen de tous les spectres sauvés
def GetAvgLIBS(listdflibs, nshots, showgraph, print_txt):
    AvgLIBS = []
    wavelengths = []
    intensities = []

    for i in range(len(listdflibs)):
        libs = listdflibs[i]
        for j in range(len(libs[0])):
            # Get la ligne à son wavelength
            wavel = libs.iloc[j][0]
            measurements = libs.iloc[j].tolist()[1:nshots[i]+1]

            # Vérifier si wavelength déjà pris et rajouter, puis add à la bonne liste dans intensities
            if wavel not in wavelengths:
                wavelengths.append(wavel)
                idx_wave = wavelengths.index(wavel)
                intensities.append([])
            else:
                idx_wave = wavelengths.index(wavel)

            # Ajouter à intensities
            for val in measurements:
                intensities[idx_wave].append(val)

    # Faire la moyenne de tous les intensities
    for i in range(len(intensities)):
        avg = np.mean(intensities[i])
        stdev = 2*np.std(intensities[i])

        newline = [wavelengths[i], avg, avg-stdev, avg+stdev]
        AvgLIBS.append(newline)

    # Transformer en dataframe et trier
    AvgLIBS = pd.DataFrame(AvgLIBS)
    AvgLIBS = AvgLIBS.sort_values(0)
    
    if showgraph:
        avgspctitle = "Averaged Spectra of " + str(nshots[0]) + " shots/position in " + str(len(nshots)) + " positions"
        PrintSingleSpectra(AvgLIBS, 0, avgspctitle, None, None, 0)
        
    if print_txt:
        title = "AveragedSpectra.txt"
        PrintAvgLIBSTxt(AvgLIBS, title)

    return AvgLIBS

# Connaître les limites de détection en longueur d'ondes des spectromètres utilisés
def IdentifyAndSplitOverlaps(intervals):
    # NB : Ne fonctionne que si les overlaps n'arrivent qu,entre deux sections, si un spectre overlap deux autres, méthode va échouer
    # Sort pour avoir en ordre croissant de premiere valeur
    intervals.sort(key=lambda x: x[0])

    # Vérifier s'il y a des overlaps
    overlaps = []
    for i in range(1, len(intervals)):
        prev_start, prev_end = intervals[i - 1]
        curr_start, curr_end = intervals[i]
        
        # Check if there is an overlap
        if curr_start < prev_end:
            overlaps.append((prev_start, prev_end, curr_start, curr_end))

    # Séparer les overlaps en trois nouveaux intervalles, la partie initiale sans overlap, l'overlap, partie finale sans overlap
    newintervals = []
    for i in range(len(overlaps)):
        firstnewint = [overlaps[i][0], overlaps[i][2]]
        secondnewint = [overlaps[i][2], overlaps[i][1]]
        thirdnewint = [overlaps[i][1], overlaps[i][3]]

        newintervals.append(firstnewint)
        newintervals.append(secondnewint)
        newintervals.append(thirdnewint)

    # Nouvelle liste avec les intervalles à garder
    FirstLambdaNewInts = list(map(list, zip(*newintervals)))[0]
    for i in range(len(intervals)):
        if intervals[i][0] not in FirstLambdaNewInts:
            newintervals.append(intervals[i])

    # Sort la nouvelle liste d'intervalles
    newintervals.sort(key=lambda x: x[0])

    
    return newintervals

# Fonction utilisée pour trouver l'intensité répétitivement
def GetIntensityFromPic(binLIBS, idxpic):
    Ibg = np.mean(binLIBS[idxpic[2]:idxpic[3]+1])
    ns = idxpic[1]-idxpic[0]+1
    sIs = np.sum(binLIBS[idxpic[0]:idxpic[1]+1])
    Is = sIs - ns*Ibg
    Is = round(Is)

    return Is

# Identifier l'intensité pour les longueurs d'ondes fournies
def PeakLambdaAndIntensityCompute(AvgLIBS, bins, thresholds, graphstoshow):
    idxpics = []
    peaksintensities = []

    
    for i in range(len(bins)):
        # Identification des pics pour le bin
        bin = bins[i]
        AvgLIBSbin = AvgLIBS.loc[(AvgLIBS[0] >= bin[0]) & (AvgLIBS[0] < bin[1])]
        binLIBS = np.array(list(map(list, zip(*np.array(AvgLIBSbin).tolist()))))

        peaks, _ = find_peaks(binLIBS[1], threshold=thresholds[i])

        if 1 in graphstoshow:
            plt.plot(binLIBS[0], binLIBS[1])
            plt.plot(binLIBS[0][peaks], binLIBS[1][peaks], "x")
            title = "Détection des pics directe sur bin avec lambda de " + str(bin[0]) + " à " + str(bin[1])
            plt.title(title)
            plt.show()

        # Pour chaque bin, identification des 4 valeurs de longueur d'onde à considérer pour calcul I
        idxpicsbin = []
        allidpicsbin = []
        allidavgbin = []
        for j in range(len(peaks)):
            idxpic = []
            peak = peaks[j]
            # Deux premiers sont sur la valeur du pic, à moitié de la hauteur de ce dernier
            Binmedian = statistics.median(binLIBS[1])
            Binmean = statistics.mean(binLIBS[1])
            roughbaseline = max(Binmedian, Binmean)
            Iini = (binLIBS[1][peak]+roughbaseline)/2
        
            # Trouver le moment où pic est à moitié avant pic
            k = 1
            while True:
                if binLIBS[1][peak-k] < Iini:
                    if k == 1:
                        idxpic.append(peak-k)
                        break
                    else:
                        idxpic.append(peak-k+1)
                        break

                else:
                    k += 1

            # Trouver le moment où pic est à moitié après pic
            k = 1
            while True:
                if binLIBS[1][peak+k] < Iini:
                    if k == 1:
                        idxpic.append(peak+k)
                        break
                    else:
                        idxpic.append(peak+k-1)
                        break

                else:
                    k += 1

            # Trouver les idx pour baseline
            boundmin = 0.15
            boundmax = 0.3
            # Test si c'est le premier pic
            if peak == peaks[0]:
                premierintervalle = [int(peak*0.7), int(peak*0.85)]
                # Vérifier pas juste un pic
                if len(peaks) == 1:
                    secondintervalle = [int((len(binLIBS[1])-peak)*boundmin+peak), int((len(binLIBS[1])-peak)*boundmax+peak)]
                else:
                    secondintervalle = [int((peaks[j+1]-peak)*boundmin+peak), int((peaks[j+1]-peak)*boundmax+peak)]

            # Test si c'est le dernier pic
            elif peak == peaks[-1]:
                premierintervalle = [int((peak-peaks[j-1])*(1-boundmax)+peaks[j-1]), int((peak-peaks[j-1])*(1-boundmin)+peaks[j-1])]
                secondintervalle = [int((len(binLIBS[1])-peak)*boundmin+peak), int((len(binLIBS[1])-peak)*(1-boundmax)+peak)]

            # Sinon
            else:
                premierintervalle = [int((peak-peaks[j-1])*(1-boundmax)+peaks[j-1]), int((peak-peaks[j-1])*(1-boundmin)+peaks[j-1])]
                secondintervalle = [int((peaks[j+1]-peak)*boundmin+peak), int((peaks[j+1]-peak)*boundmax+peak)]

            # Comparer les deux intervalles
            stdpremier = np.std(binLIBS[1][premierintervalle[0]:premierintervalle[1]+1])
            stdsecond = np.std(binLIBS[1][secondintervalle[0]:secondintervalle[1]+1])

            nbmindatapoints = 5
            if stdpremier < stdsecond:
                if premierintervalle[1]-premierintervalle[0] < nbmindatapoints:
                    idxpic.append(secondintervalle[0])
                    idxpic.append(secondintervalle[1])
                else:
                    idxpic.append(premierintervalle[0])
                    idxpic.append(premierintervalle[1])
            else:
                if secondintervalle[1]-secondintervalle[0] < nbmindatapoints:
                    idxpic.append(premierintervalle[0])
                    idxpic.append(premierintervalle[1])
                else:
                    idxpic.append(secondintervalle[0])
                    idxpic.append(secondintervalle[1])


            if 2 in graphstoshow:
                plt.plot(binLIBS[0], binLIBS[1])
                plt.plot(binLIBS[0][idxpic], binLIBS[1][idxpic], "x")
                title = "Détection des limites du pic"
                plt.title(title)
                plt.show()

            # Sauvegarder les infos obtenues
            idxpicsbin.append(idxpic)

            allidpicsbin.append(idxpic[0])
            allidpicsbin.append(idxpic[1])
            allidavgbin.append(idxpic[2])
            allidavgbin.append(idxpic[3])

            # Calculer l'intensité du pic avec les lambda maintenant définis
            Is = GetIntensityFromPic(binLIBS[1], idxpic)
            Ismin = GetIntensityFromPic(binLIBS[2], idxpic)
            Ismax = GetIntensityFromPic(binLIBS[3], idxpic)

            dI1 = Ismax - Is
            dI2 = Is-Ismin
            dI = max(dI1, dI2)

            # Sauvegarder l'I du pic
            peaksintensities.append([i+1, j+1, binLIBS[0][peak], Is, dI])

            #print('stop me here')

        # Sauvegarder pour ce bin spécifique les lambda
        idxpics.append(idxpicsbin)

        # Print graph du bin entier avec pics et averages gardés
        if 3  in graphstoshow:
            #spectra = pd.DataFrame(binLIBS)
            ax = sns.lineplot(data=AvgLIBSbin, x=0, y=1)
            ax.fill_between(AvgLIBSbin.iloc[:,0], AvgLIBSbin.iloc[:,2], AvgLIBSbin.iloc[:,3], alpha=0.2)
            plt.xlabel("Wavelength (nm)")
            plt.ylabel("Intensity")

            #plt.plot(binLIBS[0], binLIBS[1])
            plt.plot(binLIBS[0][peaks], binLIBS[1][peaks], "x", label="Peaks")
            plt.plot(binLIBS[0][allidpicsbin], binLIBS[1][allidpicsbin], "x", label="Peak Boundaries")
            plt.plot(binLIBS[0][allidavgbin], binLIBS[1][allidavgbin], "x", label="Baseline Boundaries")
            title = "Détection des pics directe sur bin avec lambda de " + str(bin[0]) + " à " + str(bin[1]) + "\navec les limites des pics et des averages"
            plt.title(title)
            plt.legend()
            plt.grid()
            plt.show()
            plt.close('all')


    # Imprimer les intensitées dans un joli tableau
    t = PrettyTable(["# Bin", "# Pic dans bin", "Longueur d'onde", "Intensité", "95% dI"])
    for pic in peaksintensities:
        t.add_row(pic)

    print(t)

    return idxpics, peaksintensities

# Identifier les pics qui sont différents pour les n premiers tirs à partir d'un test de moyenne égale
def IdentifyFirstShotsDifference(LIBS, PeakIs, ConfInterval, showgraph):
    # Initialiser, va être augmenté à chaque fois que c'est jugé nécessaire
    NTirsRemove = []
    
    # Faire la boucle dans PeakIs
    # PeakIs première colonne est #bin, seconde est #pic dans bin, troisième est la longueur d'onde
    for pic in PeakIs:
        wavelength = pic[2]
        
        # Aller chercher la wavelength dans LIBS
        # 1 dataframe par position dans LIBS, première colonne est lambda, le reste est chaque tir
        intpeak = int(LIBS[0].index[LIBS[0].iloc[:, 0] == wavelength].values[0])
        PeakIntensities = LIBS[0].iloc[intpeak].tolist()[1:-3]

        # Faire test d'ergodicité selon le pourcentage de certitude attendu par ConfInterval
        # Faire le test pour chaque pic dès le début et enlever dès que ça devient différent
        remove = GetNRemove(PeakIntensities, ConfInterval)    
        NTirsRemove.append(remove)
            
    # À partir de chaque NTirs à enlever, faire le choix du nombre de tirs à retirer globalement (voir si outlier, si oui, ignorer et passer au max suivant)
    finalCutoff = NTirsRemove[0]

    # Pour valider si bon sur le dataset       
    if showgraph:
        # Afficher dans un graphique y intensité, x # du tir et avoir une courbe par position + une asymptote verticale du nombre de tirs choisi à retirer
        PrintIntensityTime(PeakIs, LIBS, NTirsRemove, finalCutoff)

    return NTirsRemove