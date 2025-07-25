import numpy as np
import pandas as pd

from scipy.signal import find_peaks
from scipy.stats import mannwhitneyu
import statistics
from prettytable import PrettyTable




# Pour trouver le nombre de points à enlever d'une liste d'intensitées pour un pic donné
def GetNRemove(PeakIntensities, ConfInterval):
    remove = 0
    for n_remove in range(0, len(PeakIntensities)-2, 1):
        # n_remove va tester le nombre de raies au début des intensitées
        firstshots = PeakIntensities[0:n_remove+1]
        restshots = PeakIntensities[n_remove+1:-1]
            
        # Test statistique
        _, pnorm = mannwhitneyu(firstshots, restshots, method="asymptotic")
            
        # Voir si p est inférieur à l'intervalle de confiance attendu
        if pnorm < (1-ConfInterval):
            # À partir de maintenant, les deux moyennes sont jugées égales, donc on doit enlever juste à partir de l'ancien
            remove = n_remove-1

            # Sortir de la boucle plus rapidement
            break
            
        # Quoi faire si toujours p-value trop haute
        if n_remove == (len(PeakIntensities) - 1):
            remove = n_remove + 1
            
    return remove
