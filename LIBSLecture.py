# Camille Rincent
# Lecture des spectres de LIBS en format .txt vers graphique utilisable pour analyses.

from Analysis import *
from Inputs import *
from Outputs import *


# Début du code
files = "GaInSn_"
terminaison = ".txt"

fichiers, names = GrabLIBSSpectra(files, terminaison)
LIBS, LIBSrepeat, Nshots, spec_lambda = ReadLIBSSpectra(fichiers)
AvgLIBS = GetAvgLIBS(LIBS, Nshots, False, True)



# Identifier les pics
# Séparer le spectre en bin, chaque spectro ou overlap ayant son propre bin
bins = IdentifyAndSplitOverlaps(spec_lambda)

# Solution moche pour les thresholds selon le bruit de chaque spectrometre sur moyenne
#thresholds = [15, 100, 250, 15]
thresholds = [15, 300, 700, 50]

# Analyze for peaks and give table of all peaks and their wavelengths

# 2 will show peak by peak what is used
# 3 for the detail in the identification of the peaks, keep empty
# Keep list empty if no graph should be shown to show methodology
peakgraphtoshow = []
PeakIds, PeakIs = PeakLambdaAndIntensityCompute(AvgLIBS, bins, thresholds, peakgraphtoshow)

# Identifier les intensités en fonction du #tir pour chaque pic identifié
NTirsRemove = IdentifyFirstShotsDifference(LIBS, PeakIs, 0.95, True)


# Graph proof of concept plusieurs spectres ensembles
#PrintSpectrasSubplots(LIBS, Nshots, names)



print('terminé')