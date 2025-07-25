import numpy as np
import pandas as pd
import os

# Passer de la liste de str en liste de floats
def readlibslines(lines):
    clean = []
    repeat = []
    nbshots = []
    for line in lines:
        newline = [float(a) for a in line.rstrip().split()]
        nbshots.append(len(newline)-1)

        wavelength = newline[0]
        for value in newline[1:]:
            repeat.append([wavelength, value])

        avg = np.mean(newline[1:])
        stdev = 2*np.std(newline[1:])
        newline.append(avg)
        newline.append(avg-stdev)
        newline.append(avg+stdev)
        clean.append(newline)

    # Tester tous meme nombre shots
    if all(x == nbshots[0] for x in nbshots):
        nshots = nbshots[0]
    else:
        nshots = 0
        print('Pas le même nombre de tirs partout!')

    # Putting it all in a clean cute dataframe
    clean = pd.DataFrame(clean)

    # Identify boundaries of spectrometers
    signal = np.array(clean[0])
    boundaries = list(np.where(np.abs(np.diff(signal)) > 1)[0])
    spec_lambda = []
    for i in range(len(boundaries)):
        if i == 0:
            specbound = [signal[0], signal[boundaries[i]]]
            spec_lambda.append(specbound)
        else:
            specbound = [signal[boundaries[i-1]+1], signal[boundaries[i]]]
            spec_lambda.append(specbound)

            if i == len(boundaries)-1:
                lastbound = [signal[boundaries[i]+1], signal[-1]]
                spec_lambda.append(lastbound)

    # Sorting the clean spectra
    clean = clean.sort_values(0)

    return clean, pd.DataFrame(repeat), nshots, spec_lambda


# Connaître les noms des spectres
def GrabLIBSSpectra(folder, common_name, terminaison):
    # Move to the folder where the files are located
    os.chdir(folder)

    # Grab the stuff
    documents = os.listdir(os.getcwd())
    fichiers = []
    names = []
    for text in documents:
        if text.startswith(common_name):
            if text.endswith(terminaison):
                fichiers.append(text)
                names.append(text.rstrip(terminaison))
                
    # Move back in the main folder
    os.chdir('..')            

    return fichiers, names


# Read content of the libs spectra file
def ReadLIBSSpectra(folder, fichiers):
    # Go in directory with files
    os.chdir(folder)
    
    # Read the files
    LIBS = []
    LIBSrepeat = []
    Nshots = []
    spec_lambda = []
    for filename in fichiers:
        with open(filename, 'r') as fID:
            lines = fID.readlines()

        mean, repeats, nshots, spec_lambda = readlibslines(lines)

        LIBS.append(mean)
        LIBSrepeat.append(repeats)
        Nshots.append(nshots)

    # Go back to main
    os.chdir('..')    

    return LIBS, LIBSrepeat, Nshots, spec_lambda