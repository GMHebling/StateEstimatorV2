
import subprocess
import pandas as pd

def writeConfigFile(folder, system, estimator, idsim, order, flagPowerflow):
    destFile = "../config.txt"
    with open(destFile, 'w') as file:
        file.write(folder)
        file.write('\n')
        file.write(system)
        file.write('\n')
        file.write(estimator)
        file.write('\n')
        file.write(idsim)
        file.write('\n')
        file.write(order)
        file.write('\n')
        file.write(flagPowerflow)

def generateID():
    idFile = "./lastID.txt"
    lastID = 0
    with open(idFile) as file:
        for line in file:
            lastID = int(line)
    
    currID = lastID + 1
    with open(idFile, 'w') as file:
        file.write(str(currID))
    return str(currID)

def readSimList():
    fileSim = 'simulations.txt'
    simList = []
    with open(fileSim) as file:
        for line in file:
            lineS = line.split(',')
            currSys = int(lineS[0])
            currEst = int(lineS[1])
            currOrd = int(lineS[2])
            nRepeats = int(lineS[3])
            flgDMED = int(lineS[4])
            flgPowerflow = int(lines[5])
            simList.append((currSys, currEst, currOrd, nRepeats, flgDMED, flgPowerflow))
    return simList

def writeSimGuide (system, currID, _est, _ord, flgPowerflow):
    destFile = "./SimGuide.txt"
    with open(destFile, 'a') as file:
        file.write(currID)
        file.write(',')
        file.write(system)
        file.write(',')
        file.write(_est)
        file.write(', ')
        file.write(_ord)
        file.write('\n')
        file.write(flgPowerflow)
        file.write('\n')
        
def generateDMED():
    filename = '../src/refPowerflow.txt'
    df_DSIM = pd.read_csv(filename, sep = ',',header=None)
    df_DSIM.columns = ['Estado','Tipo','De','Para','Fases','Zmed','Sigma']
    pass

def runSimulations(simList):
    listSystem = ['/IEEE34', '/IEEE123', '/IEEE342SIM', '/IEEE906', '/REAL1058']
    listEstimator = [0, 1, 2, 3]
    listOrder = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    stdFolder = '../SimData'
    
    for item in simList:
        _sys = item[0]
        _est = str(item[1])
        _ord = str(item[2])
        nRepeats = item[3]
        flgDMED = item[4]
        flgPowerflow = item[5]
        if flgDMED == 1:
            #TODO: update DMED according to measurement set
            currID = generateID()
            writeConfigFile(stdFolder, listSystem[_sys], _est, currID, _ord, 1)
            writeSimGuide(listSystem[_sys], currID, _est, _ord, 1)
            subprocess.check_output('../src/testEstimator')
            #Create measurement from reference
            #save updated DMED file

            for i in range(nRepeats):
                    currID = generateID()
                    writeConfigFile(stdFolder, listSystem[_sys], _est, currID, _ord, flgPowerflow)
                    writeSimGuide(listSystem[_sys], currID, _est, _ord, flgPowerflow)
                    generateDMED()
                    subprocess.check_output('../src/testEstimator')
        else:
            for i in range(nRepeats):
                    currID = generateID()
                    writeConfigFile(stdFolder, listSystem[_sys], _est, currID, _ord, flgPowerflow)
                    writeSimGuide(listSystem[_sys], currID, _est, _ord, flgPowerflow)
                    subprocess.check_output('../src/testEstimator')




if __name__ == '__main__':
    simList = readSimList()
    runSimulations(simList)
    

    



