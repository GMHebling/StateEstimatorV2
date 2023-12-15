
import subprocess

def writeConfigFile(folder, system, estimator, idsim, order):
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
            simList.append((currSys, currEst, currOrd, nRepeats, flgDMED))
    return simList

def writeSimGuide (system, currID, _est, _ord):
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
        for i in range(nRepeats):
            if flgDMED == 1:
                pass
            else:
                currID = generateID()
                writeConfigFile(stdFolder, listSystem[_sys], _est, currID, _ord)
                writeSimGuide(listSystem[_sys], currID, _est, _ord)
                subprocess.check_output('../src/testEstimator')




if __name__ == '__main__':
    simList = readSimList()
    runSimulations(simList)
    

    



