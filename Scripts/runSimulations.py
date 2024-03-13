
import subprocess
import pandas as pd
import numpy as np

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
            flgPowerflow = int(lineS[5])
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
        
def generateDMED(system):
    filename = '../src/refPowerflow.txt'
    df_DREF = pd.read_csv(filename, sep = ',',header=None)
    df_DREF.columns = ['Estado','Tipo','De','Para','Fases','Zmed','Sigma']

    filename = '../SimData/' + system + '/DMED_BASE.csv'
    df_DMED_base = pd.read_csv(filename, sep = ',',header=None)
    df_DMED_base.columns = ['LIGADO','TIPO','DE','PARA','FASES','ZMED','SIGMA']

    dict_DMED = {'LIGADO': [],'TIPO': [],'DE': [],'PARA': [],'FASES': [],'ZMED': [],'SIGMA': []}
    for idx, row in df_DMED_base.iterrows():
        tipo = row['TIPO']
        de = row['DE']
        para = row['PARA']
        fases = row['FASES']

        z_ref = get_reference_value(df_DREF, tipo, de, para, fases)
        #precision = get_precision(z_med)
        error = calculate_meas_error(z_ref, 0.0005)
        
        z_update = z_ref + error
        
        dict_DMED['LIGADO'].append(int(row['LIGADO']))
        dict_DMED['TIPO'].append(int(row['TIPO']))
        dict_DMED['DE'].append(int(row['DE']))
        dict_DMED['PARA'].append(int(row['PARA']))
        dict_DMED['FASES'].append(int(row['FASES']))
        dict_DMED['ZMED'].append(z_update)
        dict_DMED['SIGMA'].append(row['SIGMA'])

    df_DMED_update = pd.DataFrame(dict_DMED)
    df_DMED_update.to_csv('../SimData' + system + '/DMED.csv', index=False, header=False)
    pass

def get_reference_value(DREF, tipo, de, para, fases):
    z_ref = DREF[(DREF['Tipo'] == tipo) & (DREF['De'] == de) & (DREF['Para'] == para) & (DREF['Fases'] == fases)]['Zmed'].values
    if len(z_ref) == 0:
        return None
    elif len(z_ref) > 1:
        return z_ref[0]
    else:
        return z_ref[0]
    
def calculate_meas_error(zmed, precision):
    error = (precision * abs(zmed))/(3*np.random.randn())
    return error

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
        flgPowerflow = str(item[5])
        if flgDMED == 1:
            #TODO: update DMED according to measurement set
            currID = generateID()
            writeConfigFile(stdFolder, listSystem[_sys], _est, currID, _ord, '1')
            writeSimGuide(listSystem[_sys], currID, _est, _ord, '1')
            subprocess.check_output('./testEstimator', cwd='../src')
            #Create measurement from reference
            #save updated DMED file

            for i in range(nRepeats):
                    currID = generateID()
                    writeConfigFile(stdFolder, listSystem[_sys], _est, currID, _ord, flgPowerflow)
                    writeSimGuide(listSystem[_sys], currID, _est, _ord, flgPowerflow)
                    generateDMED(listSystem[_sys])
                    subprocess.check_output('./testEstimator', cwd='../src')
        else:
            for i in range(nRepeats):
                    currID = generateID()
                    writeConfigFile(stdFolder, listSystem[_sys], _est, currID, _ord, flgPowerflow)
                    writeSimGuide(listSystem[_sys], currID, _est, _ord, flgPowerflow)
                    subprocess.check_output('../src/testEstimator')




if __name__ == '__main__':
    simList = readSimList()
    runSimulations(simList)
    
    

    



