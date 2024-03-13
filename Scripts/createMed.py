class sim_data:
    # Classe para armazenar resultados de cálculo do estimador de estado    
    nvar = 0
    nmed = 0
    
    def __init__(self, df_DREF, x, Cov_X, z, residuo):
        self.df_DREF = df_DREF
        self.x = x
        self.Cov_X = Cov_X
        self.z = z
        self.residuo = residuo

class network_data:
    # Classe para salvar dados da rede elétrica no formato de data_frame de leitura
    n_bus = 0
    n_branch = 0
        
    def __init__(self, df_DBAR, df_DLIN, df_DREG, df_DSHNT, df_DTRF):
        self.df_DBAR = df_DBAR
        self.df_DLIN = df_DLIN
        self.df_DREG = df_DREG
        self.df_DSHNT = df_DSHNT
        self.df_DTRF = df_DTRF


def PrintConfig(md,sd):
    # Imprime arquivo de configuração do estimador trifásico
    # md: pasta principal com executável e arquivos de saída
    # sd: pasta do sistema a ser simulado
    file = open(md + '/config.txt','w')
    file.write(md + '\n') 
    file.write(sd)
    file.close()
    

def LeituraDados(foldername):
    #----------------------------------------------------------------
    #Função de leitura de dados da rede elétrica e préprocessamento
    # foldername: pasta com arquivos a serem lidos
    
    # Dados de Barras
    filename = '/DBAR.csv'
    try:
        df_DBAR = pd.read_csv(foldername + filename, sep = ',',header=None)
    except:
        print("DBAR não encontrado!!!")
        df_DBAR = pd.DataFrame(columns = range(1,18))
    df_DBAR.columns = ["ID", "LIGACAO", "FASES", "TENSAO_NOM", "PNOM_A", "PNOM_B", "PNOM_C", "QNOM_A", "QNOM_B", "QNOM_C", \
            "ZIP", "VNOM_A", "VNOM_B", "VNOM_C", "ANG_VNOM_A", "ANG_VNOM_B", "ANG_VNOM_C"]
    
    # Dados de Ramais e Circuitos
    filename = '/DLIN.csv'
    try:
        df_DLIN = pd.read_csv(foldername + filename, sep = ',',header=None)
    except:
        df_DLIN = pd.DataFrame(columns = range(1,23))
    df_DLIN.columns = ["DE", "PARA", "FASES", "COMPRIMENTO", "Raa","Xaa", "Rab","Xab", "Rac","Xac", "Rbb","Xbb",\
            "Rbc","Xbc", "Rcc","Xcc",\
            "Baa", "Bab", "Bac", "Bbb", "Bbc", "Bcc" ]
    
    # Dados de Reguladores de Tensão
    filename = '/DREG.csv'
    try:
        df_DREG = pd.read_csv(foldername + filename, sep = ',',header=None)
    except:
        df_DREG = pd.DataFrame(columns = range(1,27))    
    df_DREG.columns = ["DE", "PARA", "FASES", "VNOM", "REGULACAO", "NTAPS", "S_NOMINAL",
                        "RESISTENCIA","REATANCIA","LIGACAO","RELACAO_TP","RELACAO_TC","DELTA_V","R1","X1","R2","X2","R3",
                        "X3","V1","V2","V3","TIPOCONT","TAP1","TAP2","TAP3"]
    
    # Dados de Bancos de Capacitor
    filename = '/DSHNT.csv'
    try:
        df_DSHNT = pd.read_csv(foldername + filename, sep = ',',header=None)
    except:
        df_DSHNT = pd.DataFrame(columns = range(1,9))
    df_DSHNT.columns = ["ID", "LIGACAO", "FASES","TENSAO_NOM","QNOM_A","QNOM_B","QNOM_C","TipoCont"]
    
    # Dados de Transformadores
    filename = '/DTRF.csv'
    try:
        df_DTRF = pd.read_csv(foldername + filename, sep = ',',header=None)
    except:
        df_DTRF = pd.DataFrame(columns = range(1,14))
    df_DTRF.columns = ["DE", "PARA", "FASES", "VPRI", "VSEC", "S_NOMINAL", "RESISTENCIA", "REATANCIA",\
            "LIGACAO_PRI", "LIGACAO_SEC", "DEFASAMENTO", "TAP_PRI", "TAP_SEC"]
    
    network_model = network_data(df_DBAR, df_DLIN, df_DREG, df_DSHNT, df_DTRF)
    network_model.n_bus = len(network_model.df_DBAR)
    network_model.n_branch = len(network_model.df_DLIN) + len(network_model.df_DREG) + len(network_model.df_DTRF)   
  
    
    return network_model

def ExportSystemData(md, sd, network_model):
    #----------------------------------------------------------------
    #Função que exporta arquivos da rede elétrica
    # network_model: classe com dataframes da rede elétrica
    
    network_model.df_DBAR.to_csv(md + sd + "\DBAR.csv", sep = ',', index=False, header=False, float_format='%.15f')
    network_model.df_DSHNT.to_csv(md + sd + "\DSHNT.csv", sep = ',', index=False, header=False, float_format='%.15f')
    network_model.df_DLIN.to_csv(md + sd + "\DLIN.csv", sep = ',', index=False, header=False, float_format='%.15f')
    network_model.df_DREG.to_csv(md + sd + "\DREG.csv", sep = ',', index=False, header=False, float_format='%.15f')
    network_model.df_DSWTC.to_csv(md + sd + "\DSWTC.csv", sep = ',', index=False, header=False, float_format='%.15f')
    network_model.df_DTRF.to_csv(md + sd + "\DTRF.csv", sep = ',', index=False, header=False, float_format='%.15f')
    

def ExportMeasurementSet(md, sd, filename, mode, df_DREF, locMed):
    #----------------------------------------------------------------
    #Função que exporta em arquvio DMED.csv o plano de medição e respecitvos valores medidos
    # 
    # filename: nome do arquivo a ser salvo
    # mode: modo de escrita do arquivo
    # df_DMED: Valores de referência do plano de medição
    # locMed: plano de medição, local de instalação de cada medidor e respectivo tipo
    
    if pmu_polar == 1:
        # Tipos para PMUs em coordenadas polares
        tipos = {'FPQ': [0 , 1],    
                  'IPQ': [2 , 3],
                  'V': [4],
                  'ICur': [],
                  'FCur': [6 , 7],
                  'Vp': [4 , 5]} 
    elif pmu_polar == 0:
        # Tipos para PMUs em coordenadas retangulares
        tipos = {'FPQ': [0 , 1],    
                 'IPQ': [2 , 3],
                 'V': [4],
                 'ICur': [10 , 11],
                 'FCur': [12 , 13],
                 'Vp': [8, 9]} 
    
    
    df_DMED = pd.DataFrame()
    for key in locMed:
        if len(locMed[key]) > 0:
            
            # Testa se valor do localizador é tupla (medidas em ramos) ou int (medidas em barras)
            if type(locMed[key][0]) == tuple:
                df_Aux = df_DREF[df_DREF[['De', 'Para']].apply(tuple, axis=1).isin(locMed[key])]
                df_DMED = df_DMED.append(df_Aux[df_Aux.Tipo.isin(tipos[key])], ignore_index=True)
                
            else: 
                
                df_Aux = df_DREF[df_DREF.De.isin(locMed[key])]
                df_DMED = df_DMED.append(df_Aux[df_Aux.Tipo.isin(tipos[key])], ignore_index=True)
    
    df_DMED.to_csv(md + sd + filename, sep = ',', index=False, header=False, float_format='%.15f', mode=mode)

def PowerFlow(md, sd, network_model, loading, method):
    #----------------------------------------------------------------
    #Função que calcula o fluxo de potência para a rede elétrica
    # Chama rotina externa para o cálculo
    # network_model: modelo da rede elétrica em dataframe
    # loading: cenário de carga em dataframe
    # method: esolha do método 1= Newthon Rapshon e 2= Varredura Direta/Inversa (futuro)
    
    # Plano de medição para fluxo de carga - sem redundância seguindo modelo de tipo de barra PQ, PV ou VTeta
    locMed_PF = {'IPQ': network_model.df_DBAR['ID'][1:network_model.n_bus].to_list(),
                'FPQ': [],
                'V': [network_model.df_DBAR['ID'][0]]} 
         
    #Exporta condicao de carga como plano de medição sem redundância
    ExportMeasurementSet(md, sd, "/DMED.csv", 'w', loading, locMed_PF)
    
    #Calula fluxo de potência pelo método de newton
    #subprocess.check_call('./powerflow', cwd = md, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.check_output('../powerflow')
    
    #Leitura do resultado
    filename = 'refPowerflow.txt'
    #df_DSIM = pd.read_csv(md + filename, sep = ',',header=None)
    df_DSIM = pd.read_csv(filename, sep = ',',header=None)
    df_DSIM.columns = ['Estado','Tipo','De','Para','Fases','Zmed','Sigma']
    
    filename = 'state.txt'
    x = pd.read_csv(filename, sep = '\t',header=None)
    x.columns = ['regua','val']
    
    
    Cov_X = []
    
    residuo = []
    
    z = []
        
    #Salva resultado
    sim_simul = sim_data(df_DSIM, x, Cov_X, z, residuo)
    sim_simul.nvar = len(x)
    sim_simul.nmed = len(z)
    
    return sim_simul   

def SampleMeasurementsMC(md, sd, filename, mode, df_DREF, locMed, precision):
    #----------------------------------------------------------------
    #Função que exporta em arquvio DMED.csv o plano de medição e respecitvos valores medidos com inserção de ruído
    # 
    # filename: nome do arquivo a ser salvo
    # mode: modo de escrita do arquivo
    # df_DMED: Valores de referência do plano de medição
    # locMed: plano de medição, local de instalação de cada medidor e respectivo tipo
    
    nmed = len(df_DREF)
    
    # Inseri ruído aleatório nos valores de referência
    df_DMED = df_DREF.copy()
    df_DMED['Zmed'] = df_DREF['Zmed'].values + precision * abs(df_DREF['Zmed'].values) / 3 * np.random.randn(nmed)
    
    if precision == 0:
        precision = 0.00001
    df_DMED['Sigma'] = precision
    
    #Exporta condicao de carga como plano de medição sem redundância
    ExportMeasurementSet(md, sd, filename, mode, df_DMED, locMed)

def StateEstimation(md, sd, network_model, measurement_set, method):
    #----------------------------------------------------------------
    #Função que roda o estimador de estado para um
    # Chama rotina externa para o cálculo
    # network_model: modelo da rede elétrica em dataframe
    # measurement_set: vetor de medidas em dataframe
    # method: esolha do método 1= Newthon Rapshon
    
           
    
    #Roda estimador de estado
    #subprocess.check_call('./estimator', cwd = md, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.check_output('../estimator')
    #Leitura do resultado
    filename = 'referencia.txt'
    #df_DSIM = pd.read_csv(md + filename, sep = ',',header=None)
    df_DSIM = pd.read_csv(filename, sep = ',',header=None)
    df_DSIM.columns = ['Estado','Tipo','De','Para','Fases','Zmed','Sigma']
    
    filename = 'state.txt'
    x = pd.read_csv(filename, sep = '\t',header=None)
    x.columns = ['regua','val']
    print('state estimation', len(x))
    Cov_X = []
    
    # filename = '/residuoNormalizado.txt'
    # residuo = pd.read_csv(md + filename, sep = ',',header=None)
    # residuo.columns = ['id','r','rN','ec','UI']
    residuo = []
    
    z = []    
    
    #Salva resultado
    sim_simul = sim_data(df_DSIM, x, Cov_X, z, residuo)
    sim_simul.nvar = len(x)
    sim_simul.nmed = len(z)
    
    return sim_simul