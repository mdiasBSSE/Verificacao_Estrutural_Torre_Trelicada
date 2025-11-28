import numpy as np
import csv
from Wind.Wind_Ec import Zone

def Model_csv_read(Model_csv_file):
    with open(Model_csv_file, newline='', encoding='utf-8') as f:
        leitor = list(csv.reader(f))  # transforma em lista de listas
        
    # Procurar índices dos marcadores
    for ilinha, linha in enumerate(leitor):
        if linha and linha[0] == "Nodes":
            iNodes = ilinha+1
        elif linha and linha[0] == "Lines":
            iElements = ilinha+1
        elif linha and linha[0] == "Diag_Hor_Connection":
            i_Diag_Hor = ilinha+1
        elif linha and linha[0] == "Shaft Connection - Flange":
            i_Shaft_Flange = ilinha+1
        elif linha and linha[0] == "Shaft Connection - Bolted":
            i_Shaft_Bolted = ilinha+1
        elif linha and linha[0]=="Flange Foundation Angle":
            i_Flange_Foundation_Angle=ilinha+1
        elif linha and linha[0] == "Foundation":
            i_Foundation = ilinha+1    
        elif linha and linha[0] == "Orientation":
            i_Orientation = ilinha+1


    # Criar as matrizes entre os marcadores
    Nodes_Matrix        = np.array(leitor[iNodes+1 : iElements-2])          # de Nodes até Elements
    Elements_Matrix     = np.array(leitor[iElements+1 : i_Diag_Hor-2])      # de Elements até Diag_Hor
    Dig_Hor_Matrix      = np.array( leitor[i_Diag_Hor+1 : i_Shaft_Flange-2])  # de Diag_Hor até Shaft Flange
    Shaft_Flange_Matrix = np.array( leitor[i_Shaft_Flange+1 : i_Shaft_Bolted-2])  # de Shaft Flange até Shaft Bolted
    Shaft_Bolted_Matrix = np.array( leitor[i_Shaft_Bolted+1 :i_Flange_Foundation_Angle-2])        # de Shaft Bolted até Foundation
    Flange_Foundation_Angle = np.array( leitor[i_Flange_Foundation_Angle+1:i_Foundation-2])
    Foundation_Matrix   = np.array( leitor[i_Foundation :i_Orientation-2])           #de Foundation até Orientation
    Alpha_North=float(leitor[i_Orientation+1][0])
    
    return Nodes_Matrix  ,Elements_Matrix  ,Dig_Hor_Matrix ,Shaft_Flange_Matrix,Shaft_Bolted_Matrix,Flange_Foundation_Angle,Foundation_Matrix,Alpha_North


def Structure_txt_read(txt_file):
    # with open(txt_file, 'r') as file:
    #     linhas = [linha.strip() for linha in file.readlines()]
    # Pais=linhas[0]
    # Alt=0
    # FundMethod=False
    # Terrain=linhas[1]
    # Zona=linhas[2]
    # Classe_Fiabilidade=int(linhas[3])
    # Gelo=linhas[4]
    # Altitude=float(linhas[5])
    # C0_calc=linhas[6]
    # Tipoc0=linhas[7]
    # Alt_col=float(linhas[8])
    # Lu=float(linhas[9])
    # Ld=float(linhas[10])
    # Xtopo=float(linhas[11])
    # Ac_c0=float(linhas[12])
    # A500=float(linhas[13])
    # A1000=float(linhas[14])
    # for i in range(len(linhas)):
    #     if (linhas[i]=="Fund:\n" or linhas[i]=="Fund:") and len(linhas)==i+6:
    #         FundMethod=True
    #         MfundMax=np.array(linhas[i+1].split(";"), dtype=float)
    #         RatioFundCalc=np.array(linhas[i+2].split(";"), dtype=float)
    #         SoilSelfWeight = float(linhas[i+3].split(";")[1])
    #         AllowedSoilTension_SLS = float(linhas[i+4].split(";")[1])
    #         AllowedSoilTension_ULS = float(linhas[i+5].split(";")[1])
    # if FundMethod==False:
    #         MfundMax=0
    #         RatioFundCalc=0
    #         SoilSelfWeight = 0
    #         AllowedSoilTension_SLS = 0
    #         AllowedSoilTension_ULS = 0
    with open(txt_file, "r", encoding="utf-8") as f:
        for linha in f:
            partes = linha.strip().split(",")
        Pais=partes[0]
        Calc_Wind_Auto=partes[1]
        if Calc_Wind_Auto=="Yes":
            z0=float(partes[2])
            zmin=float(partes[3])
            vb=float(partes[4])
            Classe_Fiabilidade=int(partes[5])
        else:
            Terrain=partes[2]
            vb=float(partes[3])
            _,z0,zmin=Zone("-",Terrain,Pais)
            Classe_Fiabilidade=int(partes[5])
    Alt=0
    Gelo="No"
    Altitude=0
    C0_calc="No"
    Tipoc0=0
    Alt_col=0
    Lu=0
    Ld=0
    Xtopo=0
    Ac_c0=0
    A500=0
    A1000=0
    FundMethod=False
    MfundMax=0
    RatioFundCalc=0
    SoilSelfWeight = 0
    AllowedSoilTension_SLS = 0
    AllowedSoilTension_ULS = 0
    return Pais,vb,z0,zmin,Classe_Fiabilidade,Gelo,Altitude,C0_calc,Tipoc0,Alt_col,Lu,Ld,Xtopo,Ac_c0,A500,A1000,Alt,FundMethod,MfundMax,RatioFundCalc,SoilSelfWeight,AllowedSoilTension_SLS,AllowedSoilTension_ULS


def Ant_txt_read(ant_txt_file):
    Nanalise=0
    Org=0
    with open(ant_txt_file, 'r') as file:
        Antenas = file.readlines()
    for i in range(len(Antenas)):
        if Antenas[i]=="Final;1\n" or Antenas[i]=="Final;0\n" or Antenas[i]=="Final;1" or Antenas[i]=="Final;0":
            Finalraw=Antenas[i]

            if Antenas[i]=="Final;1\n":
                Analise=(Antenas[i+1].strip()).split(';')
                Nanalise=Analise[1]
                Org=Antenas[i+3:-1]
    Final=(Finalraw.strip()).split(';')
    Final=float(Final[1])

    iLE=-1
    iELE=-1

    for i in range(len(Antenas)):
        if Antenas[i]=="LE\n":
            iLE=i
        if Antenas[i]=="ELE\n":
            iELE=i
    CabosOutside=1
    Ladder=1
    Auto=1
    if iLE!=-1:
        CabosOutside=int(Antenas[iLE+1].strip())
        Ladder=int(Antenas[iLE+2].strip())
        Auto=int(Antenas[iLE+3].strip())
        equipline_manual=Antenas[iLE+4:iELE]
        Antenas=Antenas[:iLE]

    return Antenas,equipline_manual,Auto,Ladder,CabosOutside,Final,Nanalise,Org