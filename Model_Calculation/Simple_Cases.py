import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # necessário para plot 3D
from openseespy import opensees as ops
from scipy.cluster.hierarchy import linkage, fcluster
import csv
import os
import time
from pathlib import Path
import sys
from Utilities.Utilities import encontra_vizinhos
from Utilities.Utilities_EC import ClasseFiabilidade,Equivalent_Mass
from Wind.Wind_Ec import CalcKa,calculoc0,Iv_calc,Velocidade_media,Pressure_wind,Coef_forca,CsCd_calc,Zone
from Model_Calculation.Calculation_Methods_ops import Create_Load_Case,Calc_Load_Case
from Model_Calculation.Load_Creation_Area import Area_Definition,Area_Position_Mass
from Output.Output_Lattice_Expanded import Shaft_Wind
from Output.Output_Lattice_Reduced import Wind_Ice_Out
from Wind.Wind_Lattice import Elements_Areas
def Simple_Cases(Ncasoscarga,nos_por_classe,C0_calc,Pais,zmax,Alt,Tipoc0,Alt_col,Lu,Ld,Xtopo,Ac_c0,A500,A1000,rho_ar,Dim,visco_ar,nos_vento,Troco,TrussType,z_vento,ElementosID,ExpVento,Comprimento_Barra,divisor,esp_gelo,Gelo,h_torre,munit,Massa_Extra,todos_nos,nEle,nos_por_elemento,coords,h_vetor,z_grupo,nos_grupo,dh,torre,dalpha,Posicao_Antenna,Massa_Antenna,AreaDirection,equipline,Ant,iStart,rho_gelo,vb,z0,zmin,Calc_Wind_Auto):    

    A_lin_Troco=np.zeros(len(nos_por_classe)+1)
    A_lin_Troco_Cf=np.zeros(len(nos_por_classe)+1)
    L_lin_equiv=np.zeros(len(nos_por_classe)+1)
    Z_med_vento=np.zeros(len(nos_por_classe)+1)
    if C0_calc=="Sim":
        c0_no_maximo=calculoc0(Pais,zmin,zmax+Alt,Tipoc0,Alt_col,Lu,Ld,Xtopo,Ac_c0,A500,A1000)
    else:
        c0_no_maximo=np.ones(len(zmax+Alt))
    Pressao_no_maximo=Pressure_wind(zmax+Alt,vb,z0,zmin,Alt,Pais,rho_ar,c0_no_maximo)
    v_no_maximo=np.sqrt(Pressao_no_maximo*2/rho_ar)
    Reynolds=v_no_maximo*Dim/visco_ar

    nos_vento_troco,Largura,Largura_total,Largura_media,At,Af,Ac,Acsup,As,Indice_Cheios=Elements_Areas(nos_por_classe,nos_vento,Troco,TrussType,Dim,z_vento,ElementosID,ExpVento,Comprimento_Barra,divisor,Reynolds,esp_gelo)


    for i in range(len(z_vento)-1):
        Z_med_vento[i+1]=(z_vento[i]+z_vento[i+1])/2
    if C0_calc=="Sim":
        c0_troco_medio=calculoc0(Pais,zmin,Z_med_vento+Alt,Tipoc0,Alt_col,Lu,Ld,Xtopo,Ac_c0,A500,A1000)
    else:
        c0_troco_medio=np.ones(len(Z_med_vento+Alt))

    Pressao_troco_medio=Pressure_wind(Z_med_vento+Alt,vb,z0,zmin,Alt,Pais,rho_ar,c0_troco_medio)

    Mascara= (TrussType == "Leg")
    Contagem_Montante_Troco = np.bincount(Troco[Mascara].astype(int), minlength=int(np.max(Troco)+1)) # ignora índice 0
    if C0_calc=="Sim":
        c0_Zs=calculoc0(Pais,zmin,np.array([0.6*h_torre+Alt]),Tipoc0,Alt_col,Lu,Ld,Xtopo,Ac_c0,A500,A1000)
    else:
        c0_Zs=np.ones(1)
    vm=Velocidade_media(np.array([0.6*h_torre+Alt]),vb,z0,zmin,Pais,c0_Zs)
    Iv_Zs=Iv_calc(np.array([0.6*h_torre])+Alt,z0,zmin,Alt,Pais,c0_Zs)
    if ExpVento[int(Mascara[0])]=="Circ":
        Leg_Type="Circ"
    else:
        Leg_Type="Flat"



    Massa_Troco=np.zeros(len(nos_por_classe)+1)
    ###Calculo me

    Massa_Barra=Comprimento_Barra*munit
    for i in range(len(nos_por_classe)+1):
        Massa_Troco[i] = Massa_Barra[Troco == i].sum()
    Massa_Troco=Massa_Troco*(1+Massa_Extra/100)


    Zs=0.6*h_torre
    iZs=encontra_vizinhos(z_vento[::-1], Zs)
    Largura_ref=Largura[:,iZs[0]+1]
    Larg=np.zeros(3)
    Larg[0]=Largura[1,-1]
    Larg[2]=Largura[0,1]
    if Largura_ref[0]==Largura_ref[1]:
        Larg[1]=Largura_ref[0]
    else:
        h1=z_vento[iZs[0]]
        h2=z_vento[iZs[1]]
        L1=Largura_ref[0]
        L2=Largura_ref[1]
        Larg[1]=L2+(Zs-h2)*(L1-L2)/(h1-h2)



    cfs=np.zeros((Ncasoscarga,len(nos_por_classe)+1))
    FmW=np.zeros((Ncasoscarga,len(nos_por_classe)+1))
    FTW=np.zeros((Ncasoscarga,len(nos_por_classe)+1))
    FTW_lin=np.zeros((Ncasoscarga,len(nos_por_classe)+1))
    FTW_dist=np.zeros((Ncasoscarga,len(nos_por_classe)+1))
    FTW_lin_dist=np.zeros((Ncasoscarga,len(nos_por_classe)+1))
    F_Troco=np.zeros((Ncasoscarga,len(nos_por_classe)+1))
    F_nos=np.zeros((Ncasoscarga,len(todos_nos)+2,6))
    Reactions=np.zeros((Ncasoscarga,len(todos_nos)+2,6))
    Massa_nos=np.zeros(len(todos_nos)+2)
    Massa_nos_lin=np.zeros(len(todos_nos)+2)
    F=np.zeros((Ncasoscarga,nEle+1,12))
    F_axial=np.zeros((Ncasoscarga,nEle+1))
    Desloc=np.zeros((Ncasoscarga,len(todos_nos)+2,6))
    MzVer=np.zeros(Ncasoscarga)
    Ry=np.zeros(Ncasoscarga)
    Desloc_100=np.zeros((Ncasoscarga,len(todos_nos)+2,6))
    Desloc_120=np.zeros((Ncasoscarga,len(todos_nos)+2,6))
    Desloc_aux=np.zeros((Ncasoscarga,len(todos_nos)+2))
    CsCd=np.zeros((Ncasoscarga))
    NameCaso=np.zeros((Ncasoscarga)).astype(str)
    #######Caso de carga 1 - Peso Proprio #################



    NameCaso[1]= Create_Load_Case(1+iStart-1,"Peso Proprio")

    for ele, nos in nos_por_elemento.items():
        ops.load(nos[0],0,0,-Massa_Barra[ele]*9.81/2*(1+Massa_Extra/100),0,0,0)
        ops.load(nos[1],0,0,-Massa_Barra[ele]*9.81/2*(1+Massa_Extra/100),0,0,0)
        F_nos[1,int(nos[0]),2]=F_nos[1,int(nos[0]),2]-Massa_Barra[ele]*9.81/2*(1+Massa_Extra/100)
        F_nos[1,int(nos[1]),2]=F_nos[1,int(nos[1]),2]-Massa_Barra[ele]*9.81/2*(1+Massa_Extra/100)

    Desloc[1,:,:],Reactions[1,:,:],Desloc_aux[1,:],F[1,:,:],F_axial[1,:]=Calc_Load_Case(1+iStart-1,coords,nos_por_elemento,todos_nos,nEle)
    # for tag, u in desloc.items():
    #     # arredonda cada componente a 5 casas
    #     u_rounded = [round(val, 5) for val in u]
    #     print(f"Desloc nó {tag}: {u_rounded}")

    # ops.reactions()
    # reaction_1 = ops.nodeReaction(1)
    # print(f"Reações no nó 1: Rx = {reaction_1[0]} N, Ry = {reaction_1[1]} N, Rz = {reaction_1[2]} N")
    # reaction_1 = ops.nodeReaction(25)
    # print(f"Reações no nó 25: Rx = {reaction_1[0]} N, Ry = {reaction_1[1]} N, Rz = {reaction_1[2]} N")
    # reaction_1 = ops.nodeReaction(37)
    # print(f"Reações no nó 37: Rx = {reaction_1[0]} N, Ry = {reaction_1[1]} N, Rz = {reaction_1[2]} N")


    #######Caso de carga Vento #################
    icasoscarga=1
    alpha_vector=np.arange(0, 360, dalpha)

    ivento=0 




    if C0_calc=="Sim":
        c0_Vento_Ant=calculoc0(Pais,zmin,Posicao_Antenna+Alt,Tipoc0,Alt_col,Lu,Ld,Xtopo,Ac_c0,A500,A1000)
    else:
        c0_Vento_Ant=np.ones(len(Posicao_Antenna+Alt))
    Pressao_Vento_Ant=Pressure_wind(Posicao_Antenna+Alt,vb,z0,zmin,Alt,Pais,rho_ar, c0_Vento_Ant)  
    F_Total_Ant=np.zeros((Ncasoscarga,AreaDirection.shape[0]))
    
    iAnt,x_Ant, Massa_nos,Massa_nos_lin=Area_Position_Mass(Posicao_Antenna,h_torre,z_grupo,nos_grupo,Massa_Antenna,divisor,equipline,TrussType,todos_nos,Gelo)

    for i in todos_nos:  
        ops.mass(int(i),Massa_nos[int(i)]+Massa_nos_lin[int(i)],Massa_nos[int(i)]+Massa_nos_lin[int(i)],Massa_nos[int(i)]+Massa_nos_lin[int(i)] ,0,0,0)
    # === Modelo modal
    num_modes = 10
    eigenvalues = ops.eigen(num_modes)  # retorna ω²

    # Converter para frequências e períodos
    omegas = np.sqrt(np.array(eigenvalues))
    freqs = omegas / (2*np.pi)
    # print("\n=== Resultados Modais ===")
    # for i in range(num_modes):
    #     print(f"Modo {i+1}: ω²={eigenvalues[i]:.4e}, f={freqs[i]:.4f} Hz")

    #### Caso 2 - Peso equipamentos

    icasoscarga=icasoscarga+1
    NameCaso[icasoscarga]= Create_Load_Case(icasoscarga+iStart-1,"Peso Equipamentos - Linear")

    for i in todos_nos:
        if i!=0:
            ops.load(int(i), 0,0, -Massa_nos_lin[i]*9.81,0,0,0)
            F_nos[icasoscarga,int(i),2]=F_nos[icasoscarga,int(i),2]-Massa_nos_lin[i]*9.81
            
    Desloc[icasoscarga,:,:],Reactions[icasoscarga,:,:],Desloc_aux[icasoscarga,:],F[icasoscarga,:,:],F_axial[icasoscarga,:]=Calc_Load_Case(icasoscarga+iStart-1,coords,nos_por_elemento,todos_nos,nEle)

    #### Caso 3 - Peso equipamentos

    icasoscarga=icasoscarga+1
    NameCaso[icasoscarga]= Create_Load_Case(icasoscarga+iStart-1,"Peso Equipamento - Pontual")
    for i in range(len(Posicao_Antenna)):
        for j in range(nos_grupo.shape[1]):
            ops.load(int(nos_grupo[int(iAnt[i,1]),j]), 0,0, -(1-x_Ant[i])*Massa_Antenna[i]/divisor*9.81,0,0,0)
            ops.load(int(nos_grupo[int(iAnt[i,0]),j]), 0,0,-x_Ant[i]*Massa_Antenna[i]/divisor*9.81,0,0,0)
            F_nos[icasoscarga,int(nos_grupo[int(iAnt[i,1]),j]),2]=F_nos[icasoscarga,int(nos_grupo[int(iAnt[i,1]),j]),2]-(1-x_Ant[i])*Massa_Antenna[i]/divisor*9.81
            F_nos[icasoscarga,int(nos_grupo[int(iAnt[i,0]),j]),2]=F_nos[icasoscarga,int(nos_grupo[int(iAnt[i,0]),j]),2]-x_Ant[i]*Massa_Antenna[i]/divisor*9.81
    Desloc[icasoscarga,:,:],Reactions[icasoscarga,:,:],Desloc_aux[icasoscarga,:],F[icasoscarga,:,:],F_axial[icasoscarga,:]=Calc_Load_Case(icasoscarga+iStart-1,coords,nos_por_elemento,todos_nos,nEle)

    try:
        Ant_max=max(Ant[:,0])
    except:
        Ant_max=0
    if Ant_max>h_torre:
        h_equip=np.arange(0,max(Ant[:,0])+dh,dh)
    else:
        h_equip=h_vetor
    me,equip_mass_loc,equip_mass_dist,m_vector=Equivalent_Mass(h_vetor,Posicao_Antenna,Ant,h_torre,Massa_Antenna,equipline,dh,z_vento,Massa_Troco,h_equip,"No")


    for i in range(len(z_vento)-1):
        for j in range(equipline.shape[0]):

            if z_vento[i+1]<=equipline[j,-1] and z_vento[i+1]>=equipline[j,-2]:
                if equipline[j,-1]>=z_vento[i]:
                    dL=min(z_vento[i]-z_vento[i+1],z_vento[i]-equipline[j,-2])
                else:
                    dL=equipline[j,-1]-z_vento[i+1]
                
                A_lin_Troco[i+1]=A_lin_Troco[i+1]+equipline[j,0]*dL
                A_lin_Troco_Cf[i+1]=A_lin_Troco_Cf[i+1]+equipline[j,0]*dL*equipline[j,1]
        L_lin_equiv[i+1]=A_lin_Troco[i+1]/(z_vento[i]-z_vento[i+1])



    Ka=CalcKa(A_lin_Troco,As,At,L_lin_equiv,Largura_media)
    z_vento_roll=np.roll(z_vento,1)
    FmW_lin=Pressao_troco_medio/(1+7*Iv_Zs)*A_lin_Troco_Cf*Ka
    if C0_calc=="Sim":
        c0_Vento_For=calculoc0(Pais,zmin,z_vento_roll+Alt,Tipoc0,Alt_col,Lu,Ld,Xtopo,Ac_c0,A500,A1000)
    else:
        c0_Vento_For=np.ones(len(z_vento_roll+Alt))
    ivento=icasoscarga
    Case_Out_csv=np.array([["Expanded Results","","","","","","","","","","","","","","","","","","","",""]])
    for i in range(len(alpha_vector)):
        
        ops.remove('loadPattern',int(i+icasoscarga))
        ops.wipeAnalysis()
        NameCaso[i+icasoscarga+1]= Create_Load_Case(i+icasoscarga+1+iStart-1,f"Wind Load - {alpha_vector[i]}")
        cf0f,cf0c,cf0csup,cfs0,cfs[i+icasoscarga+1,:]= Coef_forca(Af,Ac,Acsup,As,alpha_vector[i],Indice_Cheios,torre)
        CsCd[i+icasoscarga+1],CsCd_out=CsCd_calc(0.6*h_torre+Alt,h_torre,Iv_Zs,Larg,z0,zmin,freqs[0],me,cfs[i+icasoscarga+1,:],rho_ar,vm,Alt,Pais)
        FmW[i+icasoscarga+1,:]=Pressao_troco_medio/(1+7*Iv_Zs)*cfs[i+icasoscarga+1,:]*As
        FTW[i+icasoscarga+1,:]=FmW[i+icasoscarga+1,:]*(1+(1+0.2*(z_vento_roll/h_torre)**2)*((1+7*Iv_Zs)*CsCd[i+icasoscarga+1]-1)/c0_Vento_For)
        FTW_lin[i+icasoscarga+1,:]=FmW_lin*(1+(1+0.2*(z_vento_roll/h_torre)**2)*((1+7*Iv_Zs)*CsCd[i+icasoscarga+1]-1)/c0_Vento_For)
        for k in range(len(nos_vento)-1):
            Mascara = (Troco[1:] == int(k+1)) & (TrussType[1:] == "Leg")
            FTW_dist[i+icasoscarga+1,k+1]=FTW[i+icasoscarga+1,k+1]/divisor/(z_vento[k]-z_vento[k+1])
            FTW_lin_dist[i+icasoscarga+1,k+1]=FTW_lin[i+icasoscarga+1,k+1]/divisor/(z_vento[k]-z_vento[k+1])
            for ele in ElementosID[Mascara]:
                ops.eleLoad("-ele", int(ele), "-type", "beamUniform",-(FTW_dist[i+icasoscarga+1,k+1]+FTW_lin_dist[i+icasoscarga+1,k+1])*np.sin(np.radians(alpha_vector[i])),(FTW_dist[i+icasoscarga+1,k+1]+FTW_lin_dist[i+icasoscarga+1,k+1])*np.cos(np.radians(alpha_vector[i])),0)
                # nós do elemento
                nodes = ops.eleNodes(int(ele))
                # coordenada x de cada nó
                x_coords = [ops.nodeCoord(n)[2] for n in nodes]
                Ry[i+icasoscarga+1]=Ry[i+icasoscarga+1]+(FTW_dist[i+icasoscarga+1,k+1]+FTW_lin_dist[i+icasoscarga+1,k+1])*Comprimento_Barra[int(ele)]
                MzVer[i+icasoscarga+1]=MzVer[i+icasoscarga+1]+(FTW_dist[i+icasoscarga+1,k+1]+FTW_lin_dist[i+icasoscarga+1,k+1])*Comprimento_Barra[int(ele)]*np.mean(x_coords)
                #print(MzVer)
        
        F_Total_Ant[i+icasoscarga+1,:] = (Pressao_Vento_Ant * AreaDirection[:,i]) * CsCd[i+icasoscarga+1]
        for k in range(len(Posicao_Antenna)):
            for j in range(nos_grupo.shape[1]):
                ops.load(int(nos_grupo[int(iAnt[k,1]),j]), -(1-x_Ant[k])*F_Total_Ant[i+icasoscarga+1,k]*np.sin(np.radians(alpha_vector[i]))/divisor, (1-x_Ant[k])*F_Total_Ant[i+icasoscarga+1,k]*np.cos(np.radians(alpha_vector[i]))/divisor,0,0,0,0)
                ops.load(int(nos_grupo[int(iAnt[k,0]),j]), -x_Ant[k]*F_Total_Ant[i+icasoscarga+1,k]*np.sin(np.radians(alpha_vector[i]))/divisor , x_Ant[k]*F_Total_Ant[i+icasoscarga+1,k]*np.cos(np.radians(alpha_vector[i]))/divisor,0,0,0,0)
                F_nos[i+icasoscarga+1,int(nos_grupo[int(iAnt[k,1]),j]),0]=F_nos[i+icasoscarga+1,int(nos_grupo[int(iAnt[k,1]),j]),0] -(1-x_Ant[k])*F_Total_Ant[i+icasoscarga+1,k]*np.sin(np.radians(alpha_vector[i]))/divisor
                F_nos[i+icasoscarga+1,int(nos_grupo[int(iAnt[k,1]),j]),1]=F_nos[i+icasoscarga+1,int(nos_grupo[int(iAnt[k,1]),j]),1]+(1-x_Ant[k])*F_Total_Ant[i+icasoscarga+1,k]*np.cos(np.radians(alpha_vector[i]))/divisor
                F_nos[i+icasoscarga+1,int(nos_grupo[int(iAnt[k,0]),j]),0]=F_nos[i+icasoscarga+1,int(nos_grupo[int(iAnt[k,0]),j]),0]-x_Ant[k]*F_Total_Ant[i+icasoscarga+1,k]*np.sin(np.radians(alpha_vector[i]))/divisor
                F_nos[i+icasoscarga+1,int(nos_grupo[int(iAnt[k,0]),j]),1]=F_nos[i+icasoscarga+1,int(nos_grupo[int(iAnt[k,0]),j]),1]+x_Ant[k]*F_Total_Ant[i+icasoscarga+1,k]*np.cos(np.radians(alpha_vector[i]))/divisor
                x_coords = ops.nodeCoord(int(nos_grupo[int(iAnt[k,0]),j]))[2]
                MzVer[i+icasoscarga+1]=MzVer[i+icasoscarga+1]+x_Ant[k]*F_Total_Ant[i+icasoscarga+1,k]/divisor*x_coords
                x_coords = ops.nodeCoord(int(nos_grupo[int(iAnt[k,1]),j]))[2]
                MzVer[i+icasoscarga+1]=MzVer[i+icasoscarga+1]+(1-x_Ant[k])*F_Total_Ant[i+icasoscarga+1,k]/divisor*x_coords
                Ry[i+icasoscarga+1]=Ry[i+icasoscarga+1]+F_Total_Ant[i+icasoscarga+1,k]/divisor
                #print(MzVer)
        icombinacao=i+icasoscarga+1
        Desloc[i+icasoscarga+1,:,:],Reactions[i+icasoscarga+1,:,:],Desloc_aux[i+icasoscarga+1,:],F[i+icasoscarga+1,:,:],F_axial[i+icasoscarga+1,:]=Calc_Load_Case(i+icasoscarga+1+iStart-1,coords,nos_por_elemento,todos_nos,nEle)
        Case_Out_csv=Shaft_Wind(Case_Out_csv,NameCaso[i+icasoscarga+1],At,Af,Ac,Acsup,As,Indice_Cheios,cf0f,cf0c,cf0csup,cfs0,cfs[i+icasoscarga+1,:],c0_Vento_For,CsCd[i+icasoscarga+1],Pressao_troco_medio,FmW[i+icasoscarga+1,:],FTW[i+icasoscarga+1,:],L_lin_equiv,Ka,FmW_lin,FTW_lin[i+icasoscarga+1,:],nos_vento)

    

    ops.remove('loadPattern',int(icombinacao+iStart-1))
    ops.wipeAnalysis()
    if C0_calc=="Sim":
        c0_top=calculoc0(Pais,zmin,[h_torre+Alt],Tipoc0,Alt_col,Lu,Ld,Xtopo,Ac,A500,A1000)
    else:
        c0_top=np.ones(1)
    Pressure_top=Pressure_wind(np.array([h_torre+Alt]),vb,z0,zmin,Alt,Pais,rho_ar,c0_top)
    v_top=Velocidade_media(np.array([h_torre+Alt]),vb,z0,zmin,Pais,c0_top)
    #vb,_,_=Zone(Zona,Terrain,Pais)
    q0=1/2*rho_ar*vb**2
    Wind_Out_Csv,Ice_Out_Csv=Wind_Ice_Out(Pais,z0,vb,zmin,c0_top[0],q0,v_top[0],Pressure_top[0],Gelo,esp_gelo,rho_gelo,Calc_Wind_Auto)
    return Desloc,Reactions,Desloc_aux,F,F_axial,Ry,MzVer,F_nos,nos_grupo,z_grupo,divisor,x_Ant,iAnt,F_Total_Ant, FTW_dist,FTW,FTW_lin,FTW_lin_dist,CsCd,cf0f,cf0c,cf0csup,cfs0,cfs,NameCaso,ivento,Massa_nos,Massa_nos_lin,alpha_vector,icombinacao,v_no_maximo,Leg_Type,Massa_Troco,Case_Out_csv,Wind_Out_Csv,Ice_Out_Csv




def Simple_Cases_Gelo(Ncasoscarga,nos_por_classe,C0_calc,Pais,Terrain,Zona,zmax,Alt,Tipoc0,Alt_col,Lu,Ld,Xtopo,Ac_c0,A500,A1000,rho_ar,Dim,visco_ar,nos_vento,Troco,TrussType,z_vento,ElementosID,ExpVento,Comprimento_Barra,divisor,esp_gelo,Gelo,h_torre,munit,Massa_Extra,todos_nos,nEle,nos_por_elemento,coords,h_vetor,z_grupo,nos_grupo,dh,torre,dalpha,Posicao_Antenna,Massa_Antenna,AreaDirection,equipline,Ant,iStart,Esp,Area,rho_gelo,Massa_Antenna_Gelo,Massa_nos_lin,Massa_Troco,vb,z0,zmin):    

    A_lin_Troco=np.zeros(len(nos_por_classe)+1)
    A_lin_Troco_Cf=np.zeros(len(nos_por_classe)+1)
    L_lin_equiv=np.zeros(len(nos_por_classe)+1)
    Z_med_vento=np.zeros(len(nos_por_classe)+1)
    if C0_calc=="Sim":
        c0_no_maximo=calculoc0(Pais,zmin,zmax+Alt,Tipoc0,Alt_col,Lu,Ld,Xtopo,Ac_c0,A500,A1000)
    else:
        c0_no_maximo=np.ones(len(zmax+Alt))
    Pressao_no_maximo=Pressure_wind(zmax+Alt,vb,z0,zmin,Alt,Pais,rho_ar,c0_no_maximo)
    v_no_maximo=np.sqrt(Pressao_no_maximo*2/rho_ar)
    Reynolds=v_no_maximo*Dim/visco_ar

    nos_vento_troco,Largura,Largura_total,Largura_media,At,Af,Ac,Acsup,As,Indice_Cheios=Elements_Areas(nos_por_classe,nos_vento,Troco,TrussType,Dim,z_vento,ElementosID,ExpVento,Comprimento_Barra,divisor,Reynolds,esp_gelo)


    for i in range(len(z_vento)-1):
        Z_med_vento[i+1]=(z_vento[i]+z_vento[i+1])/2
    if C0_calc=="Sim":
        c0_troco_medio=calculoc0(Pais,zmin,Z_med_vento+Alt,Tipoc0,Alt_col,Lu,Ld,Xtopo,Ac_c0,A500,A1000)
    else:
        c0_troco_medio=np.ones(len(Z_med_vento+Alt))

    Pressao_troco_medio=Pressure_wind(Z_med_vento+Alt,vb,z0,zmin,Alt,Pais,rho_ar,c0_troco_medio)

    Mascara= (TrussType == "Leg")
    Contagem_Montante_Troco = np.bincount(Troco[Mascara].astype(int), minlength=int(np.max(Troco)+1)) # ignora índice 0
    if C0_calc=="Sim":
        c0_Zs=calculoc0(Pais,zmin,np.array([0.6*h_torre+Alt]),Tipoc0,Alt_col,Lu,Ld,Xtopo,Ac_c0,A500,A1000)
    else:
        c0_Zs=np.ones(1)
    vm=Velocidade_media(np.array([0.6*h_torre+Alt]),Zona,Terrain,Pais,c0_Zs)
    Iv_Zs=Iv_calc(np.array([0.6*h_torre])+Alt,z0,zmin,Alt,Pais,c0_Zs)
    if ExpVento[int(Mascara[0])]=="Circ":
        Leg_Type="Circ"
    else:
        Leg_Type="Flat"



    Massa_Troco_gelo=np.zeros(len(nos_por_classe)+1)
    
    Massa_Barra=np.zeros(nEle+1)
    iAngle=ExpVento=="Flat"
    #Massa_Barra[iAngle]= ((2*(Dim[iAngle]+2*esp_gelo)*(Esp[iAngle]+2*esp_gelo)-(Esp[iAngle]+2*esp_gelo)**2) - Area[iAngle])*Comprimento_Barra[iAngle]*rho_gelo
    Massa_Barra[iAngle]=((Dim[iAngle]+2*esp_gelo)*(Esp[iAngle]+2*esp_gelo)+(Esp[iAngle]+2*esp_gelo)*(Dim[iAngle]+2*esp_gelo-(Esp[iAngle]+2*esp_gelo))-Area[iAngle])*Comprimento_Barra[iAngle]*rho_gelo
    iPipe=ExpVento=="Circ"
    Massa_Barra[iPipe]= (np.pi/4*((Dim[iPipe]+2*esp_gelo)**2-(Dim[iPipe])**2))*Comprimento_Barra[iPipe]*rho_gelo
    for i in range(len(nos_por_classe)+1):
        Massa_Troco_gelo[i] = Massa_Barra[Troco == i].sum()
    #Massa_Troco=Massa_Troco*(1+Massa_Extra/100)
    Massa_Troco_gelo=Massa_Troco_gelo
    Zs=0.6*h_torre
    iZs=encontra_vizinhos(z_vento[::-1], Zs)
    Largura_ref=Largura[:,iZs[0]+1]
    Larg=np.zeros(3)
    Larg[0]=Largura[1,-1]
    Larg[2]=Largura[0,1]
    if Largura_ref[0]==Largura_ref[1]:
        Larg[1]=Largura_ref[0]
    else:
        h1=z_vento[iZs[0]]
        h2=z_vento[iZs[1]]
        L1=Largura_ref[0]
        L2=Largura_ref[1]
        Larg[1]=L2+(Zs-h2)*(L1-L2)/(h1-h2)


    cfs=np.zeros((Ncasoscarga,len(nos_por_classe)+1))
    FmW=np.zeros((Ncasoscarga,len(nos_por_classe)+1))
    FTW=np.zeros((Ncasoscarga,len(nos_por_classe)+1))
    FTW_lin=np.zeros((Ncasoscarga,len(nos_por_classe)+1))
    FTW_dist=np.zeros((Ncasoscarga,len(nos_por_classe)+1))
    FTW_lin_dist=np.zeros((Ncasoscarga,len(nos_por_classe)+1))
    F_Troco=np.zeros((Ncasoscarga,len(nos_por_classe)+1))
    F_nos=np.zeros((Ncasoscarga,len(todos_nos)+2,6))
    Reactions=np.zeros((Ncasoscarga,len(todos_nos)+2,6))
    F=np.zeros((Ncasoscarga,nEle+1,12))
    F_axial=np.zeros((Ncasoscarga,nEle+1))
    Desloc=np.zeros((Ncasoscarga,len(todos_nos)+2,6))
    MzVer=np.zeros(Ncasoscarga)
    Ry=np.zeros(Ncasoscarga)
    Desloc_aux=np.zeros((Ncasoscarga,len(todos_nos)+2))
    CsCd=np.zeros((Ncasoscarga))
    NameCaso=np.zeros((Ncasoscarga)).astype(str)
    #######Caso de carga 1 - Peso Proprio #################



    NameCaso[1]= Create_Load_Case(1+iStart-1,"Peso Proprio")

    for ele, nos in nos_por_elemento.items():
        # ops.load(nos[0],0,0,-Massa_Barra[ele]*9.81/2*(1+Massa_Extra/100),0,0,0)
        # ops.load(nos[1],0,0,-Massa_Barra[ele]*9.81/2*(1+Massa_Extra/100),0,0,0)
        # F_nos[1,int(nos[0]),2]=F_nos[1,int(nos[0]),2]-Massa_Barra[ele]*9.81/2*(1+Massa_Extra/100)
        # F_nos[1,int(nos[1]),2]=F_nos[1,int(nos[1]),2]-Massa_Barra[ele]*9.81/2*(1+Massa_Extra/100)
        ops.load(nos[0],0,0,-Massa_Barra[ele]*9.81/2,0,0,0)
        ops.load(nos[1],0,0,-Massa_Barra[ele]*9.81/2,0,0,0)
        F_nos[1,int(nos[0]),2]=F_nos[1,int(nos[0]),2]-Massa_Barra[ele]*9.81/2
        F_nos[1,int(nos[1]),2]=F_nos[1,int(nos[1]),2]-Massa_Barra[ele]*9.81/2

    Desloc[1,:,:],Reactions[1,:,:],Desloc_aux[1,:],F[1,:,:],F_axial[1,:]=Calc_Load_Case(1+iStart-1,coords,nos_por_elemento,todos_nos,nEle)
    # for tag, u in desloc.items():
    #     # arredonda cada componente a 5 casas
    #     u_rounded = [round(val, 5) for val in u]
    #     print(f"Desloc nó {tag}: {u_rounded}")

    # ops.reactions()
    # reaction_1 = ops.nodeReaction(1)
    # print(f"Reações no nó 1: Rx = {reaction_1[0]} N, Ry = {reaction_1[1]} N, Rz = {reaction_1[2]} N")
    # reaction_1 = ops.nodeReaction(25)
    # print(f"Reações no nó 25: Rx = {reaction_1[0]} N, Ry = {reaction_1[1]} N, Rz = {reaction_1[2]} N")
    # reaction_1 = ops.nodeReaction(37)
    # print(f"Reações no nó 37: Rx = {reaction_1[0]} N, Ry = {reaction_1[1]} N, Rz = {reaction_1[2]} N")


    #######Caso de carga Vento #################
    icasoscarga=1
    alpha_vector=np.arange(0, 360, dalpha)

    ivento=0 




    if C0_calc=="Sim":
        c0_Vento_Ant=calculoc0(Pais,zmin,Posicao_Antenna+Alt,Tipoc0,Alt_col,Lu,Ld,Xtopo,Ac_c0,A500,A1000)
    else:
        c0_Vento_Ant=np.ones(len(Posicao_Antenna+Alt))
    Pressao_Vento_Ant=Pressure_wind(Posicao_Antenna+Alt,vb,z0,zmin,Alt,Pais,rho_ar, c0_Vento_Ant)  ##### Tratar do c0
    F_Total_Ant=np.zeros((Ncasoscarga,AreaDirection.shape[0]))

    iAnt,x_Ant, Massa_nos_gelo,Massa_nos_lin_gelo=Area_Position_Mass(Posicao_Antenna,h_torre,z_grupo,nos_grupo,Massa_Antenna_Gelo,divisor,equipline,TrussType,todos_nos,Gelo)

    for i in todos_nos:  
        ops.mass(int(i),Massa_nos_gelo[int(i)]+Massa_nos_lin_gelo[int(i)]+np.abs(F_nos[1,int(i),2]/9.81),Massa_nos_gelo[int(i)]+Massa_nos_lin_gelo[int(i)]-(F_nos[1,int(i),2]/9.81),Massa_nos_gelo[int(i)]+Massa_nos_lin_gelo[int(i)]-(F_nos[1,int(i),2]/9.81) ,0,0,0)
    # === Modelo modal
    num_modes = 10
    eigenvalues = ops.eigen(num_modes)  # retorna ω²

    # Converter para frequências e períodos
    omegas = np.sqrt(np.array(eigenvalues))
    freqs = omegas / (2*np.pi)
    # print("\n=== Resultados Modais ===")
    # for i in range(num_modes):
    #     print(f"Modo {i+1}: ω²={eigenvalues[i]:.4e}, f={freqs[i]:.4f} Hz")

    #### Caso 2 - Peso equipamentos

    icasoscarga=icasoscarga+1
    NameCaso[icasoscarga]= Create_Load_Case(icasoscarga+iStart-1,"Peso Equipamentos - Linear")

    for i in todos_nos:
        if i!=0:
            ops.load(int(i), 0,0, -(Massa_nos_lin_gelo[i]-Massa_nos_lin[i])*9.81,0,0,0)
            F_nos[icasoscarga,int(i),2]=F_nos[icasoscarga,int(i),2]-(Massa_nos_lin_gelo[i]-Massa_nos_lin[i])*9.81
      
    Desloc[icasoscarga,:,:],Reactions[icasoscarga,:,:],Desloc_aux[icasoscarga,:],F[icasoscarga,:,:],F_axial[icasoscarga,:]=Calc_Load_Case(icasoscarga+iStart-1,coords,nos_por_elemento,todos_nos,nEle)

    #### Caso 3 - Peso equipamentos

    icasoscarga=icasoscarga+1
    NameCaso[icasoscarga]= Create_Load_Case(icasoscarga+iStart-1,"Peso Equipamento - Pontual")
    for i in range(len(Posicao_Antenna)):
        for j in range(nos_grupo.shape[1]):
            ops.load(int(nos_grupo[int(iAnt[i,1]),j]), 0,0, -(1-x_Ant[i])*(Massa_Antenna_Gelo[i]-Massa_Antenna[i])/divisor*9.81,0,0,0)
            ops.load(int(nos_grupo[int(iAnt[i,0]),j]), 0,0,-x_Ant[i]*(Massa_Antenna_Gelo[i]-Massa_Antenna[i])/divisor*9.81,0,0,0)
            F_nos[icasoscarga,int(nos_grupo[int(iAnt[i,1]),j]),2]=F_nos[icasoscarga,int(nos_grupo[int(iAnt[i,1]),j]),2]-(1-x_Ant[i])*(Massa_Antenna_Gelo[i]-Massa_Antenna[i])/divisor*9.81
            F_nos[icasoscarga,int(nos_grupo[int(iAnt[i,0]),j]),2]=F_nos[icasoscarga,int(nos_grupo[int(iAnt[i,0]),j]),2]-x_Ant[i]*(Massa_Antenna_Gelo[i]-Massa_Antenna[i])/divisor*9.81
    Desloc[icasoscarga,:,:],Reactions[icasoscarga,:,:],Desloc_aux[icasoscarga,:],F[icasoscarga,:,:],F_axial[icasoscarga,:]=Calc_Load_Case(icasoscarga+iStart-1,coords,nos_por_elemento,todos_nos,nEle)


    try:
        Ant_max=max(Ant[:,0])
    except:
        Ant_max=0
    if Ant_max>h_torre:
        h_equip=np.arange(0,max(Ant[:,0])+dh,dh)
    else:
        h_equip=h_vetor
    me,equip_mass_loc,equip_mass_dist,m_vector=Equivalent_Mass(h_vetor,Posicao_Antenna,Ant,h_torre,Massa_Antenna_Gelo,equipline,dh,z_vento,Massa_Troco_gelo+Massa_Troco,h_equip,Gelo)


    for i in range(len(z_vento)-1):
        for j in range(equipline.shape[0]):

            if z_vento[i+1]<=equipline[j,-1] and z_vento[i+1]>=equipline[j,-2]:
                if equipline[j,-1]>=z_vento[i]:
                    dL=min(z_vento[i]-z_vento[i+1],z_vento[i]-equipline[j,-2])
                else:
                    dL=equipline[j,-1]-z_vento[i+1]
                
                A_lin_Troco[i+1]=A_lin_Troco[i+1]+equipline[j,0]*dL
                A_lin_Troco_Cf[i+1]=A_lin_Troco_Cf[i+1]+equipline[j,0]*dL*equipline[j,1]
        L_lin_equiv[i+1]=A_lin_Troco[i+1]/(z_vento[i]-z_vento[i+1])



    Ka=CalcKa(A_lin_Troco,As,At,L_lin_equiv,Largura_media)
    z_vento_roll=np.roll(z_vento,1)
    FmW_lin=Pressao_troco_medio/(1+7*Iv_Zs)*A_lin_Troco_Cf*Ka
    if C0_calc=="Sim":
        c0_Vento_For=calculoc0(Pais,zmin,z_vento_roll+Alt,Tipoc0,Alt_col,Lu,Ld,Xtopo,Ac_c0,A500,A1000)
    else:
        c0_Vento_For=np.ones(len(z_vento_roll+Alt))
    ivento=icasoscarga
    Case_Out_csv=np.array([["Expanded Results Gelo","","","","","","","","","","","","","","","","","","","",""]])
    for i in range(len(alpha_vector)):
        
        ops.remove('loadPattern',int(i+icasoscarga))
        ops.wipeAnalysis()
        NameCaso[i+icasoscarga+1]= Create_Load_Case(i+icasoscarga+1+iStart-1,f"Wind Load - {alpha_vector[i]}")
        cf0f,cf0c,cf0csup,cfs0,cfs[i+icasoscarga+1,:]= Coef_forca(Af,Ac,Acsup,As,alpha_vector[i],Indice_Cheios,torre)
        CsCd[i+icasoscarga+1],CsCd_out=CsCd_calc(0.6*h_torre+Alt,h_torre,Iv_Zs,Larg,z0,zmin,freqs[0],me,cfs[i+icasoscarga+1,:],rho_ar,vm,Alt,Pais)
        FmW[i+icasoscarga+1,:]=Pressao_troco_medio/(1+7*Iv_Zs)*cfs[i+icasoscarga+1,:]*As
        FTW[i+icasoscarga+1,:]=FmW[i+icasoscarga+1,:]*(1+(1+0.2*(z_vento_roll/h_torre)**2)*((1+7*Iv_Zs)*CsCd[i+icasoscarga+1]-1)/c0_Vento_For)
        FTW_lin[i+icasoscarga+1,:]=FmW_lin*(1+(1+0.2*(z_vento_roll/h_torre)**2)*((1+7*Iv_Zs)*CsCd[i+icasoscarga+1]-1)/c0_Vento_For)
        for k in range(len(nos_vento)-1):
            Mascara = (Troco[1:] == int(k+1)) & (TrussType[1:] == "Leg")
            FTW_dist[i+icasoscarga+1,k+1]=FTW[i+icasoscarga+1,k+1]/divisor/(z_vento[k]-z_vento[k+1])
            FTW_lin_dist[i+icasoscarga+1,k+1]=FTW_lin[i+icasoscarga+1,k+1]/divisor/(z_vento[k]-z_vento[k+1])
            for ele in ElementosID[Mascara]:
                ops.eleLoad("-ele", int(ele), "-type", "beamUniform",-(FTW_dist[i+icasoscarga+1,k+1]+FTW_lin_dist[i+icasoscarga+1,k+1])*np.sin(np.radians(alpha_vector[i])),(FTW_dist[i+icasoscarga+1,k+1]+FTW_lin_dist[i+icasoscarga+1,k+1])*np.cos(np.radians(alpha_vector[i])),0)
                # nós do elemento
                nodes = ops.eleNodes(int(ele))
                # coordenada x de cada nó
                x_coords = [ops.nodeCoord(n)[2] for n in nodes]
                Ry[i+icasoscarga+1]=Ry[i+icasoscarga+1]+(FTW_dist[i+icasoscarga+1,k+1]+FTW_lin_dist[i+icasoscarga+1,k+1])*Comprimento_Barra[int(ele)]
                MzVer[i+icasoscarga+1]=MzVer[i+icasoscarga+1]+(FTW_dist[i+icasoscarga+1,k+1]+FTW_lin_dist[i+icasoscarga+1,k+1])*Comprimento_Barra[int(ele)]*np.mean(x_coords)
                #print(MzVer)
        F_Total_Ant[i+icasoscarga+1,:] = (Pressao_Vento_Ant * AreaDirection[:,i]) * CsCd[i+icasoscarga+1]
        
        for k in range(len(Posicao_Antenna)):
            for j in range(nos_grupo.shape[1]):
                ops.load(int(nos_grupo[int(iAnt[k,1]),j]), -(1-x_Ant[k])*F_Total_Ant[i+icasoscarga+1,k]*np.sin(np.radians(alpha_vector[i]))/divisor, (1-x_Ant[k])*F_Total_Ant[i+icasoscarga+1,k]*np.cos(np.radians(alpha_vector[i]))/divisor,0,0,0,0)
                ops.load(int(nos_grupo[int(iAnt[k,0]),j]), -x_Ant[k]*F_Total_Ant[i+icasoscarga+1,k]*np.sin(np.radians(alpha_vector[i]))/divisor , x_Ant[k]*F_Total_Ant[i+icasoscarga+1,k]*np.cos(np.radians(alpha_vector[i]))/divisor,0,0,0,0)
                F_nos[i+icasoscarga+1,int(nos_grupo[int(iAnt[k,1]),j]),0]=F_nos[i+icasoscarga+1,int(nos_grupo[int(iAnt[k,1]),j]),0] -(1-x_Ant[k])*F_Total_Ant[i+icasoscarga+1,k]*np.sin(np.radians(alpha_vector[i]))/divisor
                F_nos[i+icasoscarga+1,int(nos_grupo[int(iAnt[k,1]),j]),1]=F_nos[i+icasoscarga+1,int(nos_grupo[int(iAnt[k,1]),j]),1]+(1-x_Ant[k])*F_Total_Ant[i+icasoscarga+1,k]*np.cos(np.radians(alpha_vector[i]))/divisor
                F_nos[i+icasoscarga+1,int(nos_grupo[int(iAnt[k,0]),j]),0]=F_nos[i+icasoscarga+1,int(nos_grupo[int(iAnt[k,0]),j]),0]-x_Ant[k]*F_Total_Ant[i+icasoscarga+1,k]*np.sin(np.radians(alpha_vector[i]))/divisor
                F_nos[i+icasoscarga+1,int(nos_grupo[int(iAnt[k,0]),j]),1]=F_nos[i+icasoscarga+1,int(nos_grupo[int(iAnt[k,0]),j]),1]+x_Ant[k]*F_Total_Ant[i+icasoscarga+1,k]*np.cos(np.radians(alpha_vector[i]))/divisor
                x_coords = ops.nodeCoord(int(nos_grupo[int(iAnt[k,0]),j]))[2]
                MzVer[i+icasoscarga+1]=MzVer[i+icasoscarga+1]+x_Ant[k]*F_Total_Ant[i+icasoscarga+1,k]/divisor*x_coords
                x_coords = ops.nodeCoord(int(nos_grupo[int(iAnt[k,1]),j]))[2]
                MzVer[i+icasoscarga+1]=MzVer[i+icasoscarga+1]+(1-x_Ant[k])*F_Total_Ant[i+icasoscarga+1,k]/divisor*x_coords
                Ry[i+icasoscarga+1]=Ry[i+icasoscarga+1]+F_Total_Ant[i+icasoscarga+1,k]/divisor
                #print(MzVer)
        icombinacao=i+icasoscarga+1
        Desloc[i+icasoscarga+1,:,:],Reactions[i+icasoscarga+1,:,:],Desloc_aux[i+icasoscarga+1,:],F[i+icasoscarga+1,:,:],F_axial[i+icasoscarga+1,:]=Calc_Load_Case(i+icasoscarga+1+iStart-1,coords,nos_por_elemento,todos_nos,nEle)
        Case_Out_csv=Shaft_Wind(Case_Out_csv,NameCaso[i+icasoscarga+1],At,Af,Ac,Acsup,As,Indice_Cheios,cf0f,cf0c,cf0csup,cfs0,cfs[i+icasoscarga+1,:],c0_Vento_For,CsCd[i+icasoscarga+1],Pressao_troco_medio,FmW[i+icasoscarga+1,:],FTW[i+icasoscarga+1,:],L_lin_equiv,Ka,FmW_lin,FTW_lin[i+icasoscarga+1,:],nos_vento)
        _,Ice_Out_Csv=Wind_Ice_Out(0,0,0,0,0,0,0,0,Gelo,esp_gelo,rho_gelo)
    ops.remove('loadPattern',int(icombinacao+iStart-1))
    ops.wipeAnalysis()

    
    return Desloc,Reactions,Desloc_aux,F,F_axial,Ry,MzVer,F_nos,nos_grupo,z_grupo,divisor,x_Ant,iAnt,F_Total_Ant, FTW_dist,FTW,FTW_lin,FTW_lin_dist,CsCd,cf0f,cf0c,cf0csup,cfs0,cfs,NameCaso,ivento,Massa_nos_gelo,Massa_nos_lin_gelo,alpha_vector,icombinacao,v_no_maximo,Leg_Type,Case_Out_csv,Ice_Out_Csv