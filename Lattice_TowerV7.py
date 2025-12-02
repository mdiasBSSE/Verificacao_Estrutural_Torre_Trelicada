import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # necessário para plot 3D
from openseespy import opensees as ops
import os
from pathlib import Path
import sys
from Utilities.Utilities_EC import ClasseFiabilidade
from Equipments.Calc_Areas import AntenaOrg
from Buckling.Buckling_ops import Buckling_Function
from Connections.Bolted_Connection import Bolted_Connection_Bracing,Bolted_Connection_Shaft,Flanged_Bolted_Connection,Flanged_Foundation_Angle
from Inputs.Lattice_Tower_Input import Model_csv_read,Structure_txt_read
from Model_Generation.Node_Filters import Node_Wind,Node_Class_Troco
from Model_Generation.Model_Generation import Model_Generation
from Model_Calculation.Load_Creation_Area import Area_Definition
from Model_Calculation.Simple_Cases import Simple_Cases, Simple_Cases_Gelo
from Model_Calculation.Combined_Cases import ELS,ELU,Comb_Gelo
from Output.Output_Lattice_Reduced import Desloc_Output,Basal_Effort_Output,Lattice_csv_Output
from Output.Output_Lattice_Expanded import Output_Areas,salvar_matrizes_csv
from Fondation.Foundation import Foundation_Calc_Basic
from Wind.Wind_Ec import Zone
import time
import csv
from scipy.optimize import minimize
def main(Aincrease):    
    # === 1. Ler nós e barras ===


    # Verifica se foi passado algum argumento
    if len(sys.argv) > 1:
        # Pega o primeiro argumento passado na linha de comando
        argumento = sys.argv[1]

    if getattr(sys, 'frozen', False):
        FolderPath=Path(sys.executable).resolve().parent
    else:
        FolderPath=Path(__file__).resolve().parent
    os.chdir(FolderPath)
    # Constrói o caminho completo para o arquivo .txt
    txt_file = os.path.join(FolderPath, "wind_config", f"{argumento}_wind_config.txt")
    ant_txt_file=os.path.join(FolderPath, "antennas", f"{argumento}_antennas.txt")
    Model_csv_file=os.path.join(FolderPath, "latticesurvey", f"{argumento}_latticesurvey.csv")

    # === 2. Criar modelo OpenSees ===
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)  # 3D, 3 dof/translação por nó

    Nodes_Matrix  ,Elements_Matrix  ,Dig_Hor_Matrix ,Shaft_Flange_Matrix,Shaft_Bolted_Matrix,Flange_Foundation_Angle,Foundation_Matrix,alpha_norte = Model_csv_read(Model_csv_file)
    (Type_Profile, Size_Profile, Dim, Esp, Area, Inertiay, Inertiaz,
        Inertiav, Inertiau, munit, ExpVento, E, Sigma_u, Sigma,
        Conection_Ele, N_Bolt_Ele, Shape, Top_bars, Middle_bars,
        zmax, zmed, Comprimento_Barra, ElementosID, nos_grupo,
        z_grupo, nEle, Node_Type, torre, TrussType,nos_montante,Rigidez,coords) = Model_Generation(Nodes_Matrix, Elements_Matrix)
    Regime=np.zeros(nEle+1).astype(str)
    Travamento=np.zeros(nEle+1).astype(str)


    (Pais,vb,z0,zmin,Classe_Fiabilidade,
        Gelo,Altitude,C0_calc,Tipoc0,
        Alt_col,Lu,Ld,Xtopo,Ac_c0,A500,A1000,Alt,
        FundMethod,MfundMax,RatioFundCalc,
        SoilSelfWeight,AllowedSoilTension_SLS,
        AllowedSoilTension_ULS)=Structure_txt_read(txt_file)
    zmed_alt=Alt+zmed
    if Pais=="France":
        rho_ar=1.225
    else:
        rho_ar=1.25
    c0=np.ones(len(zmed))
    #alpha_norte=0

    GammaM1=1.0
    GammaM0=1.0
    GammaM2=1.25
    visco_ar=15*10**-6
    h_torre=np.max(Nodes_Matrix[:,4].astype(float))*10**-3
    Massa_Extra=12 ####Peso extra em percentagem (Chapa, galva, parafusos, etc)
    Massa_Extra_gelo=0
    Discretização=200
    h_vetor = np.linspace(0, h_torre, Discretização)
    dh=h_vetor[1]-h_vetor[0]
    m_vector=np.zeros(len(h_vetor))

    if Gelo=="Sim": #Ver neve para inglaterra
        if Pais in ["Portugal","Portugal-RSA","France"]:
            if Altitude<600:
                esp_gelo=0.02
                Classe_gelo="G2"
                k_gelo=0.45
            elif Altitude<1200:
                esp_gelo=0.04
                Classe_gelo="G4"
                k_gelo=0.55
            else:
                esp_gelo=0.06
            if Pais=="France":
                k_gelo=0.64
                rho_gelo=600
            else:
                rho_gelo=900
            if Pais=="France":
                if h_torre<100:
                    kh=1.2
                elif h_torre<200:
                    kh=1.4
                else:
                    kh=1.7
            else:
                kh=1
        elif Pais=="United Kingdom":
            rho_gelo=9000/9.81 ###Pag12 NA.UK
        elif Pais=="Germany":
            rho_gelo=7000/9.81
    else:
        esp_gelo=0
        kh=0
        rho_gelo=0
    Psiventogelo=0.6 #NF NA pag11




    nos_vento,z_vento,divisor=Node_Wind(Nodes_Matrix,Node_Type,torre)


    nos_por_classe,Troco,todos_nos,nos_por_elemento=Node_Class_Troco(z_vento,nEle)
    Ncasoscarga=300
    nos_vento = nos_vento[::-1, :]
    z_vento=z_vento[::-1]
    dalpha=15
    (Ant, Nant, AreaAnt,AreaDirection,AreaDirectiongelo,equipline,Line_Equipment,
    Auto,Ladder,CabosOutside,Final,Nanalise,Org,ValueLine,NameLine,
    Massa_Antenna,Massa_Antenna_Gelo,Posicao_Antenna,Antenas)=Area_Definition(ant_txt_file,h_vetor,h_torre,Gelo,esp_gelo,rho_gelo,dalpha,alpha_norte,Aincrease)


    (Desloc,Reactions,Desloc_aux,F,F_axial,Ry,MzVer,F_nos,nos_grupo,
    z_grupo,divisor,x_Ant,iAnt,F_Total_Ant, FTW_dist,FTW,FTW_lin,
    FTW_lin_dist,CsCd,cf0f,cf0c,cf0csup,cfs0,cfs,NameCaso,ivento,
    Massa_nos,Massa_nos_lin,alpha_vector,icombinacao,v_no_maximo,Leg_Type,Massa_Troco,Case_Out_csv,Wind_Out_Csv,Ice_Out_Csv)=Simple_Cases(Ncasoscarga,nos_por_classe,C0_calc,Pais,zmax,
                                                                                                    Alt,Tipoc0,Alt_col,Lu,Ld,Xtopo,Ac_c0,A500,A1000,rho_ar,Dim,
                                                                                                     visco_ar,nos_vento,Troco,TrussType,z_vento,ElementosID,
                                                                                                    ExpVento,Comprimento_Barra,divisor,0,"Nao",h_torre,munit,
                                                                                                       Massa_Extra,todos_nos,nEle,nos_por_elemento,coords,h_vetor,
                                                                                                    z_grupo,nos_grupo,dh,torre,dalpha,Posicao_Antenna,
                                                                                                    Massa_Antenna,AreaDirection,equipline,Ant,1,rho_gelo,vb,z0,zmin)
        

    N_basal=np.abs(np.sum(F_nos[:4,:,2]))
    #for i in range(ops.getNodes()):

    if Gelo=="Sim":
        Ncasoscarga_gelo=200
        (Desloc_gelo,Reactions_gelo,Desloc_aux_gelo,F_gelo,F_axial_gelo,Ry_gelo,MzVer_gelo,F_nos_gelo,nos_grupo_gelo,
        z_grupo,divisor,x_Ant,iAnt,F_Total_Ant_gelo, FTW_dist_gelo,FTW_gelo,FTW_lin_gelo,
        FTW_lin_dist_gelo,CsCd_gelo,cf0f_gelo,cf0c_gelo,cf0csup_gelo,cfs0_gelo,cfs_gelo,NameCaso_gelo,icasoscarga_gelo,
        Massa_nos_gelo,Massa_nos_lin_gelo,alpha_vector,icomb_gelo,v_no_maximo,Leg_Type,Case_Gelo_Out_csv,Ice_Out_Csv)=Simple_Cases_Gelo(Ncasoscarga,nos_por_classe,C0_calc,Pais,zmax,Alt,
                                                                                                        Tipoc0,Alt_col,Lu,Ld,Xtopo,Ac_c0,A500,A1000,rho_ar,Dim,visco_ar,
                                                                                                        nos_vento,Troco,TrussType,z_vento,ElementosID,ExpVento,
                                                                                                        Comprimento_Barra,divisor,esp_gelo,Gelo,h_torre,munit,0,
                                                                                                        todos_nos,nEle,nos_por_elemento,coords,h_vetor,z_grupo,nos_grupo,dh,
                                                                                                        torre,dalpha,Posicao_Antenna,Massa_Antenna,AreaDirectiongelo,equipline
                                                                                                        ,Ant,icombinacao+1,Esp,Area,rho_gelo,Massa_Antenna_Gelo,Massa_nos_lin,Massa_Troco,vb,z0,zmin)
        N_basal_gelo=np.abs(np.sum(F_nos_gelo[:4,:,2]))

    Gamma_Perm_Des,Gamma_Variable_Des,Gamma_Perm_Fav,Gamma_Variable_Fav = ClasseFiabilidade(Pais,Classe_Fiabilidade)

    iSLS=icombinacao+1

    for i in range(len(alpha_vector)):
        ###ELS###
        NameCaso[icombinacao+i+1]=f"ELS-Vento {alpha_vector[i]}"
        Desloc,F,F_axial,Reactions,Ry,MzVer,N_basal_ELS=ELS(Desloc,Ry,MzVer,N_basal,F,Reactions,F_axial,ivento+i+1,icombinacao+i+1)
    iholder=icombinacao+i+1
    iELU=iholder+1
    for i in range(len(alpha_vector)):
        ###ELU###
        Desloc,F,F_axial,Reactions,Ry,MzVer,N_basal_ELU_Des=ELU(Desloc,Ry,MzVer,N_basal,F,Reactions,F_axial,ivento+i+1,iholder+i+1,Gamma_Perm_Des,Gamma_Variable_Des)
        NameCaso[iholder+i+1]=f"ELU-Vento {alpha_vector[i]} Pdes"
    iholder=iholder+i+1
    for i in range(len(alpha_vector)):
        ###ELU###
        NameCaso[iholder+i+1]=f"ELU-Vento {alpha_vector[i]} Pfav"
        Desloc,F,F_axial,Reactions,Ry,MzVer,N_basal_ELU_Fav=ELU(Desloc,Ry,MzVer,N_basal,F,Reactions,F_axial,ivento+i+1,iholder+i+1,Gamma_Perm_Fav,Gamma_Variable_Des)
    iEnd=iholder+i+1
    if Gelo=="Sim":
        for i in range(len(alpha_vector)):
            #### ELS - Vento Dominante
            NameCaso_gelo[icomb_gelo+i+1]=f"ELS - Vento Dominante {alpha_vector[i]}"

            Desloc_gelo,F_gelo,F_axial_gelo,Reactions_gelo,Ry_gelo,MzVer_gelo,N_basal_ELS_gelo_1=Comb_Gelo(Desloc,Ry,MzVer,N_basal,F,F_nos,Reactions,F_axial,Desloc_gelo,Ry_gelo,MzVer_gelo,
                                                                                                           N_basal_gelo,F_gelo,Reactions_gelo,F_axial_gelo,icasoscarga_gelo+i+1,
                                                                                                           icomb_gelo+i+1,1.0,k_gelo,Psiventogelo)
        iholder=icomb_gelo+i+1
        for i in range(len(alpha_vector)):
            #### ELS - Gelo Dominante
            NameCaso_gelo[iholder+i+1]=f"ELS - Gelo Dominante {alpha_vector[i]}"
            Desloc_gelo,F_gelo,F_axial_gelo,Reactions_gelo,Ry_gelo,MzVer_gelo,N_basal_ELS_gelo_2=Comb_Gelo(Desloc,Ry,MzVer,N_basal,F,F_nos,Reactions,F_axial,Desloc_gelo,Ry_gelo,
                                                                                                           MzVer_gelo,N_basal_gelo,F_gelo,Reactions_gelo,
                                                                                                           F_axial_gelo,icasoscarga_gelo+i+1,iholder+i+1,1.0,
                                                                                                           k_gelo*Psiventogelo,1.0)
        iholder=iholder+i+1
        iELU_gelo=iholder+1
        for i in range(len(alpha_vector)):
            ### ELU-Vento Dominante ###
            NameCaso_gelo[iholder+i+1]=f"ELU-Vento Dominante {alpha_vector[i]} Pdes"
            Desloc_gelo,F_gelo,F_axial_gelo,Reactions_gelo,Ry_gelo,MzVer_gelo,N_basal_ELU_gelo_1=Comb_Gelo(Desloc,Ry,MzVer,N_basal,F,F_nos,Reactions,F_axial,Desloc_gelo,Ry_gelo,MzVer_gelo,N_basal_gelo,F_gelo,Reactions_gelo,F_axial_gelo,icasoscarga_gelo+i+1,iholder+i+1,
                                                                                                        Gamma_Perm_Des,Gamma_Variable_Des*k_gelo,Gamma_Variable_Des*Psiventogelo)
        iholder=iholder+i+1
        for i in range(len(alpha_vector)):
            ### ELU - Gelo Dominante ###
            NameCaso_gelo[iholder+i+1]=f"ELU-Gelo Dominante {alpha_vector[i]} Pdes"
            Desloc_gelo,F_gelo,F_axial_gelo,Reactions_gelo,Ry_gelo,MzVer_gelo,N_basal_ELU_gelo_2=Comb_Gelo(Desloc,Ry,MzVer,N_basal,F,F_nos,Reactions,F_axial,Desloc_gelo,Ry_gelo,MzVer_gelo,
                                                                                                           N_basal_gelo,F_gelo,Reactions_gelo,F_axial_gelo,icasoscarga_gelo+i+1,iholder+i+1,
                                                                                                        Gamma_Perm_Des,Gamma_Variable_Des*Psiventogelo*k_gelo,Gamma_Variable_Des)
        iholder=iholder+i+1
        for i in range(len(alpha_vector)):
            ### ELU - Vento Dominante ###
            NameCaso_gelo[iholder+i+1]=f"ELU-Vento Dominante {alpha_vector[i]} PFav"
            Desloc_gelo,F_gelo,F_axial_gelo,Reactions_gelo,Ry_gelo,MzVer_gelo,N_basal_ELU_gelo_3=Comb_Gelo(Desloc,Ry,MzVer,N_basal,F,F_nos,Reactions,F_axial,Desloc_gelo,Ry_gelo,MzVer_gelo,
                                                                                                           N_basal_gelo,F_gelo,Reactions_gelo,F_axial_gelo,icasoscarga_gelo+i+1,iholder+i+1,
                                                                                                            Gamma_Perm_Fav,Gamma_Variable_Des*k_gelo,Gamma_Variable_Des*Psiventogelo)
        iholder=iholder+i+1
        for i in range(len(alpha_vector)):
            ### ELU - Gelo Dominante ###
            NameCaso_gelo[iholder+i+1]=f"ELU-Gelo Dominante {alpha_vector[i]} PFav"
            Desloc_gelo,F_gelo,F_axial_gelo,Reactions_gelo,Ry_gelo,MzVer_gelo,N_basal_ELU_gelo_4=Comb_Gelo(Desloc,Ry,MzVer,N_basal,F,F_nos,Reactions,F_axial,Desloc_gelo,Ry_gelo,MzVer_gelo,
                                                                                                           N_basal_gelo,F_gelo,Reactions_gelo,F_axial_gelo,icasoscarga_gelo+i+1,iholder+i+1,
                                                                                                        Gamma_Perm_Fav,Gamma_Variable_Des*Psiventogelo*k_gelo,Gamma_Variable_Des)
        iholder=iholder+i+1


    if Gelo=="Sim":
        Basal_Ou_csv,Mfund=Basal_Effort_Output(Ry,Ry_gelo,MzVer,MzVer_gelo,Gelo,iELU,iELU_gelo,iSLS,N_basal_ELS,N_basal_ELU_Des, N_basal_ELU_Fav,
                                         N_basal_ELS_gelo_1, N_basal_ELS_gelo_2,N_basal_ELU_gelo_1, N_basal_ELU_gelo_2,N_basal_ELU_gelo_3, N_basal_ELU_gelo_4)
    else:
        Basal_Ou_csv,Mfund=Basal_Effort_Output(Ry,0,MzVer,0,Gelo,iELU,0,iSLS,N_basal_ELS,N_basal_ELU_Des, N_basal_ELU_Fav,0, 0,0, 0,0,0)



    if Gelo=="Sim":
        iAngle_Crit,Rot_SLS_top,Desloc_Out=Desloc_Output(Desloc,nos_grupo,iSLS,iELU,iELU_gelo,iEnd,ivento,z_grupo,h_torre,alpha_vector,v_no_maximo,
                                             Ncasoscarga,todos_nos,Gamma_Perm_Des,Gamma_Variable_Des,Gamma_Perm_Fav,Gelo,Desloc_gelo,
                                             Ncasoscarga_gelo,icasoscarga_gelo,k_gelo,Psiventogelo)
    else:
        iAngle_Crit,Rot_SLS_top,Desloc_Out=Desloc_Output(Desloc,nos_grupo,iSLS,iELU,0,iEnd,ivento,z_grupo,h_torre,alpha_vector,v_no_maximo,Ncasoscarga,
                                             todos_nos,Gamma_Perm_Des,Gamma_Variable_Des,Gamma_Perm_Fav,Gelo,0,0,0,0,0)

    while iAngle_Crit>len(alpha_vector):
        iAngle_Crit=iAngle_Crit-len(alpha_vector)
    if Final==1:
        AntOrg,OrgFinal,Mfinal=AntenaOrg(Antenas,np.deg2rad(alpha_vector[iAngle_Crit]),np.column_stack((np.zeros(len(AreaDirection[:,iAngle_Crit])), Massa_Antenna,AreaDirection[:,iAngle_Crit])),Org,alpha_norte)


    F_Tracao=np.max(F_axial,axis=0)
    F_Compressao=np.min(F_axial,axis=0)
    if Gelo=="Sim":
        Aux_F_axial=np.vstack((F_axial, F_axial_gelo))
        F_Tracao=np.max(Aux_F_axial,axis=0)
        F_Compressao=np.min(Aux_F_axial,axis=0)
    else:
        F_Tracao=np.max(F_axial,axis=0)
        F_Compressao=np.min(F_axial,axis=0)


    NTracRd=Sigma*Area
    Ratio_Enc,Enc_Out_csv,Enc_Out_Exp=Buckling_Function(nEle,Inertiav,Inertiau,Inertiay,Inertiaz,Area,Comprimento_Barra,Sigma,Rigidez,nos_vento,Shape,
                                            Travamento,Top_bars,Middle_bars,ExpVento,Dim,Esp,TrussType,nos_montante,Conection_Ele,N_Bolt_Ele,
                                            GammaM1,F_Compressao,Size_Profile,Elements_Matrix,NTracRd,Troco)





    #Lig_Out = np.column_stack((TrussOut,DimOut,EspOut,e1_gusset,p1_gusset,e2_gusset,p2_gusset,t_gusset,Bolt,N_Bolt_Lig,d_Bolt, d_0, As_bolt, A_bolt,alpha_v_Bolt,FvRd,Ratio_bolt_lig,Check_Bolt,fu_Lig,fy_Lig,alpha_d_end,alpha_d_inner,k1_edge,alpha_b_end,alpha_b_inner,FbRd_Edge_Inner,FbRd_Edge_End,FbRd,Ratio_Esmag_Lig,Check_Esmag,Anv_block,Ant_block,Veff1Rd,Veff2Rd,beta_Cant,Anet,NuRd,Ratio_Block_Lig,Check_Block))

    Ratio_Lig,Lig_Out_csv,Lig_Out_Exp=Bolted_Connection_Bracing(Dig_Hor_Matrix,F_Tracao,F_Compressao,TrussType,Dim,GammaM2,ExpVento,Esp,Area,Leg_Type,Sigma_u,Sigma,GammaM0,Troco,nos_vento)

    Ratio_Lig_B,B_Out_csv,B_Out_Exp=Bolted_Connection_Shaft(Shaft_Bolted_Matrix,TrussType,F_Tracao,F_Compressao,GammaM2,Sigma_u,Sigma,GammaM0,Area,Esp,Troco,nos_vento,divisor)


    #B_Out = np.column_stack((e1_B,p1_B,e2_B,p2_B,t_B,Bolt_B,N_Bolt_B,d_Bolt_B, d_0_B, As_bolt_B, A_bolt_B,alpha_v_Bolt_B,FvRd_B,Ratio_bolt_B,Check_Bolt_B,fu_B,fy_B,alpha_d_end_B,alpha_d_inner_B,k1_edge_B,alpha_b_end_B,alpha_b_inner_B,FbRd_Edge_Inner_B,FbRd_Edge_End_B,FbRd_B,Ratio_Esmag_B,Check_Esmag_B,Anv_block_B,Ant_block_B,Veff1Rd_B,Veff2Rd_B,beta_Cant_B,Anet_B,NuRd_B,Ratio_Block_B,Check_Block_B))


    #F_Out = np.column_stack((Ext_Per_F,D_between_Bolt_F,t_F,Steel_F,Bolt_F,N_Bolt_F,Bolt_Steel_F,Stiffeners_F,N_Stiff_F,t_Stiff_F,L_Stiff_F,h_Stiff_F,D_Ext_F,D_Bolt_Flange_F,D_int_F,t_Tube_F,aws_F, e1_F,e2_F,alpha_r0_F,d_Bolt_F,Lb_F,leff,alpha_r_F,epsilon_F,fub_F,fyb_F,n_F,m_F,elinha_F,FtRd,Rb_F,R_F,Re_F,k3_F,NT1Rd,NT2Rd,NT3Rd,NT4Rd,BpRd,NTRd,Ratio_F,Check_F))

    Ratio_F,F_Out_csv,F_Out_Exp=Flanged_Bolted_Connection(Shaft_Flange_Matrix,TrussType,F_Tracao,Dim,Esp,Sigma_u,Sigma,GammaM2,GammaM0,Area,nos_vento,divisor)

    RatioFund,Fund_Out_Geo_Csv,Fund_Out_Ratio_Csv=Foundation_Calc_Basic(Foundation_Matrix,FundMethod,MfundMax,RatioFundCalc,SoilSelfWeight,AllowedSoilTension_SLS,AllowedSoilTension_ULS,Mfund)

    Ratio_Flange_Fund, Flange_Fund_out=Flanged_Foundation_Angle(Flange_Foundation_Angle,TrussType,F_Tracao)






    timestamp = int(time.time())
    Ratio_max=np.round(np.nanmax([np.nanmax(Ratio_Enc),np.nanmax(Ratio_Lig),np.nanmax(Ratio_Lig_B),np.nanmax(Ratio_F),RatioFund]),2)
    Ratio_max_out=np.round(np.nanmax([np.nanmax(Ratio_Enc),np.nanmax(Ratio_Lig),np.nanmax(Ratio_Lig_B),np.nanmax(Ratio_F),np.nanmax(Ratio_Flange_Fund)]),2)
    Ratio_max_optimization=np.nanmax([np.nanmax(Ratio_Enc),np.nanmax(Ratio_Lig),np.nanmax(Ratio_Lig_B),np.nanmax(Ratio_F),RatioFund,np.nanmax(Ratio_Flange_Fund)])
    if Final==1:
        if Aincrease==0:
            timestamp = int(time.time())
            OutputEstrutural=os.path.join(FolderPath, f"resultados/{Nanalise}_{timestamp}_result_collocation.csv")
            Lattice_csv_Output(OutputEstrutural,Enc_Out_csv,Lig_Out_csv,B_Out_csv,F_Out_csv,Desloc_Out,Basal_Ou_csv,OrgFinal,NameLine,ValueLine,Fund_Out_Geo_Csv,Fund_Out_Ratio_Csv,Wind_Out_Csv,Ice_Out_Csv,Classe_Fiabilidade,Flange_Fund_out)
            OutputEstruturalExp=os.path.join(FolderPath, f"resultados/{Nanalise}_{timestamp}_Expanded_Results.csv")
            Area_Out_Exp=Output_Areas(AreaDirection,NameCaso[4:icombinacao+1])
            if Gelo=="Sim":
                Area_Gelo_Out_Exp=Output_Areas(AreaDirectiongelo,NameCaso_gelo[4:icomb_gelo+1])
            else:
                Area_Gelo_Out_Exp=[]
                Case_Gelo_Out_csv=[]
            salvar_matrizes_csv(OutputEstruturalExp,Case_Out_csv,Area_Out_Exp,Case_Gelo_Out_csv,Area_Gelo_Out_Exp,Enc_Out_Exp,Lig_Out_Exp,B_Out_Exp,F_Out_Exp)
            if FundMethod==True:
                print(f"{Ratio_max_out*100};{np.round(Rot_SLS_top,3)};{RatioFund*100};resultados/{Nanalise}_{timestamp}_result_collocation.csv")
            else:
                print(f"{Ratio_max_out*100};{np.round(Rot_SLS_top,3)};0;resultados/{Nanalise}_{timestamp}_result_collocation.csv")
        else:
            OutputEstrutural="0"
    else:
        if FundMethod==True:
            print(f"{Ratio_max_out*100};{np.round(Rot_SLS_top,3)};{RatioFund*100};resultados/xx_{timestamp}_result_collocation.csv")
        else:
            print(f"{Ratio_max_out*100};{np.round(Rot_SLS_top,3)};0;resultados/xx_{timestamp}_result_collocation.csv")   
        OutputEstrutural="0"
    #print("Limites Troço 1: \n",nos_vento_troco[:,:,1])
    #print("Limites Troço 2: \n",nos_vento_troco[:,:,2])
    #print("Limites Troço 3: \n",nos_vento_troco[:,:,3])
    #print("Limites Troço 4: \n",nos_vento_troco[:,:,4])
    #print("Limites Troço 5: \n",nos_vento_troco[:,:,5])
    #print("Limites Troço 6: \n",nos_vento_troco[:,:,6])

    #print("Barras Limite Troço 1: \n",Barras_vento_troco[:,:,1])
    #print("Barras Limite Troço 2: \n",Barras_vento_troco[:,:,2])
    #print("Barras Limite Troço 3: \n",Barras_vento_troco[:,:,3])
    #print("Barras Limite Troço 4: \n",Barras_vento_troco[:,:,4])
    #print("Barras Limite Troço 5: \n",Barras_vento_troco[:,:,5])
    #print("Barras Limite Troço 6: \n",Barras_vento_troco[:,:,6])

    # print("Largura: ",Largura)
    # print("Largura total: ", Largura_total)
    # print("Largura_media: ",Largura_media)
    # print("Area Bruta:", At)
    # print("Area Cantoneiras:", Af)
    # print("Area Subcritica:", Ac)
    # print("Area Supercritica:", Acsup)
    # print("Area bruta: ",At)
    # print("Indice Cheios: ",Indice_Cheios)

    # print("cf0f: ", cf0f)
    # print("cf0c: ",cf0c)
    # print("cf0csup: ",cf0csup)
    # print("cfs0: ",cfs0)
    # print("cfs: ",cfs)
    # print("Massa Barra", Massa_Barra)
    # print("Massa Troco", Massa_Troco)
    # print("Massa Vector", m_vector)
    # print("Z_vento",z_vento)
    # print("me",me)
    # print("CsCd", CsCd)
    #print("Z_med_vento: ",Z_med_vento )
    #print("Pressao_troco_medio: ",Pressao_troco_medio )
    #print("Contagem_Montante_Troco: ",Contagem_Montante_Troco)
    # print("FmW", FmW)
    # print("FTW", FTW)

    # with  pd.ExcelWriter("Desloc.xlsx", engine="openpyxl") as writer:
    #     for i in range(Desloc.shape[0]):
    #         df = pd.DataFrame(Desloc[i, 1:, :])  # fatia i -> 2D
    #         df.to_excel(writer, sheet_name=NameCaso[i], index=False, header=False)
    # with pd.ExcelWriter("F.xlsx", engine="openpyxl") as writer:
    #     for i in range(F.shape[0]):
    #         df = pd.DataFrame(F[i, 1:, :])  # fatia i -> 2D
    #         df.to_excel(writer, sheet_name=NameCaso[i], index=False, header=False)
    # with pd.ExcelWriter("F_axial.xlsx", engine="openpyxl") as writer:
    #     for i in range(F_axial.shape[0]):
    #         df = pd.DataFrame(F_axial[i, 1:])  # fatia i -> 2D
    #         df.to_excel(writer, sheet_name=NameCaso[i], index=False, header=False)
    # with pd.ExcelWriter("FTW.xlsx", engine="openpyxl") as writer:
    #     for i in range(FTW.shape[0]):
    #         df = pd.DataFrame(FTW[i, 1:])  # fatia i -> 2D
    #         df.to_excel(writer, sheet_name=NameCaso[i], index=False, header=False)
    # with pd.ExcelWriter("FTW_dist.xlsx", engine="openpyxl") as writer:
    #     for i in range(FTW_dist.shape[0]):
    #         df = pd.DataFrame(FTW_dist[i, 1:])  # fatia i -> 2D
    #         df.to_excel(writer, sheet_name=NameCaso[i], index=False, header=False)
    # with pd.ExcelWriter("FTW_lin.xlsx", engine="openpyxl") as writer:
    #     for i in range(FTW_lin.shape[0]):
    #         df = pd.DataFrame(FTW_lin[i, 1:])  # fatia i -> 2D
    #         df.to_excel(writer, sheet_name=NameCaso[i], index=False, header=False)
    # with pd.ExcelWriter("FTW_lin_dist.xlsx", engine="openpyxl") as writer:
    #     for i in range(FTW_lin_dist.shape[0]):
    #         df = pd.DataFrame(FTW_lin_dist[i, 1:])  # fatia i -> 2D
    #         df.to_excel(writer, sheet_name=NameCaso[i], index=False, header=False)
    # with pd.ExcelWriter("F_nos.xlsx", engine="openpyxl") as writer:
    #     for i in range(F_nos.shape[0]):
    #         df = pd.DataFrame(F_nos[i, 1:,:])  # fatia i -> 2D
    #         df.to_excel(writer, sheet_name=NameCaso[i], index=False, header=False)
    # with pd.ExcelWriter("Reactions.xlsx", engine="openpyxl") as writer:
    #     for i in range(Reactions.shape[0]):
    #         df = pd.DataFrame(Reactions[i, 1:,:])  # fatia i -> 2D
    #         df.to_excel(writer, sheet_name=NameCaso[i], index=False, header=False)

    # with pd.ExcelWriter("CsCd.xlsx", engine="openpyxl") as writer:
    #     df = pd.DataFrame(CsCd)  # fatia i -> 2D
    #     df.to_excel(writer, sheet_name=f"CsCd", index=False, header=False)
    # with pd.ExcelWriter("Cfs.xlsx", engine="openpyxl") as writer:
    #     df = pd.DataFrame(cfs[1:,])  # fatia i -> 2D
    #     df.to_excel(writer, sheet_name=f"Cfs", index=False, header=False)
    # with pd.ExcelWriter("F_Tracao.xlsx", engine="openpyxl") as writer:
    #     df = pd.DataFrame(F_Tracao[1:])  # fatia i -> 2D
    #     df.to_excel(writer, sheet_name=f"F_Tracao", index=False, header=False)
    # with pd.ExcelWriter("F_Compressao.xlsx", engine="openpyxl") as writer:
    #     df = pd.DataFrame(F_Compressao[1:])  # fatia i -> 2D
    #     df.to_excel(writer, sheet_name=f"F_Compressao", index=False, header=False)

    # if Gelo=="Sim":
    #     with pd.ExcelWriter("Desloc_gelo.xlsx", engine="openpyxl") as writer:
    #         for i in range(Desloc_gelo.shape[0]):
    #             df = pd.DataFrame(Desloc_gelo[i, 1:, :])  # fatia i -> 2D
    #             df.to_excel(writer, sheet_name=NameCaso_gelo[i], index=False, header=False)
    #     with pd.ExcelWriter("F_gelo.xlsx", engine="openpyxl") as writer:
    #         for i in range(F_gelo.shape[0]):
    #             df = pd.DataFrame(F_gelo[i, 1:, :])  # fatia i -> 2D
    #             df.to_excel(writer, sheet_name=NameCaso_gelo[i], index=False, header=False)
    #     with pd.ExcelWriter("F_axial_gelo.xlsx", engine="openpyxl") as writer:
    #         for i in range(F_axial_gelo.shape[0]):
    #             df = pd.DataFrame(F_axial_gelo[i, 1:])  # fatia i -> 2D
    #             df.to_excel(writer, sheet_name=NameCaso_gelo[i], index=False, header=False)
    #     with pd.ExcelWriter("FTW_gelo.xlsx", engine="openpyxl") as writer:
    #         for i in range(FTW_gelo.shape[0]):
    #             df = pd.DataFrame(FTW_gelo[i, 1:])  # fatia i -> 2D
    #             df.to_excel(writer, sheet_name=NameCaso_gelo[i], index=False, header=False)
    #     with pd.ExcelWriter("FTW_dist_gelo.xlsx", engine="openpyxl") as writer:
    #         for i in range(FTW_dist_gelo.shape[0]):
    #             df = pd.DataFrame(FTW_dist_gelo[i, 1:])  # fatia i -> 2D
    #             df.to_excel(writer, sheet_name=NameCaso_gelo[i], index=False, header=False)
    #     with pd.ExcelWriter("FTW_lin_gelo.xlsx", engine="openpyxl") as writer:
    #         for i in range(FTW_lin_gelo.shape[0]):
    #             df = pd.DataFrame(FTW_lin_gelo[i, 1:])  # fatia i -> 2D
    #             df.to_excel(writer, sheet_name=NameCaso_gelo[i], index=False, header=False)
    #     with pd.ExcelWriter("FTW_lin_dist_gelo.xlsx", engine="openpyxl") as writer:
    #         for i in range(FTW_lin_dist_gelo.shape[0]):
    #             df = pd.DataFrame(FTW_lin_dist_gelo[i, 1:])  # fatia i -> 2D
    #             df.to_excel(writer, sheet_name=NameCaso_gelo[i], index=False, header=False)
    #     with pd.ExcelWriter("F_nos_gelo.xlsx", engine="openpyxl") as writer:
    #         for i in range(F_nos_gelo.shape[0]):
    #             df = pd.DataFrame(F_nos_gelo[i, 1:,:])  # fatia i -> 2D
    #             df.to_excel(writer, sheet_name=NameCaso_gelo[i], index=False, header=False)
    #     with pd.ExcelWriter("Reactions_gelo.xlsx", engine="openpyxl") as writer:
    #         for i in range(Reactions_gelo.shape[0]):
    #             df = pd.DataFrame(Reactions_gelo[i, 1:,:])  # fatia i -> 2D
    #             df.to_excel(writer, sheet_name=NameCaso_gelo[i], index=False, header=False)

    #     with pd.ExcelWriter("CsCd_gelo.xlsx", engine="openpyxl") as writer:
    #         df = pd.DataFrame(CsCd_gelo)  # fatia i -> 2D
    #         df.to_excel(writer, sheet_name=f"CsCd_gelo", index=False, header=False)
    #     with pd.ExcelWriter("Cfs_gelo.xlsx", engine="openpyxl") as writer:
    #         df = pd.DataFrame(cfs_gelo[1:,])  # fatia i -> 2D
    #         df.to_excel(writer, sheet_name=f"Cfs_gelo", index=False, header=False)


    # with pd.ExcelWriter("Encurvadura.xlsx", engine="openpyxl") as writer:
    #     df = pd.DataFrame(Enc_Out[1:,:])  # fatia i -> 2D
    #     df.to_excel(
    #         writer,
    #         sheet_name="F_Encurvadura",
    #         index=False,
    #         header=[
    #             "TrussType","Dim","Esp","Comprimento_Enc", "Ratio_Class", "Lim_Class", "Class_Enc",
    #             "Lambda", "Lambda_Lim", "Lambda_Check", "k_enc","Lambda_Enc_eff",
    #             "alpha_Enc", "Phi_Enc", "Chi_Enc", "Aeff", "Sigma", "Red",
    #             "eta_Enc", "F_Compressao", "NbRd", "Ratio_Enc", "Check_Enc"
    #         ]
    #     )
    # with pd.ExcelWriter("Ligacao_Bolted.xlsx", engine="openpyxl") as writer:
    #     df = pd.DataFrame(Lig_Out[1:,:])  # fatia i -> 2D
    #     df.to_excel(
    #         writer,
    #         sheet_name="Ligacao",
    #         index=False,
    #         header=[
    #         "TrussType","L", "t","e1_gusset","p1_gusset","e2_gusset","p2_gusset","t_gusset","Bolt",
    #         "N_Bolt_Lig","d_Bolt","d_0","As_bolt","A_bolt","alpha_v_Bolt","FvRd","Ratio_Bolt_Lig","Check_Bolt",
    #         "fu_Lig","fy_Lig","alpha_d_end","alpha_d_inner","k1_edge",
    #         "alpha_b_end","alpha_b_inner","FbRd_Edge_Inner","FbRd_Edge_End","FbRd","Ratio_Esmag_Lig","Check_Esmag","Anv_block","Ant_block",
    #         "Veff1Rd","Veff2Rd","beta_Cant","Anet","NuRd","Ratio_Block_Lig","Check_Block"
    #         ]
    #     )
    # with pd.ExcelWriter("Shaft_Bolted.xlsx", engine="openpyxl") as writer:
    #     df = pd.DataFrame(B_Out[1:,:])  # fatia i -> 2D
    #     df.to_excel(
    #         writer,
    #         sheet_name="Shaft - Bolted",
    #         index=False,
    #         header=[
    #         "e1_B","p1_B","e2_B","p2_B","t_B","Bolt_B",
    #         "N_Bolt_B","d_Bolt_B","d_0_B","As_bolt_B","A_bolt_B","alpha_v_Bolt_B","FvRd_B","Ratio_Bolt_Lig_B","Check_Bolt_B",
    #         "fu_B","fy_B","alpha_d_end_B","alpha_d_inner_B","k1_edge_B",
    #         "alpha_b_end_B","alpha_b_inner_B","FbRd_Edge_Inner_B","FbRd_Edge_End_B","FbRd_B","Ratio_Esmag_B","Check_Esmag_B","Anv_block_B","Ant_block_B",
    #         "Veff1Rd_B","Veff2Rd_B","beta_Cant_B","Anet_B","NuRd_B","Ratio_Block_B","Check_Block_B"
    #         ]
    #     )

    # with pd.ExcelWriter("Shaft_Flanged.xlsx", engine="openpyxl") as writer:
    #     df = pd.DataFrame(F_Out[1:, :])
    #     df.to_excel(
    #         writer,
    #         sheet_name="Shaft - Flanged",
    #         index=False,
    #         header=[
    #             "Ext_Per_F", "D_between_Bolt_F", "t_F", "Steel_F", "Bolt_F",
    #             "N_Bolt_F", "Bolt_Steel_F", "Stiffeners_F", "N_Stiff_F",
    #             "t_Stiff_F", "L_Stiff_F", "h_Stiff_F", "D_Ext_F",
    #             "D_Bolt_Flange_F", "D_int_F", "t_Tube_F", "aws_F",
    #             "e1_F", "e2_F", "alpha_r0_F", "d_Bolt_F", "Lb_F", "leff",
    #             "alpha_r_F", "epsilon_F", "fub_F", "fyb_F", "n_F", "m_F",
    #             "elinha_F", "FtRd", "Rb_F", "R_F", "Re_F", "k3_F",
    #             "NT1Rd", "NT2Rd", "NT3Rd", "NT4Rd","BpRd", "NTRd","Ratio_F","Check_F"
    #         ])


        # === 9. Plot deformada estática ===
    #scale_def = 100
    # fig = plt.figure(figsize=(8,6))
    # ax = fig.add_subplot(111, projection='3d')
    # for i in range(1,ID_Lig_Bolted+1):
    #     n1, n2 = ops.eleNodes(int(Ele_Lig[i]))

    #     xi, yi, zi =  ops.nodeCoord(n1)
    #     xj, yj, zj = ops.nodeCoord(n2)
    #     uxi, uyi, uzi = desloc[ni]
    #    uxj, uyj, uzj = desloc[nj]
    #     ax.plot([xi, xj], [yi, yj], [zi, zj], 'b--', label='Original')
    #    ax.plot([xi+uxi*scale_def, xj+uxj*scale_def],
    #            [yi+uyi*scale_def, yj+uyj*scale_def],
    #            [zi+uzi*scale_def, zj+uzj*scale_def], 'r-', label='Deformada' if _ == 0 else "")

    # ax.set_title("Deformada estática (ampliada)")
    # ax.set_xlabel("X")
    # ax.set_ylabel("Y")
    # ax.set_zlabel("Z")

    # plt.show()
    # ops.reactions()
    # reaction_1 = ops.nodeReaction(1)
    # print(f"Reações no nó 1: Rx = {reaction_1[0]} N, Ry = {reaction_1[1]} N, Rz = {reaction_1[2]} N")
    # reaction_1 = ops.nodeReaction(25)
    # print(f"Reações no nó 25: Rx = {reaction_1[0]} N, Ry = {reaction_1[1]} N, Rz = {reaction_1[2]} N")
    # reaction_1 = ops.nodeReaction(37)
    # print(f"Reações no nó 37: Rx = {reaction_1[0]} N, Ry = {reaction_1[1]} N, Rz = {reaction_1[2]} N")
    ops.wipe()
    return np.abs(Ratio_max_optimization-1), Final, OutputEstrutural,Ratio_max_optimization


f_obj, Final,OutputEstrutural,Ratio_max_optimization=main(0)


if Final==1 and f_obj>10**-4 and Ratio_max_optimization<1:

    def func_to_minimize(x):
        f_obj, Final, a,b = main(x[0])
        
        return f_obj

    res = minimize(
        func_to_minimize,
        [0.3],
        method='Nelder-Mead',
        options={
            'xatol': 1e-2,   # tolerância em x (antes era 1e-6)
            'fatol': 1e-4,   # tolerância em f(x)
            'maxiter': 500,
        }
    )
    Aincrease = res.x[0]
    Ratio_opt, Final_opt, Out_opt,b = main(Aincrease)
    
else:
    Aincrease=0
if Final==1:
    with open(OutputEstrutural, mode='a', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(["Area Left", np.round(Aincrease,2)])
 