import numpy as np
from Wind.Wind_Ec import ZoneInverse
def Desloc_Output(Desloc,nos_grupo,iSLS,iELU,iELU_gelo,iEnd,ivento,z_grupo,h_torre,alpha_vector,v_no_maximo,Ncasoscarga,todos_nos,Gamma_Perm_Des,Gamma_Variable_Des,Gamma_Perm_Fav,Gelo,Desloc_gelo,Ncasoscarga_gelo,icasoscarga_gelo,k_gelo,Psiventogelo): 
    Desloc_aux=np.sqrt(Desloc[:,:,0]**2+Desloc[:,:,1]**2)


    ##Topo
    Desloc_top=np.mean(Desloc_aux[:,nos_grupo[0,:].astype(int)], axis=1)

    Desloc_top2=np.mean(Desloc_aux[:,nos_grupo[1,:].astype(int)], axis=1)
    Desloctop_SLS=np.max(Desloc_top[iSLS:iELU])
    Desloctop2_SLS=np.max(Desloc_top2[iSLS:iELU])

    Rottop_SLS=np.rad2deg(np.arctan((Desloctop_SLS-Desloctop2_SLS)/(z_grupo[0]-z_grupo[1])))
    Desloctop_ELU=np.max(Desloc_top[iELU:iEnd])
    iAngle_Crit=np.argmax(Desloc_top[iELU:iEnd])

    Desloctop2_ELU=np.max(Desloc_top2[iELU:iEnd])

    Rottop_ELU=np.rad2deg(np.arctan((Desloctop_ELU-Desloctop2_ELU)/(z_grupo[0]-z_grupo[1])))


    ##2/3
    # índices de elementos maiores e menores que o valor
    superiores = np.where(2/3*h_torre > z_grupo)[0]  # elementos maiores que valor
    inferiores = np.where(2/3*h_torre  <= z_grupo)[0]  # elementos menores que valor
    # pegar os mais próximos
    idx_superior = superiores[0]  # último que é maior
    idx_inferior = inferiores[-1]   # primeiro que é menor
    Desloc_23=np.mean(Desloc_aux[:,nos_grupo[idx_superior,:].astype(int)], axis=1)
    Desloc_232=np.mean(Desloc_aux[:,nos_grupo[idx_inferior,:].astype(int)], axis=1)
    Desloc23_SLS=np.max(Desloc_23[iSLS:iELU])
    Desloc232_SLS=np.max(Desloc_232[iSLS:iELU])
    Rot23_SLS=np.rad2deg(np.arctan((Desloc23_SLS-Desloc232_SLS)/(z_grupo[idx_superior]-z_grupo[idx_inferior])))

    Desloc23_ELU=np.max(Desloc_23[iELU:iEnd])
    Desloc232_ELU=np.max(Desloc_232[iELU:iEnd])

    Rot23_ELU=np.rad2deg(np.arctan((Desloc23_ELU-Desloc232_ELU)/(z_grupo[idx_superior]-z_grupo[idx_inferior])))
    Desloc_100=np.zeros((Ncasoscarga,len(todos_nos)+2,6))
    Desloc_120=np.zeros((Ncasoscarga,len(todos_nos)+2,6))
    Ratio100=100*1000/3600/np.max(v_no_maximo)
    Ratio120=120*1000/3600/np.max(v_no_maximo)
    Desloc_100[:4,:,:]=Desloc[:4,:,:]
    Desloc_120[:4,:,:]=Desloc[:4,:,:]
    Desloc_100[4:iSLS,:,:]=Desloc[4:iSLS,:,:]*Ratio100**2
    Desloc_120[4:iSLS,:,:]=Desloc[4:iSLS,:,:]*Ratio120**2
    i100=iSLS
    for i in range(len(alpha_vector)):
        Desloc_100[i100+i,:,:]=1.0*Desloc_100[1,:,:]+1.0*Desloc_100[2,:,:]+1.0*Desloc_100[3,:,:]+1.0*Desloc_100[ivento+i+1,:,:]
        Desloc_120[i100+i,:,:]=1.0*Desloc_120[1,:,:]+1.0*Desloc_120[2,:,:]+1.0*Desloc_120[3,:,:]+1.0*Desloc_120[ivento+i+1,:,:]
    iholder_100=i100+i   
    for i in range(len(alpha_vector)):
        Desloc_100[iholder_100+i+1,:,:]=Gamma_Perm_Des*(Desloc_100[1,:,:]+Desloc_100[2,:,:]+Desloc_100[3,:,:])+Gamma_Variable_Des*Desloc_100[ivento+i+1,:,:]
        Desloc_120[iholder_100+i+1,:,:]=Gamma_Perm_Des*(Desloc_120[1,:,:]+Desloc_120[2,:,:]+Desloc_120[3,:,:])+Gamma_Variable_Des*Desloc_120[ivento+i+1,:,:]
    iholder_100=iholder_100+i   
    for i in range(len(alpha_vector)):
        Desloc_100[iholder_100+i+1,:,:]=Gamma_Perm_Fav*(Desloc_100[1,:,:]+Desloc_100[2,:,:]+Desloc_100[3,:,:])+Gamma_Variable_Des*Desloc_100[ivento+i+1,:,:]
        Desloc_120[iholder_100+i+1,:,:]=Gamma_Perm_Fav*(Desloc_120[1,:,:]+Desloc_120[2,:,:]+Desloc_120[3,:,:])+Gamma_Variable_Des*Desloc_120[ivento+i+1,:,:]


    Desloc_100_aux=np.sqrt(Desloc_100[:,:,0]**2+Desloc_100[:,:,1]**2)
    Desloc_120_aux=np.sqrt(Desloc_120[:,:,0]**2+Desloc_120[:,:,1]**2)

    Desloc_100_top=np.mean(Desloc_100_aux[:,nos_grupo[0,:].astype(int)], axis=1)
    Desloc_120_top=np.mean(Desloc_120_aux[:,nos_grupo[0,:].astype(int)], axis=1)
    Desloc_100_top_SLS=np.max(Desloc_100_top[iSLS:iELU])
    Desloc_120_top_SLS=np.max(Desloc_120_top[iSLS:iELU])

    Desloc_100_top2=np.mean(Desloc_100_aux[:,nos_grupo[1,:].astype(int)], axis=1)
    Desloc_120_top2=np.mean(Desloc_120_aux[:,nos_grupo[1,:].astype(int)], axis=1)

    Desloc_100_top2_SLS=np.max(Desloc_100_top2[iSLS:iELU])
    Desloc_120_top2_SLS=np.max(Desloc_120_top2[iSLS:iELU])
    Rot_100_top_SLS=np.rad2deg(np.arctan((Desloc_100_top_SLS-Desloc_100_top2_SLS)/(z_grupo[0]-z_grupo[1])))
    Rot_120_top_SLS=np.rad2deg(np.arctan((Desloc_120_top_SLS-Desloc_120_top2_SLS)/(z_grupo[0]-z_grupo[1])))

    Desloc_100_top_ULS=np.max(Desloc_100_top[iELU:])
    Desloc_120_top_ULS=np.max(Desloc_120_top[iELU:])
    Desloc_100_top2_ULS=np.max(Desloc_100_top2[iELU:])
    Desloc_120_top2_ULS=np.max(Desloc_120_top2[iELU:])
    Rot_100_top_ULS=np.rad2deg(np.arctan((Desloc_100_top_ULS-Desloc_100_top2_ULS)/(z_grupo[0]-z_grupo[1])))
    Rot_120_top_ULS=np.rad2deg(np.arctan((Desloc_120_top_ULS-Desloc_120_top2_ULS)/(z_grupo[0]-z_grupo[1])))

    if Gelo=="Sim":
        Desloc_100_gelo=np.zeros((Ncasoscarga_gelo,len(todos_nos)+2,6))
        Desloc_120_gelo=np.zeros((Ncasoscarga_gelo,len(todos_nos)+2,6))
        Desloc_aux_gelo=np.sqrt(Desloc_gelo[:,:,0]**2+Desloc_gelo[:,:,1]**2)


        ##Topo
        Desloc_top_gelo=np.mean(Desloc_aux_gelo[:,nos_grupo[0,:].astype(int)], axis=1)

        Desloc_top2_gelo=np.mean(Desloc_aux_gelo[:,nos_grupo[1,:].astype(int)], axis=1)
        Desloctop_SLS_gelo=np.max(Desloc_top_gelo[iSLS:iELU_gelo])
        Desloctop2_SLS_gelo=np.max(Desloc_top2_gelo[iSLS:iELU_gelo])

        Rottop_SLS_gelo=np.rad2deg(np.arctan((Desloctop_SLS_gelo-Desloctop2_SLS_gelo)/(z_grupo[0]-z_grupo[1])))
        Desloctop_ELU_gelo=np.max(Desloc_top_gelo[iELU_gelo:])
        Desloctop2_ELU_gelo=np.max(Desloc_top2_gelo[iELU_gelo:])
        Rottop_ELU_gelo=np.rad2deg(np.arctan((Desloctop_ELU_gelo-Desloctop2_ELU_gelo)/(z_grupo[0]-z_grupo[1])))

        iAngle_Crit_gelo=np.argmax(Desloc_top_gelo[iELU_gelo:])
        if Desloctop_ELU_gelo>Desloctop_ELU:
            iAngle_Crit=iAngle_Crit_gelo

        ##2/3
        # índices de elementos maiores e menores que o valor
        superiores = np.where(2/3*h_torre > z_grupo)[0]  # elementos maiores que valor
        inferiores = np.where(2/3*h_torre  <= z_grupo)[0]  # elementos menores que valor

        # pegar os mais próximos
        idx_superior = superiores[0]  # último que é maior
        idx_inferior = inferiores[-1]   # primeiro que é menor
        Desloc_23_gelo=np.mean(Desloc_aux_gelo[:,nos_grupo[idx_superior,:].astype(int)], axis=1)

        Desloc_232_gelo=np.mean(Desloc_aux_gelo[:,nos_grupo[idx_inferior,:].astype(int)], axis=1)

        Desloc23_SLS_gelo=np.max(Desloc_23_gelo[iSLS:iELU_gelo])
        Desloc232_SLS_gelo=np.max(Desloc_232_gelo[iSLS:iELU_gelo])

        Rot23_SLS_gelo=np.rad2deg(np.arctan((Desloc23_SLS_gelo-Desloc232_SLS_gelo)/(z_grupo[idx_superior]-z_grupo[idx_inferior])))
        Desloc23_ELU_gelo=np.max(Desloc_23_gelo[iELU_gelo:])
        Desloc232_ELU_gelo=np.max(Desloc_232_gelo[iELU_gelo:])

        Rot23_ELU_gelo=np.rad2deg(np.arctan((Desloc23_ELU_gelo-Desloc232_ELU_gelo)/(z_grupo[idx_superior]-z_grupo[idx_inferior])))


        Ratio100=100*1000/3600/np.max(v_no_maximo)

        Ratio120=120*1000/3600/np.max(v_no_maximo)
        Desloc_100_gelo[:4,:,:]=Desloc_gelo[:4,:,:]
        Desloc_120_gelo[:4,:,:]=Desloc_gelo[:4,:,:]
        Desloc_100_gelo[4:iSLS,:,:]=Desloc_gelo[4:iSLS,:,:]*Ratio100**2
        Desloc_120_gelo[4:iSLS,:,:]=Desloc_gelo[4:iSLS,:,:]*Ratio120**2
        i100=iSLS
        for i in range(len(alpha_vector)):
            Desloc_100_gelo[i100+i,:,:]=1.0*(Desloc[1,:,:]+Desloc[2,:,:]+Desloc[3,:,:])+k_gelo*Desloc_100_gelo[icasoscarga_gelo+i+1,:,:]+Psiventogelo*(Desloc_100_gelo[1,:,:]+Desloc_100_gelo[2,:,:]+Desloc_100_gelo[3,:,:])
            Desloc_120_gelo[i100+i,:,:]=1.0*(Desloc[1,:,:]+Desloc[2,:,:]+Desloc[3,:,:])+k_gelo*Desloc_120_gelo[icasoscarga_gelo+i+1,:,:]+Psiventogelo*(Desloc_120_gelo[1,:,:]+Desloc_120_gelo[2,:,:]+Desloc_120_gelo[3,:,:])
        iholder_100=i100+i  
        for i in range(len(alpha_vector)):
            #### ELS - Gelo Dominante
            Desloc_100_gelo[iholder_100+i+1,:,:]=1*(Desloc[1,:,:]+Desloc[2,:,:]+Desloc[3,:,:])+k_gelo*Psiventogelo*Desloc_100_gelo[icasoscarga_gelo+i+1,:,:]+(Desloc_100_gelo[1,:,:]+Desloc_100_gelo[2,:,:]+Desloc_100_gelo[3,:,:])
            Desloc_120_gelo[iholder_100+i+1,:,:]=1*(Desloc[1,:,:]+Desloc[2,:,:]+Desloc[3,:,:])+k_gelo*Psiventogelo*Desloc_120_gelo[icasoscarga_gelo+i+1,:,:]+(Desloc_120_gelo[1,:,:]+Desloc_120_gelo[2,:,:]+Desloc_120_gelo[3,:,:])
        iholder_100=iholder_100+i   
        for i in range(len(alpha_vector)):
            ### ELU-Vento Dominante ###
            Desloc_100_gelo[iholder_100+i+1,:,:]=Gamma_Perm_Des*(Desloc[1,:,:]+Desloc[2,:,:]+Desloc[3,:,:])+Gamma_Variable_Des*k_gelo*Desloc_100_gelo[icasoscarga_gelo+i+1,:,:]+Gamma_Variable_Des*Psiventogelo*(Desloc_100_gelo[1,:,:]+Desloc_100_gelo[2,:,:]+Desloc_100_gelo[3,:,:])
            Desloc_120_gelo[iholder_100+i+1,:,:]=Gamma_Perm_Des*(Desloc[1,:,:]+Desloc[2,:,:]+Desloc[3,:,:])+Gamma_Variable_Des*k_gelo*Desloc_120_gelo[icasoscarga_gelo+i+1,:,:]+Gamma_Variable_Des*Psiventogelo*(Desloc_120_gelo[1,:,:]+Desloc_120_gelo[2,:,:]+Desloc_120_gelo[3,:,:])
        iholder_100=iholder_100+i+1
        for i in range(len(alpha_vector)):
            ### ELU - Gelo Dominante ###
            Desloc_100_gelo[iholder_100+i+1,:,:]=Gamma_Perm_Des*(Desloc[1,:,:]+Desloc[2,:,:]+Desloc[3,:,:])+Gamma_Variable_Des*Psiventogelo*k_gelo*Desloc_100_gelo[icasoscarga_gelo+i+1,:,:]+Gamma_Variable_Des*(Desloc_100_gelo[1,:,:]+Desloc_100_gelo[2,:,:]+Desloc_100_gelo[3,:,:])
            Desloc_120_gelo[iholder_100+i+1,:,:]=Gamma_Perm_Des*(Desloc[1,:,:]+Desloc[2,:,:]+Desloc[3,:,:])+Gamma_Variable_Des*Psiventogelo*k_gelo*Desloc_120_gelo[icasoscarga_gelo+i+1,:,:]+Gamma_Variable_Des*(Desloc_120_gelo[1,:,:]+Desloc_120_gelo[2,:,:]+Desloc_120_gelo[3,:,:])
        iholder_100=iholder_100+i+1
        for i in range(len(alpha_vector)):
            ### ELU-Vento Dominante ###
            Desloc_100_gelo[iholder_100+i+1,:,:]=Gamma_Perm_Fav*(Desloc[1,:,:]+Desloc[2,:,:]+Desloc[3,:,:])+Gamma_Variable_Des*k_gelo*Desloc_100_gelo[icasoscarga_gelo+i+1,:,:]+Gamma_Variable_Des*Psiventogelo*(Desloc_100_gelo[1,:,:]+Desloc_100_gelo[2,:,:]+Desloc_100_gelo[3,:,:])
            Desloc_120_gelo[iholder_100+i+1,:,:]=Gamma_Perm_Fav*(Desloc[1,:,:]+Desloc[2,:,:]+Desloc[3,:,:])+Gamma_Variable_Des*k_gelo*Desloc_120_gelo[icasoscarga_gelo+i+1,:,:]+Gamma_Variable_Des*Psiventogelo*(Desloc_120_gelo[1,:,:]+Desloc_120_gelo[2,:,:]+Desloc_120_gelo[3,:,:])
        iholder_100=iholder_100+i+1
        for i in range(len(alpha_vector)):
            ### ELU - Gelo Dominante ###
            Desloc_100_gelo[iholder_100+i+1,:,:]=Gamma_Perm_Fav*(Desloc[1,:,:]+Desloc[2,:,:]+Desloc[3,:,:])+Gamma_Variable_Des*Psiventogelo*k_gelo*Desloc_100_gelo[icasoscarga_gelo+i+1,:,:]+Gamma_Variable_Des*(Desloc_100_gelo[1,:,:]+Desloc_100_gelo[2,:,:]+Desloc_100_gelo[3,:,:])
            Desloc_120_gelo[iholder_100+i+1,:,:]=Gamma_Perm_Fav*(Desloc[1,:,:]+Desloc[2,:,:]+Desloc[3,:,:])+Gamma_Variable_Des*Psiventogelo*k_gelo*Desloc_120_gelo[icasoscarga_gelo+i+1,:,:]+Gamma_Variable_Des*(Desloc_120_gelo[1,:,:]+Desloc_120_gelo[2,:,:]+Desloc_120_gelo[3,:,:])
        iholder_100=iholder_100+i+1  
        Desloc_100_aux_gelo=np.sqrt(Desloc_100_gelo[:,:,0]**2+Desloc_100_gelo[:,:,1]**2)
        Desloc_120_aux_gelo=np.sqrt(Desloc_120_gelo[:,:,0]**2+Desloc_120_gelo[:,:,1]**2)

        Desloc_100_top_gelo=np.mean(Desloc_100_aux_gelo[:,nos_grupo[0,:].astype(int)], axis=1)
        Desloc_120_top_gelo=np.mean(Desloc_120_aux_gelo[:,nos_grupo[0,:].astype(int)], axis=1)
        Desloc_100_top_SLS_gelo=np.max(Desloc_100_top_gelo[iSLS:iELU_gelo])
        Desloc_120_top_SLS_gelo=np.max(Desloc_120_top_gelo[iSLS:iELU_gelo])

        Desloc_100_top2_gelo=np.mean(Desloc_100_aux_gelo[:,nos_grupo[1,:].astype(int)], axis=1)
        Desloc_120_top2_gelo=np.mean(Desloc_120_aux_gelo[:,nos_grupo[1,:].astype(int)], axis=1)

        Desloc_100_top2_SLS_gelo=np.max(Desloc_100_top2_gelo[iSLS:iELU_gelo])
        Desloc_120_top2_SLS_gelo=np.max(Desloc_120_top2_gelo[iSLS:iELU_gelo])
        Rot_100_top_SLS_gelo=np.rad2deg(np.arctan((Desloc_100_top_SLS_gelo-Desloc_100_top2_SLS_gelo)/(z_grupo[0]-z_grupo[1])))
        Rot_120_top_SLS_gelo=np.rad2deg(np.arctan((Desloc_120_top_SLS_gelo-Desloc_120_top2_SLS_gelo)/(z_grupo[0]-z_grupo[1])))

        Desloc_100_top_ULS_gelo=np.max(Desloc_100_top_gelo[iELU_gelo:])
        Desloc_120_top_ULS_gelo=np.max(Desloc_120_top_gelo[iELU_gelo:])
        Desloc_100_top2_ULS_gelo=np.max(Desloc_100_top2_gelo[iELU_gelo:])
        Desloc_120_top2_ULS_gelo=np.max(Desloc_120_top2_gelo[iELU_gelo:])
        Rot_100_top_ULS_gelo=np.rad2deg(np.arctan((Desloc_100_top_ULS_gelo-Desloc_100_top2_ULS_gelo)/(z_grupo[0]-z_grupo[1])))
        Rot_120_top_ULS_gelo=np.rad2deg(np.arctan((Desloc_120_top_ULS_gelo-Desloc_120_top2_ULS_gelo)/(z_grupo[0]-z_grupo[1])))
        Desloc_Out = np.array([
            ["Case", "Degree", "Minute", "Distance"],
            ["SLS - Top",
                np.round(np.max([Rottop_SLS, Rottop_SLS_gelo]), 3),
                np.round(np.max([Rottop_SLS, Rottop_SLS_gelo]) * 60, 3),
                np.round(np.max([Desloctop_SLS, Desloctop_SLS_gelo]) * 1e3, 3)
            ],
            ["SLS - 2/3",
                np.round(np.max([Rot23_SLS, Rot23_SLS_gelo]), 3),
                np.round(np.max([Rot23_SLS, Rot23_SLS_gelo]) * 60, 3),
                np.round(np.max([Desloc232_SLS, Desloc232_SLS_gelo]) * 1e3, 3)
            ],
            ["ULS - Top",
                np.round(np.max([Rottop_ELU, Rottop_ELU_gelo]), 3),
                np.round(np.max([Rottop_ELU, Rottop_ELU_gelo]) * 60, 3),
                np.round(np.max([Desloctop_ELU, Desloctop_ELU_gelo]) * 1e3, 3)
            ],
            ["ULS - 2/3",
                np.round(np.max([Rot23_ELU, Rot23_ELU_gelo]), 3),
                np.round(np.max([Rot23_ELU, Rot23_ELU_gelo]) * 60, 3),
                np.round(np.max([Desloc232_ELU, Desloc232_ELU_gelo]) * 1e3, 3)
            ],
            ["SLS - 100 km/h",
                np.round(np.max([Rot_100_top_SLS, Rot_100_top_SLS_gelo]), 3),
                np.round(np.max([Rot_100_top_SLS, Rot_100_top_SLS_gelo]) * 60, 3),
                np.round(np.max([Desloc_100_top_SLS, Desloc_100_top_SLS_gelo]) * 1e3, 3)
            ],

            ["SLS - 120 km/h",
                np.round(np.max([Rot_120_top_SLS, Rot_120_top_SLS_gelo]), 3),
                np.round(np.max([Rot_120_top_SLS, Rot_120_top_SLS_gelo]) * 60, 3),
                np.round(np.max([Desloc_120_top_SLS, Desloc_120_top_SLS_gelo]) * 1e3, 3)
            ]
        ], dtype=object)
        # ["ULS - 120 km/h",
        #         np.round(np.max([Rot_120_top_ULS, Rot_120_top_ULS_gelo]), 3),
        #         np.round(np.max([Rot_120_top_ULS, Rot_120_top_ULS_gelo]) * 60, 3),
        #         np.round(np.max([Desloc_120_top_ULS, Desloc_120_top_ULS_gelo]) * 1e3, 3)
        #     ]
        #         ["ULS - 100 km/h",
        #     np.round(np.max([Rot_100_top_ULS, Rot_100_top_ULS_gelo]), 3),
        #     np.round(np.max([Rot_100_top_ULS, Rot_100_top_ULS_gelo]) * 60, 3),
        #     np.round(np.max([Desloc_100_top_ULS, Desloc_100_top_ULS_gelo]) * 1e3, 3)
        # ],
    else:
        Desloc_Out= np.array([["Case", "Degree", "Minute", "Distance"],
                        ["SLS - Top", np.round(Rottop_SLS,3)  ,np.round(Rottop_SLS*60,3) ,np.round(Desloctop_SLS*10**3,1) ],
                        ["SLS - 2/3", np.round(Rot23_SLS,3)  ,np.round(Rot23_SLS*60,3) ,np.round(Desloc232_SLS*10**3,1) ],
                        ["ULS - Top", np.round(Rottop_ELU,3)  ,np.round(Rottop_ELU*60,3) ,np.round(Desloctop_ELU*10**3,1) ],
                        ["ULS - 2/3", np.round(Rot23_ELU,3)  ,np.round(Rot23_ELU*60,3) ,np.round(Desloc232_ELU*10**3,1) ],
                        ["SLS - 100 km/h", np.round(Rot_100_top_SLS,3)  ,np.round(Rot_100_top_SLS*60,3) ,np.round(Desloc_100_top_SLS*10**3,1) ],
                        ["SLS - 120 km/h", np.round(Rot_120_top_SLS,3)  ,np.round(Rot_120_top_SLS*60,3) ,np.round(Desloc_120_top_SLS*10**3,1) ]], dtype=object)
        #["ULS - 120 km/h", np.round(Rot_120_top_ULS,3)  ,np.round(Rot_120_top_ULS*60,3) ,np.round(Desloc_120_top_ULS*10**3,1) ]
        #["ULS - 100 km/h", np.round(Rot_100_top_ULS,3)  ,np.round(Rot_100_top_ULS*60,3) ,np.round(Desloc_100_top_ULS*10**3,1) ],
    if Gelo=="Sim":
        Rot_Top_SLS=np.max([Rottop_SLS, Rottop_SLS_gelo])
        Rot_23_SLS=np.max([Rot23_SLS, Rot23_SLS_gelo])
        Rot_Top_ELU=np.max([Rottop_ELU, Rottop_ELU_gelo])
        Rot_23_ELU=np.max([Rot23_ELU, Rot23_ELU_gelo])
        Rot_Top_100=np.max([Rot_100_top_SLS, Rot_100_top_SLS_gelo])
        Rot_Top_120=np.max([Rot_120_top_SLS, Rot_120_top_SLS_gelo])
    else:
        Rot_Top_SLS=Rottop_SLS
        Rot_23_SLS=np.max([Rot23_SLS])
        Rot_Top_ELU=np.max([Rottop_ELU])
        Rot_23_ELU=np.max([Rot23_ELU])
        Rot_Top_100=np.max([Rot_100_top_SLS])
        Rot_Top_120=np.max([Rot_120_top_SLS])
    return iAngle_Crit,Rot_Top_SLS,Rot_23_SLS,Rot_Top_ELU,Rot_23_ELU,Rot_Top_100,Rot_Top_120,Desloc_Out

def Buckling_Output(Class_Enc,nos_vento,TrussType,Truss_Out,Size_Profile,Elements_Matrix,Lambda_Enc_eff,Lambda_1,NTracRd,NbRd,Troco,Ratio_Enc,F_Compressao):
    Enc_Out_csv=np.array(["Member", "Profile", "Steel", "Slenderness","NtRd","NbRd","Ned","Ratio","Result"])
    ID_Ele_vec= np.arange(0, len(Class_Enc) + 1)

    for j in range(3):
        for i in range(1, nos_vento.shape[0]):
            # Agrupa quando j==2 ("Horizontal Bar")
            if j == 2:
                mask_type = (TrussType == "Horizontal Bar") | (TrussType == "External Manual Bar") | (TrussType == "Internal Manual Bar")
            else:
                mask_type = (TrussType == Truss_Out[j])

            i_Out_Zone = (Troco == i) & mask_type
            indices = np.where(i_Out_Zone)[0]

            if len(indices) != 0:
                iCrit_local = np.argmax(Ratio_Enc[indices])
                iCrit_abs = indices[iCrit_local]
                if Ratio_Enc[iCrit_abs]==0:
                    Check="NA"
                elif np.round(Ratio_Enc[iCrit_abs], 2)>1:
                    Check="KO"
                else:
                    Check="OK"

                # nome do grupo de saída
                Name_Enc_Out = f"{Truss_Out[j]} - Shaft { nos_vento.shape[0] -i}"
                row = np.array([
                    Name_Enc_Out,
                    Size_Profile[iCrit_abs],
                    Elements_Matrix[iCrit_abs - 1, 7],
                    round(Lambda_Enc_eff[iCrit_abs] * Lambda_1[iCrit_abs], 1),
                    round(NTracRd[iCrit_abs] * 1e-3, 1),
                    round(NbRd[iCrit_abs] * 1e-3, 1),
                    round(np.abs(F_Compressao[iCrit_abs]) * 1e-3, 1),
                    np.round(Ratio_Enc[iCrit_abs], 2),
                    Check
                ], dtype=object)

                Enc_Out_csv = np.vstack((Enc_Out_csv, row))
    return Enc_Out_csv

def Connection_Bracing_Output(TrussType,Troco,nos_vento, Truss_Out,Ratio_Lig,t_gusset,e1_gusset,e2_gusset,p1_gusset,p2_gusset,fy_Lig,N_Bolt_Lig,Bolt_Steel,Bolt,Fele,FLig_Rd,Ele_Lig,Check_Lig,Calc_Block,Calc_Esmag,Calc_Bolt,Calc_Nu):
    Lig_Out_csv=np.array(["Member","Plate thickness","e1","e2","p1","p2", "Steel Plate","Number of bolts","Bolts Size","Bolts steel","Ned","NRd","Ratio","Result"])
    TrussType_Lig_Out = TrussType[TrussType != "Leg"]
    Troco_Lig_Out = Troco[TrussType != "Leg"]
    for j in [1,2]:
        for i in range(1, nos_vento.shape[0]):
            # Agrupa quando j==2 ("Horizontal Bar")
            if j == 2:
                mask_type = (TrussType_Lig_Out == "Horizontal Bar") | (TrussType_Lig_Out == "External Manual Bar") | (TrussType_Lig_Out== "Internal Manual Bar")
            else:
                mask_type = (TrussType_Lig_Out == Truss_Out[j])

            i_Out_Zone = (Troco_Lig_Out == i) & mask_type
            indices = np.where(i_Out_Zone)[0]
            if len(indices) != 0:
                iCrit_local = np.argmax(Ratio_Lig[indices])
                iCrit_abs = indices[iCrit_local]
                Name_Enc_Out = f"{Truss_Out[j]} - Shaft { nos_vento.shape[0] -i}"
                if Check_Lig[iCrit_abs]!="NA":

                # if Ratio_Lig[iCrit_abs]==0:
                #     Check="NA"
                # elif np.round(Ratio_Lig[iCrit_abs], 2)>1:
                #     Check="KO"
                # else:
                #     Check="OK"
                    if Calc_Block[iCrit_abs]=="Nao" and Calc_Bolt[iCrit_abs]=="Sim" and Calc_Esmag[iCrit_abs]=="Nao" and Calc_Nu[iCrit_abs]=="Nao":
                        row = np.array([
                            Name_Enc_Out,
                            "--",
                            "--",
                            "--",
                            "--",
                            "--",
                            "--",
                            int(N_Bolt_Lig[iCrit_abs]),
                            Bolt_Steel[iCrit_abs],
                            Bolt[iCrit_abs],
                            np.round(Fele[int(Ele_Lig[iCrit_abs])]*10**-3,1),
                            np.round(FLig_Rd[iCrit_abs]*10**-3,1),
                            np.round(Ratio_Lig[iCrit_abs],2),
                            Check_Lig[iCrit_abs]
                            ], dtype=object)
                    elif Check_Lig[iCrit_abs]=="Welded":
                        row = np.array([
                            Name_Enc_Out,
                            "--",
                            "--",
                            "--",
                            "--",
                            "--",
                            "--",
                            "--",
                            "--",
                            "--",
                            "--",
                            "--",
                            "--",
                            Check_Lig[iCrit_abs]
                            ], dtype=object)
                    else:
                        row = np.array([
                            Name_Enc_Out,
                            np.round(t_gusset[iCrit_abs]*10**3,1),
                            np.round(e1_gusset[iCrit_abs]*10**3,1),
                            np.round(e2_gusset[iCrit_abs]*10**3,1),
                            np.round(p1_gusset[iCrit_abs]*10**3,1),
                            np.round(p2_gusset[iCrit_abs]*10**3,1),
                            f"S{int(fy_Lig[iCrit_abs]*10**-6)}",
                            int(N_Bolt_Lig[iCrit_abs]),
                            Bolt_Steel[iCrit_abs],
                            Bolt[iCrit_abs],
                            np.round(Fele[int(Ele_Lig[iCrit_abs])]*10**-3,1),
                            np.round(FLig_Rd[iCrit_abs]*10**-3,1),
                            np.round(Ratio_Lig[iCrit_abs],2),
                            Check_Lig[iCrit_abs]
                            ], dtype=object)
                    
                else:
                    row = np.array([
                        Name_Enc_Out,
                        "--",
                        "--",
                        "--",
                        "--",
                        "--",
                        "--",
                        "--",
                        "--",
                        "--",
                        "--",
                        "--",
                        "--",
                        Check_Lig[iCrit_abs]
                        ], dtype=object)
                Lig_Out_csv = np.vstack((Lig_Out_csv, row))
    return Lig_Out_csv

def Connection_Shaft_Output(TrussType,Troco,Ratio_Lig_B,nos_vento,divisor,t_B,e1_B,e2_B,p1_B,p2_B,Corte_B,fy_B,N_Bolt_B,Bolt_Steel_B,Bolt_B,F_Ele_B,FLig_Rd_B,Check_Lig_B,Calc_Bolt_B,Calc_Esmag_B,Calc_Block_B,Calc_Nu_B):

    TrussType_Lig_Out = TrussType[TrussType == "Leg"]
    Troco_Lig_Out = Troco[TrussType == "Leg"]


    B_Out_csv=np.array([["Member","Plate thickness","e1","e2","p1","p2","Number of plates", "Steel Plate","Number of bolts","Bolts Size","Bolts steel","NRd","Ned","Ratio","Result"]])
    if len(Ratio_Lig_B)!=1:
        for i in range(1, nos_vento.shape[0]):
            inicio = (i-1) * (divisor)+1
            fim = inicio + divisor
            bloco = Ratio_Lig_B[inicio:fim]
            i_local = np.argmax(bloco)           # posição dentro do bloco
            i_global = inicio + i_local   
            Name_Enc_Out = f"Leg - Shaft {i}"
            if Check_Lig_B[i_global]!="NA":
                if Calc_Block_B[i_global]=="Nao" and Calc_Bolt_B[i_global]=="Sim" and Calc_Esmag_B[i_global]=="Nao" and Calc_Nu_B[i_global]=="Nao":
                    row = np.array([
                        Name_Enc_Out,
                        "--",
                        "--",
                        "--",
                        "--",
                        "--",
                        int(Corte_B[i_global]),
                        "--",
                        int(N_Bolt_B[i_global]),
                        Bolt_Steel_B[i_global],
                        Bolt_B[i_global],
                        np.round(F_Ele_B[i_global]*10**-3,1),
                        np.round(FLig_Rd_B[i_global]*10**-3,1),
                        np.round(Ratio_Lig_B[i_global],2),
                        Check_Lig_B[i_global]
                        ], dtype=object)
                else:
                    row = np.array([
                        Name_Enc_Out,
                        np.round(t_B[i_global]*10**3,1),
                        np.round(e1_B[i_global]*10**3,1),
                        np.round(e2_B[i_global]*10**3,1),
                        np.round(p1_B[i_global]*10**3,1),
                        np.round(p2_B[i_global]*10**3,1),
                        int(Corte_B[i_global]),
                        f"S{int(fy_B[i_global]*10**-6)}",
                        int(N_Bolt_B[i_global]),
                        Bolt_Steel_B[i_global],
                        Bolt_B[i_global],
                        np.round(F_Ele_B[i_global]*10**-3,1),
                        np.round(FLig_Rd_B[i_global]*10**-3,1),
                        np.round(Ratio_Lig_B[i_global],2),
                        Check_Lig_B[i_global]
                        ], dtype=object)
            else:
                row = np.array([
                    Name_Enc_Out,
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    Check_Lig_B[i_global]
                    ], dtype=object)
            B_Out_csv = np.vstack((row,B_Out_csv))
        B_Out_csv=np.flip(B_Out_csv,axis=0)
        B_Out_csv= np.vstack((B_Out_csv[0], B_Out_csv[1:][::-1]))

    return B_Out_csv


def Flanged_Shaft_Output(Ratio_F,Check_F,nos_vento,divisor,D_Ext_F,t_F,Steel_F,D_Bolt_Flange_F,N_Bolt_F,Bolt_F,Bolt_Steel_F,NTRd,F_Ele_F):
    F_Out_csv=np.array([["Member","Flange Diameter","Flange thickness","Flange Steel","Bolts Inscription","Number of bolts","Bolts Diameter","Bolts steel", "Ned","Nrd","Ratio","Result"]])
    if len(Ratio_F)!=1:
        for i in range(1, nos_vento.shape[0]+1):
            inicio = (i-1) * (divisor)+1
            fim = inicio + divisor
            bloco = Ratio_F[inicio:fim]
            i_local = np.argmax(bloco)           # posição dentro do bloco
            i_global = inicio + i_local   
            if i==1:
                Name_Enc_Out="Base Flange"
            else:
                Name_Enc_Out = f"Leg - Shaft {i-1}"
            # if np.round(Ratio_F[i_global], 2)>1:
            #     Check="KO"
            # else:
            #     Check="OK"
            if Check_F[i_global]=="NA":
                row = np.array([
                    Name_Enc_Out,
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    Check_F[i_global]
                    ], dtype=object)

            else:
                row = np.array([
                    Name_Enc_Out,
                    np.round(D_Ext_F[i_global]*10**3,0),
                    np.round(t_F[i_global]*10**3,0),
                    Steel_F[i_global],
                    np.round(D_Bolt_Flange_F[i_global]*10**3,0),
                    N_Bolt_F[i_global],
                    Bolt_F[i_global],
                    Bolt_Steel_F[i_global],
                    np.round(F_Ele_F[i_global]*10**-3,1),
                    np.round(NTRd[i_global]*10**-3,1),
                    np.round(Ratio_F[i_global], 2),
                    Check_F[i_global]
                    ], dtype=object)
            F_Out_csv = np.vstack((row,F_Out_csv))
        F_Out_csv=np.flip(F_Out_csv,axis=0)
        F_Out_csv= np.vstack((F_Out_csv[0],F_Out_csv[1:][::-1]))
    
    return F_Out_csv

def Basal_Effort_Output(Ry,Ry_gelo,MzVer,MzVer_gelo,Gelo,iELU,iELU_gelo,iSLS,N_basal_ELS,N_basal_ELU_Des, N_basal_ELU_Fav,N_basal_ELS_gelo_1, N_basal_ELS_gelo_2,N_basal_ELU_gelo_1, N_basal_ELU_gelo_2,N_basal_ELU_gelo_3, N_basal_ELU_gelo_4):
    Ry = np.abs(Ry)
    MzVer = np.abs(MzVer)

    if Gelo == "Sim":
        Ry_gelo = np.abs(Ry_gelo)
        MzVer_gelo = np.abs(MzVer_gelo)

        # Máximos para ELU
        Ry_ELU = max(np.max(Ry[iELU:]), np.max(Ry_gelo[iELU_gelo:]))
        MzVer_ELU = max(np.max(MzVer[iELU:]), np.max(MzVer_gelo[iELU_gelo:]))
        # Máximos para ELS
        Ry_ELS = max(np.max(Ry[iSLS:iELU]), np.max(Ry_gelo[iSLS:iELU_gelo]))
        MzVer_ELS = max(np.max(MzVer[iSLS:iELU]), np.max(MzVer_gelo[iSLS:iELU_gelo]))
        # Esforço axial basal
        N_basal_ELU = np.max([
            N_basal_ELU_Des, N_basal_ELU_Fav,
            N_basal_ELS_gelo_1, N_basal_ELS_gelo_2,
            N_basal_ELU_gelo_1, N_basal_ELU_gelo_2,
            N_basal_ELU_gelo_3, N_basal_ELU_gelo_4
        ])
        N_basal_ELS=np.max([N_basal_ELS,N_basal_ELS_gelo_1,N_basal_ELS_gelo_2])
    else:
        Ry_ELU = np.max(Ry[iELU:])
        MzVer_ELU = np.max(MzVer[iELU:])
        Ry_ELS = np.max(Ry[iSLS:iELU])
        MzVer_ELS = np.max(MzVer[iSLS:iELU])
        N_basal_ELU = np.max([N_basal_ELU_Des, N_basal_ELU_Fav])


    Basal_Ou_csv= np.array([["Case", "Ned", "Ted", "Med"],
                    ["SLS", np.round(N_basal_ELS*10**-3,3)  ,np.round(Ry_ELS*10**-3,3) ,np.round(MzVer_ELS*10**-3,3) ],
                    ["ULS", np.round(N_basal_ELU*10**-3,3)  ,np.round(Ry_ELU*10**-3,3) ,np.round(MzVer_ELU*10**-3,3) ]
            ], dtype=object)
    return Basal_Ou_csv,MzVer_ELS

def Lattice_csv_Output(Output_Filename,Enc_Out_csv,Lig_Out_csv,B_Out_csv,F_Out_csv,Desloc_Out,Basal_Ou_csv,OrgFinal,NameLine,ValueLine,Fund_Out_Geo_Csv,Fund_Out_Ratio_Csv,Wind_Out_Csv,Ice_Out_Csv,Classe_Fiabilidade,Flange_Fund_out):
    with open(Output_Filename, "w", encoding="utf-8") as f:
        f.write("Bars" + "\n")

        f.write("\n")
        # --- primeira matriz ---
        for linha in Enc_Out_csv:
            f.write(",".join(map(str, linha)) + "\n")
        
        # --- linha em branco ---
        f.write("\n")
        
        # --- palavra chave ---
        f.write("Connections" + "\n")
        
        # --- outra linha em branco ---
        f.write("\n")
        
        # --- segunda matriz ---
        for linha in Lig_Out_csv:
            f.write(",".join(map(str, linha)) + "\n")

        f.write("\n")

        f.write("Shaft connections - Bolted " + "\n")

        f.write("\n")

        for linha in B_Out_csv:
            f.write(",".join(map(str, linha)) + "\n")

        f.write("\n")

        f.write("Shaft connections - Flanged" + "\n")

        f.write("\n")

        for linha in F_Out_csv:
            f.write(",".join(map(str, linha)) + "\n")

        
        f.write("\n")

        f.write("Tilt" + "\n")

        f.write("\n")

        for linha in Desloc_Out:
            f.write(",".join(map(str, linha)) + "\n")

        f.write("\n")

        f.write("Basal Efforts" + "\n")

        f.write("\n")

        for linha in Basal_Ou_csv:
            f.write(",".join(map(str, linha)) + "\n")

        f.write("\n")

        f.write("Area - Equipments" + "\n")

        f.write("\n")

        f.write(",".join(map(str, OrgFinal)) + "\n")

        f.write("\n")

        f.write("Area - Linear equipments" + "\n")

        f.write("\n")
        f.write("Name,Widht,Shape Coefficient,Weight,Start,End\n")
        if NameLine.size == 0:
            f.write(" ,0,0,0,0,0\n")
        else:
            ValueLine = ValueLine.astype(float)
            for i in range(ValueLine.shape[0]):
                
                f.write(f"{NameLine[i]},{np.round(ValueLine[i,0],3)},{np.round(ValueLine[i,1],1)},{np.round(ValueLine[i,3],3)},{np.round(ValueLine[i,5],1)},{np.round(ValueLine[i,6],1)}\n")

        f.write("\n")

        f.write("Foundation" + "\n")

        f.write("\n")

        for linha in Fund_Out_Geo_Csv:
            f.write(",".join(map(str, linha)) + "\n")

        f.write("\n")

        f.write("Ratio Foundation" + "\n")

        f.write("\n")

        for linha in Fund_Out_Ratio_Csv:
            f.write(",".join(map(str, linha)) + "\n")

        f.write("\n")
        f.write("Wind Load" + "\n")

        f.write("\n")

        for linha in Wind_Out_Csv:
            f.write(",".join(map(str, linha)) + "\n")

        f.write("\n")

        f.write("Ice Load" + "\n")

        f.write("\n")

        for linha in Ice_Out_Csv:
            f.write(",".join(map(str, linha)) + "\n")

        f.write("\n")

        f.write(f"Reliability Class,{int(Classe_Fiabilidade)}" + "\n")

        f.write("\n")
        f.write("Flanged Fundation" + "\n")

        f.write("\n")

        for linha in Flange_Fund_out:
            f.write(",".join(map(str, linha)) + "\n")
        f.write("\n")
def Wind_Ice_Out(Pais,z0,vref,zmin,c0_top,q0,v_top,Pressure_top,Gelo,esp_gelo,rho_gelo, Calc_Wind_Auto):
    if Calc_Wind_Auto=="Yes":
        Zona,Terreno=ZoneInverse(z0, zmin, vref, Pais)
        Wind_Out_Csv= np.array([
            ["Wind Load",""],
            ["Country",Pais],
            ["vb (m/s)" , f"{Zona} - {int(vref)}"],
            ["z0", z0],
            ["zmin", zmin],
            ["Orographie Coefficient", c0_top],
            ["q0 (Pa)",  np.round(q0,1) ],
            ["Average Wind speed (top) (m/s)", np.round(v_top,2)],
            ["Dynamic Wind Pressure (top) (Pa)",np.round(Pressure_top)]
        ], dtype=object)
    else:
        if Pais in ["Portugal","Spain","France"]:
            Zona,Terreno=ZoneInverse(z0, zmin, vref, Pais)
            Wind_Out_Csv= np.array([
            ["Wind Load",""],
            ["Country",Pais],
            ["Wind Zone", Zona],
            ["vb (m/s)" , vref],
            ["Terrain Category", Terreno],
            ["Orographie Coefficient", c0_top],
            ["q0 (Pa)",  np.round(q0,1) ],
            ["Average Wind speed (top) (m/s)", np.round(v_top,2)],
            ["Dynamic Wind Pressure (top) (Pa)",np.round(Pressure_top)]
        ], dtype=object)
        elif Pais in ["Portugal-RSA"]:
            Zona,Terreno=ZoneInverse(z0, zmin, vref, Pais)
            Wind_Out_Csv= np.array([
            ["Wind Load",""],
            ["Country",Pais],
            ["Wind Zone", Zona],
            ["vb (m/s)" , vref],
            ["Terrain Category", Terreno],
            ["Orographie Coefficient", c0_top],
            ["q0 (Pa)",  np.round(q0,1) ],
            ["Average Wind speed (top) (m/s)", np.round(v_top,2)],
            ["Dynamic Wind Pressure (top) (Pa)",np.round(Pressure_top)]
            ], dtype=object)

    if Gelo=="Sim":
        Ice_Out_Csv=np.array([
        ["Snow Load",""],
        ["Radial ice thickness",esp_gelo],
        ["Increased volumetric weight of ice (kg/m3)", rho_gelo]
    ], dtype=object)
    else:
        Ice_Out_Csv=np.array([
        ["Snow Load",""],
        ["Radial ice thickness",0],
        ["Increased volumetric weight of ice (kg/m3)", 0]
    ], dtype=object)

    return Wind_Out_Csv,Ice_Out_Csv        