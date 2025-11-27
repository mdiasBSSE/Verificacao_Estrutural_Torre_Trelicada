from Utilities.Utilities_EC import Combined_Cases
import numpy as np
def ELS(Desloc,Ry,MzVer,N_basal,F,Reactions,F_axial,iv,icomb):

    Desloc[icomb,:,:]=Combined_Cases([Desloc[1,:,:],Desloc[2,:,:],[Desloc[3,:,:]],[Desloc[iv,:,:]]],[1.0,1.0,1.0,1.0])
    F[icomb,:,:]=Combined_Cases([F[1,:,:],F[2,:,:],[F[3,:,:]],[F[iv,:,:]]],[1.0,1.0,1.0,1.0])
    F_axial[icomb,:]=Combined_Cases([F_axial[1,:],F_axial[2,:],[F_axial[3,:]],[F_axial[iv,:]]],[1.0,1.0,1.0,1.0])
    Ry[icomb]=1.0*Ry[iv]
    MzVer[icomb]=1.0*MzVer[iv]
    N_basal_ELS=1.0*N_basal
    Reactions[icomb]=Combined_Cases([Reactions[1,:],Reactions[2,:],[Reactions[3,:]],[Reactions[iv,:]]],[1.0,1.0,1.0,1.0])

    return Desloc,F,F_axial,Reactions,Ry,MzVer,N_basal_ELS


def ELU(Desloc,Ry,MzVer,N_basal,F,Reactions,F_axial,iv,icomb,Gamma_Perm_Des,Gamma_Variable_Des):

    Desloc[icomb,:,:]=Combined_Cases([Desloc[1,:,:],Desloc[2,:,:],Desloc[3,:,:],Desloc[iv,:,:]],[Gamma_Perm_Des,Gamma_Perm_Des,Gamma_Perm_Des,Gamma_Variable_Des])
    F[icomb,:,:]=Combined_Cases([F[1,:,:],F[2,:,:],[F[3,:,:]],[F[iv,:,:]]],[Gamma_Perm_Des,Gamma_Perm_Des,Gamma_Perm_Des,Gamma_Variable_Des])
    F_axial[icomb,:]=Combined_Cases([F_axial[1,:],F_axial[2,:],[F_axial[3,:]],[F_axial[iv,:]]],[Gamma_Perm_Des,Gamma_Perm_Des,Gamma_Perm_Des,Gamma_Variable_Des])
    Ry[icomb]=Gamma_Variable_Des*Ry[iv]
    MzVer[icomb]=Gamma_Variable_Des*MzVer[iv]
    N_basal_ELU=Gamma_Perm_Des*N_basal
    Reactions[icomb]=Combined_Cases([Reactions[1,:],Reactions[2,:],[Reactions[3,:]],[Reactions[iv,:]]],[Gamma_Perm_Des,Gamma_Perm_Des,Gamma_Perm_Des,Gamma_Variable_Des])

    return Desloc,F,F_axial,Reactions,Ry,MzVer,N_basal_ELU


def Comb_Gelo(Desloc,Ry,MzVer,N_basal,F,F_nos,Reactions,F_axial,Desloc_gelo,Ry_gelo,MzVer_gelo,N_basal_gelo,F_gelo,Reactions_gelo,F_axial_gelo,iv,icomb,Factor_Weight,Factor_Wind,Factor_Snow):

    Desloc_gelo[icomb,:,:]=Combined_Cases([Desloc[1,:,:],Desloc[2,:,:],Desloc[3,:,:],Desloc_gelo[iv,:,:],Desloc_gelo[1,:,:],Desloc_gelo[2,:,:],Desloc_gelo[3,:,:]],[Factor_Weight,Factor_Weight,Factor_Weight,Factor_Wind,Factor_Snow,Factor_Snow,Factor_Snow])
    F_gelo[icomb,:,:]=Combined_Cases([F[1,:,:],F[2,:,:],F[3,:,:],F_gelo[iv,:,:],F_gelo[1,:,:],F_gelo[2,:,:],F_gelo[3,:,:]],[Factor_Weight,Factor_Weight,Factor_Weight,Factor_Wind,Factor_Snow,Factor_Snow,Factor_Snow])
    F_axial_gelo[icomb,:]=Combined_Cases([F_axial[1,:],F_axial[2,:],F_axial[3,:],F_axial_gelo[iv,:],F_axial_gelo[1,:],F_axial_gelo[2,:],F_axial_gelo[3,:]],[Factor_Weight,Factor_Weight,Factor_Weight,Factor_Wind,Factor_Snow,Factor_Snow,Factor_Snow])
    Ry_gelo[icomb]=Factor_Wind*Ry_gelo[iv]
    MzVer_gelo[icomb]=Factor_Wind*MzVer_gelo[iv]
    N_basal_ELS_gelo_1=Factor_Weight*N_basal+Factor_Snow*N_basal_gelo
    Reactions_gelo[icomb,:,:]=Combined_Cases([Reactions[1,:,:],Reactions[2,:,:],Reactions[3,:,:],Reactions_gelo[iv,:,:],Reactions_gelo[1,:,:],Reactions_gelo[2,:,:],Reactions_gelo[3,:,:]],[Factor_Weight,Factor_Weight,Factor_Weight,Factor_Wind,Factor_Snow,Factor_Snow,Factor_Snow])

    return Desloc_gelo,F_gelo,F_axial_gelo,Reactions_gelo,Ry_gelo,MzVer_gelo,N_basal_ELS_gelo_1