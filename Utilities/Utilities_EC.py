from Utilities.Utilities import integral,find_closest_value_position,loaddistribution
import numpy as np
def ClasseFiabilidade(Pais,Classe_Fiabilidade):
    if Pais in ["Geral","Portugal","France","Spain"]:
        if Classe_Fiabilidade==3:
            Gamma_Perm_Des=1.2
            Gamma_Variable_Des=1.6
            Gamma_Perm_Fav=1
            Gamma_Variable_Fav=0
        elif Classe_Fiabilidade==2:
            Gamma_Perm_Des=1.1
            Gamma_Variable_Des=1.4
            Gamma_Perm_Fav=1
            Gamma_Variable_Fav=0
        elif Classe_Fiabilidade==1:
            Gamma_Perm_Des=1
            Gamma_Variable_Des=1.2
            Gamma_Perm_Fav=1
            Gamma_Variable_Fav=0
    elif Pais in ["Portugal-RSA"]:
        Gamma_Perm_Des=1.35
        Gamma_Variable_Des=1.5
        Gamma_Perm_Fav=1
        Gamma_Variable_Fav=0
    elif Pais in ["Italy"]:
        Gamma_Perm_Des=1.35
        Gamma_Variable_Des=1.5
        Gamma_Perm_Fav=1
        Gamma_Variable_Fav=0        
    
    return  Gamma_Perm_Des,Gamma_Variable_Des,Gamma_Perm_Fav,Gamma_Variable_Fav    

def Equivalent_Mass(h_vetor,Posicao_Antenna,Ant,h_torre,Massa_Antenna,equipline,dh,z_vento,Massa_Troco,h_equip,Gelo):
    m_vector=np.zeros(len(h_vetor))
    equip_mass_loc=np.zeros(len(h_vetor))
    equip_mass_dist=np.zeros(len(h_vetor))
    for i in range(len(h_vetor)):
        for j in range(len(z_vento)-1):
            if h_vetor[i]<z_vento[j] and h_vetor[i]>=z_vento[j+1]:
                m_vector[i]=Massa_Troco[j+1]/(z_vento[j]-z_vento[j+1])*dh
    for i in range(len(Posicao_Antenna)):
        if Ant[i,0]<=h_torre:
            j=find_closest_value_position(h_vetor,Posicao_Antenna[i])
            if h_vetor[j]<Ant[i,0]:
                j1=j
                j2=j+1
            else:
                j1=j-1
                j2=j
            aux1,aux2=loaddistribution(Posicao_Antenna[i],h_vetor[j1],h_vetor[j2],Massa_Antenna[i])
            equip_mass_loc[j1]=equip_mass_loc[j1]+aux1
            equip_mass_loc[j2]=equip_mass_loc[j2]+aux2
            #aux1,aux2=loaddistribution(Posicao_Antenna[i],h_vetor[j1],h_vetor[j2],Ant[i,2])
            #equip_area_loc[j1]=equip_area_loc[j1]+aux1
            #equip_area_loc[j2]=equip_area_loc[j2]+aux2
            #Aset=Aset+Wind_Pressure[j1]*aux1*h_vetor[j1]/(Wind_Pressure[-1]*h_torre)+Wind_Pressure[j2]*aux2*h_vetor[j2]/(Wind_Pressure[-1]*h_torre)
        else:
            equip_mass_loc[-1]=Ant[i,1]+equip_mass_loc[-1]
            j=find_closest_value_position(h_equip,Ant[i,0])
            #equip_area_loc[-1]=equip_area_loc[-1]+Ant[i,2]*Wind_Pressure_equip[j]/Wind_Pressure[-1]  #Aumento da area associada ao aumento de pressao.
            ##Não tem em consideração o aumento de momento causado pela carga, apenas a força
            #Aset=Aset+Wind_Pressure_equip[j]*Ant[i,2]*h_equip[j]/(Wind_Pressure[-1]*h_torre)
    for i in range(equipline.shape[0]):
        idist=(h_vetor>=equipline[i,-2]) & (h_vetor<=equipline[i,-1])
        if Gelo=="Sim":
            equip_mass_dist[idist]=equip_mass_dist[idist]+equipline[i,4]*dh
        else:
            equip_mass_dist[idist]=equip_mass_dist[idist]+equipline[i,3]*dh
    phi_me=(h_vetor/h_torre)**2.5 #EC vento F.3
    phi2_me=phi_me*phi_me
    m_vector=m_vector+equip_mass_loc+equip_mass_dist

    me=integral(h_vetor,phi2_me*m_vector/dh)/integral(h_vetor,phi2_me)

    return me,equip_mass_loc,equip_mass_dist,m_vector

def Combined_Cases(vetores, constantes):
    if len(vetores) != len(constantes):
        raise ValueError("Número de vetores e constantes deve ser igual.")
    
    # Garante que tudo é array 2D (remove dimensões de tamanho 1)
    vetores = [np.squeeze(np.array(v, dtype=float)) for v in vetores]

    # Verifica se todos têm a mesma forma
    shapes = [v.shape for v in vetores]
    if len(set(shapes)) != 1:
        raise ValueError(f"As formas dos vetores não coincidem: {shapes}")
    
    resultado = np.zeros_like(vetores[0], dtype=float)
    for v, c in zip(vetores, constantes):
        resultado += c * v
    return resultado