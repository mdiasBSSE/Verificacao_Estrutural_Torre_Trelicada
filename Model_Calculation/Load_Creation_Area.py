import numpy as np
from Equipments.Calc_Areas import Antenatxt
from Inputs.Lattice_Tower_Input import Ant_txt_read
from openseespy import opensees as ops
from Utilities.Utilities import encontra_vizinhos,agrupar_indices
from Utilities.Utilities_ops import distancia_entre_nos,comprimentos_elementos_array,no_mais_proximo,elementos_do_no,elementos_do_no_filt,elementos_entre,elementos_acima,elementos_abaixo,get_nodes_same_elevation
from Utilities.Utilities_EC import ClasseFiabilidade
def Area_Definition(ant_txt_file,h_vetor,h_torre,Gelo,esp_gelo,rho_gelo,dalpha,alpha_norte,Aincrease):
    Antenas,equipline_manual,Auto,Ladder,CabosOutside,Final,Nanalise,Org=Ant_txt_read(ant_txt_file)
    # equipline=np.array([[0.125,1,0,7.5,0,0,30],[0.1,1,0,5,0,0,30]])
    # #equipline=np.array([[0.1225,1.2,0,7.5,0,0,24],[0.2,1.2,0,10,0,0,24]])
    # #equipline=np.array([[0.18,1,0,25.49,0,0,30.4]])
    # #equipline=np.array([[0.18,1,0,25.49,0,0,30.4]])
    # #equipline=np.array([[0.1,1,0.15,6,8,0,32],[0.2,1,0.3,8,10,0,32]])
    # #equipline=np.array([[0.1,1.2,0,15,0,0,20]])
    # #equipline=np.array([])
    if Aincrease!=0:
        Antenas.append(f"8;0;{h_torre};{Aincrease};{Aincrease*150}")

    Ant, Nant, AreaAnt,AreaDirection,AreaDirectiongelo = Antenatxt(Antenas,h_vetor,h_torre,Gelo,esp_gelo,rho_gelo,dalpha,alpha_norte)
    
    NameLine=[]

    NameLine=np.array(NameLine,dtype='<U64')
    MomentArea=0
    MomentoPeso=0
    if Auto==1:      
        Line_Equipment=np.empty((1,7))
        poscabo=[]
        if Ladder==1:
            Line_Equipment=np.vstack((Line_Equipment,[0.1225,1.2,0.1225-16*esp_gelo**2+7.0296*esp_gelo,8,8+(8.8168*esp_gelo**2+0.5248*esp_gelo)*rho_gelo,0,h_torre]))
            NameLine = np.append(NameLine,"Ladder") 
        for i in range(len(Nant)-1):
            if Nant[i]!=Nant[i+1]:
                poscabo.append([Nant[i],h_vetor[i]])
        poscabo.append([Nant[-1],h_vetor[-1]])
        
        
        ValueLine=Line_Equipment
        poscabo=np.array(poscabo)
        for i in range(poscabo.shape[0]):
            if poscabo[i,1]>h_torre:
                poscabo[i,1]=h_torre
            if CabosOutside==1:
                if i==0:
                    Line_Equipment=np.vstack((Line_Equipment,[poscabo[i,0]*0.03,1.2,poscabo[i,0]*(0.03+2*esp_gelo),poscabo[i,0]*0.5,poscabo[i,0]*(0.5+(0.7854*esp_gelo**2+0.0188*esp_gelo)*rho_gelo),0,poscabo[i,1]])) 
                else:
                    Line_Equipment=np.vstack((Line_Equipment,[poscabo[i,0]*0.03,1.2,poscabo[i,0]*(0.03+2*esp_gelo),poscabo[i,0]*0.5,poscabo[i,0]*(0.5+(0.7854*esp_gelo**2+0.0188*esp_gelo)*rho_gelo),poscabo[i-1,1],poscabo[i,1]]))
            else:
                if i==0:
                    Line_Equipment=np.vstack((Line_Equipment,[0,1.2,0,poscabo[i,0]*0.5,poscabo[i,0]*(0.5+(0.7854*esp_gelo**2+0.0188*esp_gelo)*rho_gelo),0,poscabo[i,1]])) 
                else:
                    Line_Equipment=np.vstack((Line_Equipment,[0,1.2,0,poscabo[i,0]*0.5,poscabo[i,0]*(0.5+(0.7854*esp_gelo**2+0.0188*esp_gelo)*rho_gelo),poscabo[i-1,1],poscabo[i,1]]))
            if i==0:
                DeltaH=poscabo[i,1]-0
            else:
                DeltaH=poscabo[i,1]-poscabo[i-1,1]
            MomentArea=MomentArea+Line_Equipment[-1,0]*(DeltaH)*(Line_Equipment[-1,-2]+0.6*DeltaH)
            MomentoPeso=MomentoPeso+Line_Equipment[-1,3]*(DeltaH)*(Line_Equipment[-1,-2]+0.6*DeltaH)
        Line_Equipment=Line_Equipment.astype(float)
        ValueLine=np.vstack((ValueLine,[MomentArea/(0.6*h_torre**2),1.2,0,MomentoPeso/(0.6*h_torre**2),0,0,h_torre]))
        NameLine = np.append(NameLine, "Cables") 
        ValueLine=np.delete(ValueLine,0,axis=0)
        Line_Equipment=np.delete(Line_Equipment,0,axis=0)

    else:
        Line_Equipment=np.empty((len(equipline_manual),7),dtype='<U64')
        for i in range(len(Line_Equipment)):
            Aux_equip=(equipline_manual[i].strip()).split(';')
            Line_Equipment[i,0]=float(Aux_equip[3])
            Line_Equipment[i,1]=1.2
            Line_Equipment[i,2]=float(Aux_equip[3])+2*esp_gelo
            Line_Equipment[i,3]=float(Aux_equip[4])
            Line_Equipment[i,4]=float(Aux_equip[4])+(0.7854*esp_gelo**2+0.0188*esp_gelo)*rho_gelo
            Line_Equipment[i,5]=float(Aux_equip[1])
            Line_Equipment[i,6]=float(Aux_equip[2])
            NameLine = np.append(NameLine,Aux_equip[0])  
        if Ladder==1:
            Flag_escada=0
            for i in range(len(equipline_manual)):
                aux_esc=np.array((equipline_manual[i].strip()).split(';'))
                if aux_esc[0]=="Ladder":
                    Flag_escada=1 
            if Flag_escada==0:
                Line_Equipment=np.vstack((Line_Equipment,[0.1225,1.2,0.1225-16*esp_gelo**2+7.0296*esp_gelo,8,8+(8.8168*esp_gelo**2+0.5248*esp_gelo)*rho_gelo,0,h_torre]))
                NameLine = np.append(NameLine,"Ladder") 
        if CabosOutside==0:
            for i in range(len(equipline_manual)):
                aux_ladder=np.array((equipline_manual[i].strip()).split(';'))
                if aux_ladder[0]=="Cable":
                    Line_Equipment[i,0]=0
                    Line_Equipment[i,2]=0
        ValueLine=Line_Equipment 
        Line_Equipment=np.array(Line_Equipment,dtype=float)

    equipline=Line_Equipment



    Massa_Antenna=Ant[:,1]
    Massa_Antenna_Gelo=Ant[:,3]
    Posicao_Antenna=Ant[:,0]
    return Ant, Nant, AreaAnt,AreaDirection,AreaDirectiongelo,equipline,Line_Equipment,Auto,Ladder,CabosOutside,Final,Nanalise,Org,ValueLine,NameLine,Massa_Antenna,Massa_Antenna_Gelo,Posicao_Antenna,Antenas

def Area_Position_Mass(Posicao_Antenna,h_torre,z_grupo,nos_grupo,Massa_Antenna,divisor,equipline,TrussType,todos_nos,Gelo):
    ##Posicionamento da antena e distribuição da massa nos nós
    iAnt=np.zeros((len(Posicao_Antenna),2))
    x_Ant=np.zeros(len(Posicao_Antenna))
    Massa_nos=np.zeros(len(todos_nos)+2)
    Massa_nos_lin=np.zeros(len(todos_nos)+2)

    for i in range(len(Posicao_Antenna)):
        if Posicao_Antenna[i]<=h_torre:
            iAnt[i,:]=encontra_vizinhos(z_grupo[::-1], Posicao_Antenna[i])
            h1_Ant=z_grupo[int(iAnt[i,1])]
            h2_Ant=z_grupo[int(iAnt[i,0])]
            x_Ant[i]=(Posicao_Antenna[i]-h1_Ant)/(h2_Ant-h1_Ant)
        else:
            iAnt[i,0]=0
            iAnt[i,1]=0
            
            h1_Ant=z_grupo[int(iAnt[i,1])]
            h2_Ant=z_grupo[int(iAnt[i,0])]
            x_Ant[i]=0
        for j in range(nos_grupo.shape[1]):
            Massa_nos[int(nos_grupo[int(iAnt[i,1]),j])]=Massa_nos[int(nos_grupo[int(iAnt[i,1]),j])]+(1-x_Ant[i])*Massa_Antenna[i]/divisor
            Massa_nos[int(nos_grupo[int(iAnt[i,0]),j])]=Massa_nos[int(nos_grupo[int(iAnt[i,0]),j])]+x_Ant[i]*Massa_Antenna[i]/divisor
    if Gelo=="Sim":
        iP=int(4)
    else:
        iP=int(3)
    for i in range(equipline.shape[0]):
        elem_completo=elementos_entre(equipline[i,-2],equipline[i,-1])
        res = elementos_acima(equipline[i,-2])

        if res.size>0:  # se não estiver vazia
            elem_acima, Per_cima = zip(*res)  # desempacota elementos e percentagens
            elem_acima = np.array(list(elem_acima))
            Per_cima = np.array(list(Per_cima))
        else:
            elem_acima =  np.empty(0)
            Per_cima = np.empty(0)
        res = elementos_abaixo(equipline[i,-1])

        if res.size>0:  # se a lista não estiver vazia
            elem_abaixo, Per_baixo = zip(*res)  # desempacota elementos e percentagens
            elem_abaixo = np.array(list(elem_abaixo))
            Per_baixo = np.array(list(Per_baixo))
        else:
            elem_abaixo = np.empty(0)
            Per_baixo = np.empty(0)
        if elem_completo.size>0:  # só entra se a lista não estiver vazia
            for j in range(len(elem_completo)):
                elem = elem_completo[j]
                
                if TrussType[elem] == "Leg":
                    ni, nj = ops.eleNodes(int(elem))

                    # obter coordenadas z dos nós
                    zi = ops.nodeCoord(ni)[2]
                    zj = ops.nodeCoord(nj)[2]
                    dz = abs(zj - zi)

                    # distribuir massa nos nós (acumular, não substituir!)
                    Massa_nos_lin[ni] =  Massa_nos_lin[ni]  + equipline[i, iP] * dz / 2/divisor
                    Massa_nos_lin[nj] = Massa_nos_lin[nj] + equipline[i, iP] * dz / 2/divisor
        
        if elem_abaixo.size>0:
            for j in range(len(elem_abaixo)):
                elem = elem_abaixo[j]
                
                if TrussType[int(elem)] == "Leg":
                    n1, n2 = ops.eleNodes(int(elem))
                    # obter coordenadas z
                    z1 = ops.nodeCoord(n1)[2]
                    z2 = ops.nodeCoord(n2)[2]
                    # obter coordenadas z dos nós
                    zi = ops.nodeCoord(ni)[2]
                    zj = ops.nodeCoord(nj)[2]
                    # garantir que ni tem z menor e nj z maior
                    if z1 <= z2:
                        ni, nj = n1, n2
                        zi, zj = z1, z2
                    else:
                        ni, nj = n2, n1
                        zi, zj = z2, z1
                    dz = abs(zj - zi)
                    zc=(equipline[i,-1]-zi)/2

                    # distribuir massa nos nós (acumular, não substituir!)
                    Massa_nos_lin[int(ni)] = Massa_nos_lin[int(ni)] + equipline[i, iP] * Per_baixo[j]*(zj-zc-zi)/divisor
                    Massa_nos_lin[int(nj)] = Massa_nos_lin[int(nj)] + equipline[i, iP] * Per_baixo[j]*(zc-zi+zi)/divisor
        if elem_acima.size>0:
            for j in range(len(elem_acima)):
                elem = elem_acima[j]
                
                if TrussType[int(elem)] == "Leg":
                    n1, n2 = ops.eleNodes(int(elem))
                    # obter coordenadas z
                    z1 = ops.nodeCoord(n1)[2]
                    z2 = ops.nodeCoord(n2)[2]
                    # obter coordenadas z dos nós
                    zi = ops.nodeCoord(ni)[2]
                    zj = ops.nodeCoord(nj)[2]
                    # garantir que ni tem z menor e nj z maior
                    if z1 <= z2:
                        ni, nj = n1, n2
                        zi, zj = z1, z2
                    else:
                        ni, nj = n2, n1
                        zi, zj = z2, z1
                    dz = abs(zj - zi)
                    zc=(zj-equipline[i,-2])/2


                    # distribuir massa nos nós (acumular, não substituir!)
                    Massa_nos_lin[int(nj)] =  Massa_nos_lin[int(nj)]  + equipline[i, iP] * Per_cima[j] * (zj-zc-zi)/divisor
                    Massa_nos_lin[int(ni)] = Massa_nos_lin[int(ni)] + equipline[i, iP] * Per_cima[j]*(zc-zi+zi)/divisor
    return iAnt,x_Ant, Massa_nos,Massa_nos_lin