import numpy as np
from Utilities.Utilities import find_closest_value_position

def Antenatxt(Antenas,h_vector,h_torre,neve,esp_gelo,rho_gelo,dalpha,alpha_norte):
    Ant=np.zeros((len(Antenas),6))
    Antout=np.empty((len(Antenas),1),dtype='U64')
    #Used=np.zeros(len(Antenas))
    Ocupation=np.zeros(len(h_vector))
    
    Cf_frente=0.95
    Cf_lado=0.85
    Cf_tras=1.1
    Cf_RRU=1.4
    Cf_MW=1.4
    Nant=np.zeros(len(h_vector))
    RRU=np.zeros(len(Antenas))
    alpha_norte=alpha_norte*np.pi/180
    Alpha_wind=np.arange(0,360,dalpha)
    Alpha_wind=Alpha_wind*np.pi/180
    AreaDirection=np.zeros((len(Antenas),len(Alpha_wind)))
    AreaDirectiongelo=np.zeros((len(Antenas),len(Alpha_wind)))
    for i in range(len(Antenas)):
        AntenasAux=np.array((Antenas[i].strip()).split(';'))
        EquipType=int(AntenasAux[0])
        ConfigType=int(AntenasAux[1])
        Height=float(AntenasAux[2])
        ih=find_closest_value_position(h_vector,Height)
        Ant[i,0]=Height #Posição
        if EquipType!=8:
            azimute=float(AntenasAux[3])*np.pi/180
            dalpha=Alpha_wind-(azimute - alpha_norte)
            dalpha[np.where(dalpha>2*np.pi)]=dalpha[np.where(dalpha>2*np.pi)]-2*np.pi
            dalpha[np.where(dalpha<0)]=dalpha[np.where(dalpha<0)]+2*np.pi
        #Calculo das areas
        Af=0
        Al=0
        At=0
        Peso=0
        mgelo=0
        Afgelo=0
        Algelo=0
        Atgelo=0
        if EquipType==1: #Antenna from a Catalog
            if ConfigType<4:
                for j in range(ConfigType):
                    Af=Af+float(AntenasAux[8*j+10])
                    Al=Al+float(AntenasAux[8*j+11])
                    At=At+float(AntenasAux[8*j+12])
                    L1=float(AntenasAux[8*j+7])*10**-3
                    L2=float(AntenasAux[8*j+8])*10**-3
                    L3=float(AntenasAux[8*j+9])*10**-3
                    Peso=Peso+float(AntenasAux[8*j+13])
                    indicehmax=find_closest_value_position(h_vector,Height+L1/2)
                    indicehmin=find_closest_value_position(h_vector,Height-L1/2)
                    Ocupation[indicehmin:indicehmax+1]=np.ones(indicehmax+1-indicehmin)
                    if neve=="Sim":
                        Vant=L1*L2*L3
                        L1gelo=L1+2*esp_gelo
                        L2gelo=L2+2*esp_gelo
                        L3gelo=L3+2*esp_gelo
                        Vantgelo=(L1gelo)*(L2gelo)*(L3gelo)
                        Vgelo=Vantgelo-Vant
                        mgelo=mgelo+Vgelo*rho_gelo+float(AntenasAux[8*j+13])
                        Afront=L1*L2
                        Alateral=L1*L3
                        Atraseira=L1*L2
                        Afgelo=Afgelo+L1gelo*L2gelo*float(AntenasAux[8*j+10])/Afront #Assumindo que o Cf não muda quando temos gelo
                        Algelo=Algelo+L1gelo*L3gelo*float(AntenasAux[8*j+11])/Alateral
                        Atgelo=Atgelo+L1gelo*L2gelo*float(AntenasAux[8*j+12])/Atraseira
                Nant[0:ih+1]=Nant[0:ih+1]+np.ones(ih+1)*ConfigType
        elif EquipType==2: #RRU from a Catalog            
            j=np.arange(ConfigType)
            Af=max(AntenasAux[8*j+10].astype(float))
            Al=sum(AntenasAux[8*j+11].astype(float))
            At=max(AntenasAux[8*j+12].astype(float))
            Peso=sum(AntenasAux[8*j+13].astype(float)) #Massa
            if neve=="Sim":
                L1=AntenasAux[8*j+7].astype(float)*10**-3
                L2=AntenasAux[8*j+8].astype(float)*10**-3
                L3=AntenasAux[8*j+9].astype(float)*10**-3
                L1gelo=L1+2*esp_gelo
                L2gelo=L2+2*esp_gelo
                L3gelo=L3+2*esp_gelo
                Vant=L1*L2*L3
                Vantgelo=L1gelo*L2gelo*L3gelo
                Vgelo=Vantgelo-Vant
                mgelo=sum(Vgelo*rho_gelo)+Peso
                Afront=L1*L2
                Alateral=L2*L3
                Atraseira=L1*L2
                Afgelo=max(L1gelo*L2gelo*AntenasAux[8*j+10].astype(float)/Afront)
                Algelo=sum(L1gelo*L3gelo*AntenasAux[8*j+11].astype(float)/Alateral)
                Atgelo=max(L1gelo*L2gelo*AntenasAux[8*j+12].astype(float)/Atraseira) 
            RRU[i]=1
        elif EquipType==3: #Microwave form catalog
            if ConfigType==0:
                Af=AntenasAux[9].astype(float)
                Al=AntenasAux[10].astype(float)
                At=AntenasAux[11].astype(float)
                Peso=AntenasAux[12].astype(float)
                if neve=="Sim":
                    D1=AntenasAux[7].astype(float)*10**-3
                    L1=AntenasAux[8].astype(float)*10**-3
                    D1gelo=D1+2*esp_gelo
                    L1gelo=L1+2*esp_gelo
                    Vant=(np.pi*D1**2/4)*L1
                    Vantgelo=(np.pi*D1gelo**2/4)*L1gelo
                    Vgelo=Vantgelo-Vant
                    mgelo=Peso+Vgelo*rho_gelo
                    Afront=np.pi*D1**2/4
                    Alateral=D1*L1
                    Atraseira=np.pi*D1**2/4
                    Afgelo=(np.pi*D1gelo**2/4)*Af/Afront
                    Algelo=D1gelo*L1gelo*Al/Alateral
                    Atgelo=(np.pi*D1gelo**2/4)*At/Atraseira
        elif EquipType==4: #Other Equipment
            if ConfigType<4:
                j=np.arange(ConfigType)
                L1=AntenasAux[5*j+7].astype(float)*10**-3
                L2=AntenasAux[5*j+8].astype(float)*10**-3
                L3=AntenasAux[5*j+9].astype(float)*10**-3
                Af=sum(L1*L2)*Cf_frente
                Al=sum(L1*L3)*Cf_lado
                At=sum(L1*L2)*Cf_tras
                Peso=sum(AntenasAux[5*j+10].astype(float))
                if neve=="Sim":
                    L1gelo=L1+2*esp_gelo
                    L2gelo=L2+2*esp_gelo
                    L3gelo=L3+2*esp_gelo
                    Vant=L1*L2*L3
                    Vantgelo=L1gelo*L2gelo*L3gelo


                    Vgelo=Vantgelo-Vant
                    mgelo=sum(Vgelo*rho_gelo)+Peso
                    Afgelo=sum(L1gelo*L2gelo)*Cf_frente
                    Algelo=sum(L1gelo*L3gelo)*Cf_lado
                    Atgelo=sum(L1gelo*L3gelo)*Cf_tras
                for jaux in j:
                    L1aux=float(AntenasAux[5*jaux+7])
                    indicehmax=find_closest_value_position(h_vector,Height+L1aux/2)
                    indicehmin=find_closest_value_position(h_vector,Height-L1aux/2)
                    Ocupation[indicehmin:indicehmax+1]=np.ones(indicehmax+1-indicehmin)
                Nant[0:ih+1]=Nant[0:ih+1]+np.ones(ih+1)*ConfigType
        elif EquipType==5:
            if ConfigType<4:
                j=np.arange(ConfigType)
                L1=AntenasAux[5*j+7].astype(float)*10**-3
                L2=AntenasAux[5*j+8].astype(float)*10**-3
                L3=AntenasAux[5*j+9].astype(float)*10**-3
                Af=max(L1*L2)*Cf_RRU
                Al=sum(L1*L3)*Cf_RRU
                At=max(L1*L2)*Cf_RRU
                Peso=sum(AntenasAux[5*j+10].astype(float))
                if neve=="Sim":
                    L1gelo=L1+2*esp_gelo
                    L2gelo=L2+2*esp_gelo
                    L3gelo=L3+2*esp_gelo
                    Vant=L1*L2*L3
                    Vantgelo=L1gelo*L2gelo*L3gelo
                    Vgelo=Vantgelo-Vant
                    mgelo=sum(Vgelo*rho_gelo)+Peso
                    Afgelo=max(L1gelo*L2gelo)*Cf_RRU
                    Algelo=sum(L1gelo*L3gelo)*Cf_RRU
                    Atgelo=max(L1gelo*L3gelo)*Cf_RRU            
                RRU[i]=1
        elif EquipType==6: #Other Microwave
            if ConfigType==0:
                D1=AntenasAux[7].astype(float)*10**-3
                L1=AntenasAux[8].astype(float)*10**-3
                Af=np.pi*(D1)**2*Cf_MW/4
                Al=D1*L1*Cf_MW
                At=Af
                Peso=AntenasAux[9].astype(float)
                if neve=="Sim":
                    D1gelo=D1+2*esp_gelo
                    L1gelo=L1+2*esp_gelo
                    Afgelo=np.pi*D1gelo**2*Cf_MW/4
                    Algelo=D1gelo*L1gelo*Cf_MW
                    Atgelo=Afgelo
                    Vant=np.pi*D1**2*L1
                    Vantgelo=np.pi*D1gelo**2*L1gelo
                    Vgelo=Vantgelo-Vant
                    mgelo=Peso+Vgelo*rho_gelo
        elif EquipType==7: #Working Plataform
            if ConfigType==1:
                Coef_Plat=0.373 #Circular
                if neve=="Sim":
                    Coef_Plat_gelo=0.0147*esp_gelo*10**3+0.373
                    Coef_gelo_peso=4.152*10**-5*(esp_gelo*10**3)**2+2.193*10**-3*esp_gelo*10**3
            elif ConfigType==2:
                Coef_Plat=0.466 #Triangular
                
                if neve=="Sim":
                    Coef_Plat_gelo=0.0147*esp_gelo*10**3+0.466
                    Coef_gelo_peso=3.965*10**-5*(esp_gelo*10**3)**2+2.094*10**-3*esp_gelo*10**3
                    
            elif ConfigType==3:
                Coef_Plat=0.659 #Quadrangular
                if neve=="Sim":
                    Coef_Plat_gelo=0.0147*esp_gelo*10**3+0.659
                    Coef_gelo_peso=5.287*10**-5*(esp_gelo*10**3)**2+2.792*10**-3*esp_gelo*10**3
            Af=AntenasAux[4].astype(float)*10**-3*Coef_Plat
            Al=Af
            At=Af
            Peso=AntenasAux[5].astype(float)
            if neve=="Sim":
                Afgelo=AntenasAux[4].astype(float)*10**-3*Coef_Plat_gelo
                Algelo=Afgelo
                Atgelo=Afgelo
                mgelo=Peso+Coef_gelo_peso*rho_gelo

        elif EquipType==8:
            if ConfigType==0:
                Af=AntenasAux[3].astype(float)
                Al=Af
                At=Af
                Peso=AntenasAux[4].astype(float)
                if neve=="Sim":
                    Afgelo=(1+19.506*esp_gelo**2+11.55*esp_gelo)*Af #Curva criada considerando o aumento causado por algumas antenas 
                    Algelo=Afgelo
                    Atgelo=Atgelo
                    Vgeloporagelo=29.297*esp_gelo**3+22.255*esp_gelo**2+4.9875*esp_gelo
                    Vgelo=Vgeloporagelo*Afgelo
                    mgelo=Peso+Vgelo*rho_gelo
        elif EquipType==9:
            if ConfigType==1:
                Coef_Patim=0.176 #Circular
                if neve=="Sim":
                    Coef_Patim_gelo=0.004*esp_gelo*10**3+0.176
                    Coef_gelo_peso=4.152*10**-5*(esp_gelo*10**3)**2+2.193*10**-3*esp_gelo*10**3
            elif ConfigType==2:
                Coef_Patim=0.220 #Triangular
                if neve=="Sim":
                    Coef_Patim_gelo=0.004*esp_gelo*10**3+0.220
                    Coef_gelo_peso=3.965*10**-5*(esp_gelo*10**3)**2+2.094*10**-3*esp_gelo*10**3
            elif ConfigType==3:
                Coef_Patim=0.311 #Quadrangular
                if neve=="Sim":
                    Coef_Patim_gelo=0.004*esp_gelo*10**3+0.311
                    Coef_gelo_peso=5.287*10**-5*(esp_gelo*10**3)**2+2.792*10**-3*esp_gelo*10**3
            Af=AntenasAux[4].astype(float)*10**-3*Coef_Patim
            Al=Af
            At=Af
            
            Peso=AntenasAux[5].astype(float)
            if neve=="Sim":
                Afgelo=AntenasAux[4].astype(float)*10**-3*Coef_Patim_gelo
                Algelo=Afgelo
                Atgelo=Afgelo
                mgelo=Peso+Coef_gelo_peso*rho_gelo
            
        Ant[i,1]=Peso
        Ant[i,3]=mgelo
        #Projeção das areas 
        if float(AntenasAux[0])==8:
            AreaDirection[i,:]=np.ones(len(AreaDirection[i,:]))*Af
            if neve=="Sim":
                AreaDirectiongelo[i,:]=np.ones(len(AreaDirection[i,:]))*Afgelo
        else:
            for j in range(len(Alpha_wind)):
                if dalpha[j]>=270*np.pi/180 or  dalpha[j]<=90*np.pi/180:
                    AreaDirection[i,j]=Af*(np.cos(dalpha[j]))**2+Al*(np.sin(dalpha[j]))**2
                    if neve=="Sim":
                        AreaDirectiongelo[i,j]=Afgelo*(np.cos(dalpha[j]))**2+Algelo*(np.sin(dalpha[j]))**2
                else:
                    AreaDirection[i,j]=At*(np.cos(dalpha[j]))**2+Al*(np.sin(dalpha[j]))**2
                    if neve=="Sim":
                        AreaDirectiongelo[i,j]=Atgelo*(np.cos(dalpha[j]))**2+Algelo*(np.sin(dalpha[j]))**2
        #Antout[i,0]=AntenasAux[0] Nome das antenas
    AreaAnt=np.zeros(len(Alpha_wind))
    for i in range(len(Antenas)):
        AntenasAux=np.array((Antenas[i].strip()).split(';'))
        indiceAnt=find_closest_value_position(h_vector,float(AntenasAux[2]))
        if RRU[i]==1 and Ocupation[indiceAnt]==1:
            AreaDirection[i,:]=AreaDirection[i,:]*0.5
            if neve=="Sim":
                AreaDirectiongelo[i,:]=AreaDirectiongelo[i,:]*0.5
    #Momento de area
    for i in range(len(Alpha_wind)):
        AreaAnt[i]=sum(AreaDirection[:,i]*Ant[:,0]/h_torre)
    #indice_angle=np.argmax(AreaAnt)
    #Ant[:,2]=AreaDirection[:,indice_angle]
    #if neve=="Sim":
        #Ant[:,4]=AreaDirectiongelo[:,indice_angle]
    Ant[:,5]=np.ones(len(Ant[:,5]))
    Ant=np.round(Ant,5)
    Antout = np.hstack((Antout, Ant[:,:-1]))

    return Ant, Nant, AreaAnt,AreaDirection,AreaDirectiongelo

def AntenaOrg(Antenas,Angle,Area,Org,alpha_norte):    
    Cf_frente=0.95
    Cf_lado=0.85
    Cf_tras=1.1
    Cf_RRU=1.4
    Cf_MW=1.4
    AntOrg=[]
    MOrg=[]
    alpha_norte=alpha_norte*np.pi/180
    for i in range(len(Antenas)):
        AntenasAux=np.array((Antenas[i].strip()).split(';'))
        EquipType=int(AntenasAux[0])
        ConfigType=int(AntenasAux[1])
        if EquipType!=8:
            azimute=float(AntenasAux[3])*np.pi/180
            dalpha=Angle-(azimute - alpha_norte)
            if dalpha<0:
                dalpha=dalpha+2*np.pi
            if dalpha>2*np.pi:
                dalpha=dalpha-2*np.pi

        if EquipType==1: #Antenna from a Catalog
            if ConfigType<4:
                for j in range(ConfigType):
                    Af=float(AntenasAux[8*j+10])
                    Al=float(AntenasAux[8*j+11])
                    At=float(AntenasAux[8*j+12])
                    if dalpha>=270*np.pi/180 and  dalpha<=90*np.pi/180:
                        AntOrg.append(Af*(np.cos(dalpha))**2+Al*(np.sin(dalpha))**2)
                    else:
                        AntOrg.append(At*(np.cos(dalpha))**2+Al*(np.sin(dalpha))**2)
                    MOrg.append(float(AntenasAux[8*j+13]))                                
        elif EquipType==2: #RRU from a Catalog            
            for j in range(ConfigType):
                AntOrg.append(Area[i,2]/ConfigType)
                MOrg.append(Area[i,1]/ConfigType)
        elif EquipType==3:
            AntOrg.append(Area[i,2])
            MOrg.append(Area[i,1])
        elif EquipType==4: #Other Equipment
            if ConfigType<4:
                for j in range(ConfigType):
                    L1=AntenasAux[5*j+7].astype(float)*10**-3
                    L2=AntenasAux[5*j+8].astype(float)*10**-3
                    L3=AntenasAux[5*j+9].astype(float)*10**-3
                    Af=L1*L2*Cf_frente
                    Al=L1*L3*Cf_lado
                    At=L1*L2*Cf_tras
                    if dalpha>=-90*np.pi/180 and  dalpha<=90*np.pi/180:
                        AntOrg.append(Af*(np.cos(dalpha))**2+Al*(np.sin(dalpha))**2)
                    else:
                        AntOrg.append(At*(np.cos(dalpha))**2+Al*(np.sin(dalpha))**2)   
                    MOrg.append(float(AntenasAux[5*j+10]))                 
        elif EquipType==5:
            for j in range(ConfigType):
                AntOrg.append(Area[i,2]/ConfigType)
                MOrg.append(Area[i,1]/ConfigType)
        elif EquipType==6: #Other Microwave
            AntOrg.append(Area[i,2])
            MOrg.append(Area[i,1])
        elif EquipType==7: #Working Plataform
            AntOrg.append(Area[i,2])
            MOrg.append(Area[i,1])
        elif EquipType==8:
            AntOrg.append(Area[i,2])
            MOrg.append(Area[i,1])
        elif EquipType==9:
            AntOrg.append(Area[i,2])
            MOrg.append(Area[i,1])
    OrgFinal=np.zeros(len(Org))
    Mfinal=np.zeros(len(Org))
    for i in range(len(Org)):
        auxOrg = np.array([int(x) for x in Org[i].strip().split(';')])
        for j in range(len(auxOrg)):
            OrgFinal[i]=AntOrg[auxOrg[j]-1]+OrgFinal[i]
            Mfinal[i]=MOrg[auxOrg[j]-1]+Mfinal[i]
        OrgFinal[i]=np.round(OrgFinal[i],3)
        Mfinal[i]=np.round(Mfinal[i],3)
    return AntOrg,OrgFinal,Mfinal