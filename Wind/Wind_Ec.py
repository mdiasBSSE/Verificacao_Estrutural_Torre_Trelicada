import numpy as np
from Utilities.Utilities import find_closest_value_position

def CalcKa(A_lin,As,Ab,L_equiv,L_med):
    Ka=np.ones(len(A_lin))
    iKa= (A_lin<As)& (A_lin<0.5*Ab)& (L_equiv<1.1*L_med)
    # print("Alin:", A_lin, "<As",As)
    # print((A_lin<As))
    # print("Alin:", A_lin, "<0.5Ab",0.5*Ab)
    # print((A_lin<0.5*Ab))
    # print("L_equiv:", L_equiv, "<1.1*L_med",1.1*L_med)
    # print((L_equiv<L_med))
    Ka[iKa]=0.8
    return Ka


def calculoc0(Pais,zmin,h_vector,Tipo,Alt_colina,Lu,Ld,Xtopo,Ac,A500,A1000):
    #vb,z0,zmin=Zone(Zona,Terrain,Pais)
    if Pais in ["Portugal","Portugal-RSA","United Kingdom","Germany"] and Tipo!="Colina-Metodo 1":
        phi=Alt_colina/Lu
        if phi<0.05:
            c0=np.ones(len(h_vector))
        else:
            c0=np.zeros(len(h_vector))
            for i in range(len(h_vector)):
                if phi<0.3:
                    Le=Lu
                else:
                    Le=Alt_colina/0.3
                ratio=h_vector[i]/Le
                if Xtopo/Lu<=0 and Xtopo/Lu>=-1.5 and h_vector[i]/Le>=0 and h_vector[i]/Le<=2.0:
                    A=0.1552*(h_vector[i]/Le)**4-0.8575*(h_vector[i]/Le)**3+1.8133*(h_vector[i]/Le)**2-1.9115*(h_vector[i]/Le)+1.0124
                    B=0.3542*(h_vector[i]/Le)**2-1.0577*(h_vector[i]/Le)+2.6456
                    s=A*np.exp(B*Xtopo/Lu)
                elif Xtopo/Lu<-1.5 or h_vector[i]/Le>2.0:
                    s=0
                elif Tipo=="Falésia e escarpas":
                    if ratio<0.1:
                        ratio=0.1

                    if Xtopo/Le>=0.1 and Xtopo/Le<=3.5 and ratio<=2.0:                    
                    
                        A=-1.342*(np.log10(ratio))**3-0.8222*(np.log10(ratio))**2+0.4609*np.log10(ratio)-0.0791
                        B=-1.0196*(np.log10(ratio))**3-0.8910*(np.log10(ratio))**2+0.5343*np.log10(ratio)-0.1156
                        C=0.8030*(np.log10(ratio))**3+0.4236*(np.log10(ratio))**2-0.5738*np.log10(ratio)+0.1606
                        s=A*(np.log10(Xtopo/Le))**2+B*(np.log10(Xtopo/Le))+C
                    elif Xtopo/Le<0.1 and ratio<2.0:
                        saux1=0.1552*(ratio)**4-0.8575*(ratio)**3+1.8133*(ratio)**2-1.9115*(ratio)+1.0124
                        A=-1.342*(np.log10(ratio))**3-0.8222*(np.log10(ratio))**2+0.4609*np.log10(ratio)-0.0791
                        B=-1.0196*(np.log10(ratio))**3-0.8910*(np.log10(ratio))**2+0.5343*np.log10(ratio)-0.1156
                        C=0.8030*(np.log10(ratio))**3+0.4236*(np.log10(ratio))**2-0.5738*np.log10(ratio)+0.1606
                        saux2=A*(np.log10(Xtopo/Le))**2+B*(np.log10(Xtopo/Le))+C
                        s=saux2+(Xtopo/Le-0.1)*(saux1-saux2)/(0-0.1)
                    else:
                        s=0
                    
                elif Tipo=="Colina em cadeia" or Tipo=="Colina isolada":
                    if Xtopo/Ld>=0 and Xtopo/Ld<=2.0  and ratio>=0 and ratio<2.0:
                        A=0.1552*(ratio)**4-0.8575*(ratio)**3+1.8133*(ratio)**2-1.9115*(ratio)+1.0124
                        B=-0.3056*(ratio)**2+1.0212*(ratio)-1.7637
                        
                        s=A*np.exp(B*Xtopo/Ld)
                        
                    elif  Xtopo/Ld>2.0 or ratio>2:
                        s=0
                    
                if phi<0.3:
                    c0[i]=1+2*s*phi
                else:
                    c0[i]=1+0.6*s
            
    elif Pais=="France":
        if Tipo=="Colina-Metodo 1":
            Am=(2*Ac+4*A500+4*A1000)/10
            DeltaAc=Ac-Am
            indice10=h_vector<10
            c0=1+0.004*DeltaAc*np.exp(-0.014*(h_vector-10*np.ones(len(h_vector))))
            c0[indice10]=1+0.004*DeltaAc
            indicec0=c0<1
            c0[indicec0]=1
        elif Alt_colina/Lu<0.05:
            c0=np.ones(len(h_vector))
        else:
            phi=Alt_colina/Lu
            if phi<0.25:
                L=Lu/2
            else:
                L=2*Alt_colina
            if Tipo=="Colina em cadeia":
                Smax=2.2*Alt_colina/L
                alpha=3
                kred=1.5
            elif Tipo=='Falesia e escarpas':
                Smax=1.3*Alt_colina/L
                alpha=2.5
                if Xtopo<0:
                    kred=1.5
                else:
                    kred=4
            elif Tipo=="Colina isolada":
                Smax=1.6*Alt_colina/L
                alpha=4
                kred=1.5
            h_vector=np.array(h_vector)
            c0=1+Smax*(1-np.abs(Xtopo)/(kred*L))*np.exp(-alpha*h_vector/L)
    indicec0=c0<1
    c0[indicec0]=1
    iZs=find_closest_value_position(h_vector,zmin)
    indiceZs=h_vector<zmin
    #Considera-se c0 constante abaixo do zmin
    c0min=c0[iZs]
    c0[indiceZs]=c0min
    return c0

def ZoneInverse(z0, zmin, vb, Pais):

    # --- PORTUGAL ---
    if Pais == "Portugal":
        zonas = {"A": 27, "B": 30}
        tipos = {
            "I":   (0.005, 1),
            "II":  (0.05, 3),
            "III": (0.3, 8),
            "IV":  (1, 15),
        }

        zona = next((z for z, v in zonas.items() if v == vb), None)
        tipo = next((t for t, vals in tipos.items() if vals == (z0, zmin)), None)
        return zona, tipo

    # --- PORTUGAL-RSA ---
    elif Pais == "Portugal-RSA":
        zonas = {"A": 27, "B": 30}
        tipos = {
            "I":  (0, 15),
            "II": (0, 10),
        }

        zona = next((z for z, v in zonas.items() if v == vb), None)
        tipo = next((t for t, vals in tipos.items() if vals == (z0, zmin)), None)
        return zona, tipo

    # --- FRANCE ---
    elif Pais == "France":
        zonas = {
            "1": 22, "2": 24, "3": 26, "4": 28,
            "Guadeloupe": 36, "Guyane": 17,
            "Martinique": 32, "Réunion": 34, "Mayotte": 30
        }

        tipos = {
            "0":   (0.005, 1),
            "II":  (0.05, 2),
            "IIIa": (0.2, 5),
            "IIIb": (0.5, 9),
            "IV":  (1, 15),
        }

        zona = next((z for z, v in zonas.items() if v == vb), None)
        tipo = next((t for t, vals in tipos.items() if vals == (z0, zmin)), None)
        return zona, tipo

    # --- UNITED KINGDOM ---
    elif Pais == "United Kingdom":
        # vb não é usado (é sempre 0)
        tipos = {
            "Sea-0":      (0.003, 1),
            "Country-I":  (0.01, 1),
            "Country-II": (0.05, 2),
            "Town-III":   (0.3, 5),
            "Town-IV":    (1, 10)
        }

        zona = "-"  # não existe zona no UK
        tipo = next((t for t, vals in tipos.items() if vals == (z0, zmin)), None)
        return zona, tipo

    # --- SPAIN ---
    elif Pais == "Spain":
        zonas = {"A": 26, "B": 27, "C": 29}

        tipos = {
            "I":   (0.01, 1),
            "II":  (0.05, 2),
            "III": (0.3, 5),
            "IV":  (1, 10),
            "0":   (0.003, 1),
        }

        zona = next((z for z, v in zonas.items() if v == vb), None)
        tipo = next((t for t, vals in tipos.items() if vals == (z0, zmin)), None)
        return zona, tipo

    return None, None

def Zone(Zona,Tipo,Pais):
    if Pais=="Portugal":
        if Zona=="A": #Pag155
            vb=27
        elif Zona=="B":
            vb=30
        if Tipo=="I": #Pag156
            z0=0.005
            zmin=1
        elif Tipo=="II":
            z0=0.05
            zmin=3
        elif Tipo=="III":
            z0=0.3
            zmin=8
        elif Tipo=="IV":
            z0=1
            zmin=15
    elif Pais=="Portugal-RSA":
        if Zona=="A":
            vb=27
        elif Zona=="B":
            vb=30
        if Tipo=="I":
            zmin=15
            z0=0
        elif Tipo=="II":
            zmin=10
            z0=0
    elif Pais=="France":#Pag6 Anexo
        if Zona=="1":
            vb=22
        elif Zona=="2":
            vb=24
        elif Zona=="3":  
            vb=26
        elif Zona=="4":
            vb=28
        elif Zona=="Guadeloupe":
            vb=36
        elif Zona=="Guyane":
            vb=17
        elif Zona=="Martinique":
            vb=32
        elif Zona=="Réunion":
            vb=34
        elif Zona=="Mayotte":
            vb=30
        if Tipo=="0":#Pag16 Anexo
            z0=0.005
            zmin=1
        elif Tipo=="II":
            z0=0.05
            zmin=2
        elif Tipo=="IIIa":
            z0=0.2
            zmin=5
        elif Tipo=="IIIb":
            z0=0.5
            zmin=9
        elif Tipo=="IV":
            z0=1
            zmin=15
    elif Pais=="United Kingdom":
        vb=0
        if Tipo=="Sea-0":
            z0=0.003
            zmin=1
        elif Tipo=="Country-I": #Pag22
            z0=0.01
            zmin=1
        elif Tipo=="Country-II":
            z0=0.05
            zmin=2
        elif Tipo=="Town-III":
            z0=0.3
            zmin=5
        elif Tipo=="Town-IV":
            z0=1
            zmin=10
    elif Pais=="Spain":
        if Zona=="A": #Pag155
            vb=26
        elif Zona=="B":
            vb=27
        elif Zona=="C":
            vb=29
        if Tipo=="I": #Pag156
            z0=0.01
            zmin=1
        elif Tipo=="II":
            z0=0.05
            zmin=2
        elif Tipo=="III":
            z0=0.3
            zmin=5
        elif Tipo=="IV":
            z0=1
            zmin=10
        elif Tipo=="0":
            z0=0.003
            zmin=1
    if Zona=="-":
        vb=0
    return vb,z0,zmin
def Iv_calc(z,z0,zmin,Alt,Pais,c0):
    #Confirmar valores para portugal
    #Props do terreno
    #vb,z0,zmin=Zone(zona,tipo,Pais)
    z02=0.05
    
    Iv=np.zeros(len(z))
    ki=1
    if Pais!="Portugal-RSA":
        Iv=np.zeros(len(z))
        mask=z>0
        Iv[mask]=ki/(c0[mask]*np.log(z[mask]/z0))
        Iv[np.where(z<zmin)]= ki/(c0[np.where(z<zmin)]*np.log(zmin/z0))
    return Iv   
def Velocidade_media(z,vb,z0,zmin,Pais,c0):
    #Props do terreno
    #Confirmar valores Portugal
    #vb,z0,zmin=Zone(zona,tipo,Pais)
    if Pais!="Portugal-RSA":
        z02=0.05
        Cr=np.zeros(len(z))
        kr=0.19*(z0/z02)**0.07
        Cr=np.zeros(len(z))
        mask=z>0
        Cr[mask]=kr*np.log(z[mask]/z0)
        Cr[np.where(z<zmin)]=kr*np.log(zmin/z0)
        vm=Cr*c0*vb
    else:
        vm=np.zeros(len(z))
        #if z0>=0.5:
        if zmin==15:
            m=18
            n=0.28
        #elif z0<0.5:
        elif zmin==10:
            m=25
            n=0.2
        vmmin=m*(zmin/10)**n
        for i in range(len(z)):
            if z[i]<zmin:
                vm[i]=vmmin
            else:
                vm[i]=m*(z[i]/10)**n
        if vb==30:
            vm=1.1*vm
    return vm
def Pressure_wind(h_vector,vb,z0,zmin,Alt,Pais,rho_ar,c0):
    #Portugal_Eurocode
    Wind_Pressure=np.zeros(len(h_vector))
    vm=Velocidade_media(h_vector,vb,z0,zmin,Pais,c0)
    if Pais!="Portugal-RSA":
        IV=Iv_calc(h_vector, z0,zmin,Alt,Pais,c0)
        Wind_Pressure=(1+7*IV)*0.5*rho_ar*vm**2
    else:
        if vb==30:
            vraj=(vm/1.1+14)*1.1
        else:
            vraj=vm+14
        Wind_Pressure=0.613*vraj**2
    return Wind_Pressure


def K2_func(phi):
    phi = np.asarray(phi)  # garante que funcione para listas, arrays ou valores únicos
    
    condicoes = [
        (phi >= 0) & (phi <= 0.2),
        (phi >= 0.8) & (phi <= 1.0),
        (phi > 0.2) & (phi <= 0.5),
        (phi > 0.5) & (phi < 0.8)
    ]
    
    valores = [
        0.2,
        0.2,
        phi,
        1 - phi
    ]
    
    return np.select(condicoes, valores, default=np.nan)  
def Coef_forca(Af,Ac,Acsup,As,theta,Indice_Cheios,torre):
    #BS 1993-3-1-2006
    if torre=="Triangular":
        C1=1.9
        C2=1.4
    elif torre=="Quadrangular":
        C1=2.25
        C2=1.5
    cf0f=1.76*C1*(1-C2*Indice_Cheios+Indice_Cheios**2)
    cf0c=C1*(1-C2*Indice_Cheios)+(C1+0.875)*Indice_Cheios**2
    cf0csup=1.9-np.sqrt((1-Indice_Cheios)*(2.8-1.14*C1+Indice_Cheios))
    cfs0=np.zeros(len(cf0csup))
    K1=np.zeros(len(cf0csup))
    K2=np.zeros(len(cf0csup))
    Ktheta=np.zeros(len(cf0csup))
    cfs0[1:]=cf0f[1:]*Af[1:]/As[1:]+cf0c[1:]*Ac[1:]/As[1:]+cf0csup[1:]*Acsup[1:]/As[1:]
    K1[1:]=0.55*Af[1:]/As[1:]+0.8*(Ac[1:]+Acsup[1:])/As[1:]
    K2[1:]=K2_func(Indice_Cheios[1:])
    if torre=="Triangular":
        Ktheta[1:]=(Ac[1:]+Acsup[1:])/As[1:]+Af[1:]/As[1:]*(1-0.1*(np.sin(np.radians(1.5*theta)))**2)
    elif torre=="Quadrangular":
        Ktheta[1:]=1.0+K1[1:]*K2[1:]*(np.sin(np.radians(2*theta)))**2
    cfs=Ktheta*cfs0
    return cf0f,cf0c,cf0csup,cfs0,cfs

def CsCd_calc(Zs,h_torre,IV,Tube,z0,zmin,n1,me,cf,rho_ar,vm,Alt,Pais):
    ###Tube é um vetor em que 0 é a largura da base, 1 largura a 0.6h e 2 (-1) a largura no topo
    if Pais!="Portugal-RSA":
        T=600 #Pag 109
        # vb,z0,zmin=Zone(zone,Terrain,Pais)
        Lt=300
        Zt=200
        alpha=0.67+0.05*np.log(z0)
        L=Lt*((Zs)/Zt)**alpha
        Lmin=Lt*(zmin/Zt)**alpha
        if L<Lmin:
            L=Lmin
        fl=n1*L/vm #Pag 107
        me=me
    
        Sl=6.8*fl/((1+10.2*fl)**(5/3))
        B2=(1+0.9*((Tube[1]+h_torre)/L)**0.63)**-1 #Pag 108
        niuh=4.6*h_torre*fl/L
        niub=4.6*Tube[1]*fl/L
        Rh=1/niuh-1/(2*niuh**2)*(1-np.exp(-2*niuh))
        Rb=1/niub-1/(2*niub**2)*(1-np.exp(-2*niub))
        delta=LogDecaiment(cf,rho_ar,vm,n1,me,Tube[-1])
        R2=np.pi**2/(2*delta)*Sl*Rh*Rb #Pag 109
        v=n1*np.sqrt(R2/(B2+R2))
    
        Kp=np.sqrt(2*np.log(v*T))+0.6/(np.sqrt(2*np.log(v*T))) #Pag 109
        if Kp<3:
            Kp=3
    
        
        Cs=(1+7*IV*np.sqrt(B2))/(1+7*IV)
        Cd=(1+2*Kp*IV*np.sqrt(B2+R2))/(1+7*IV*np.sqrt(B2))
        CsCd=Cs*Cd

        n1=n1
        CsCd_Out=[["n1",np.round(n1,3)],["Zs",np.round(Zs,3)],["vm(Zs)",np.round(vm,3)],["Iv(Zs)",np.round(IV,3)],["L(Zs)",np.round(L,3)],["Sl(Zs,n1)",np.round(Sl,3)],["fl(Zs,n1,x)",np.round(fl,3)],
              ["B^2",np.round(B2,3)],["R^2",np.round(R2,3)],["niuh",np.round(niuh,3)],["niub",np.round(niub,3)],["Rh",np.round(Rh,3)],["Rb",np.round(Rb,3)],["kp",np.round(Kp,3)],["v",np.round(v,3)],["Delta",np.round(delta,3)],["me",np.round(me,3)],
              ["cf",np.round(np.nanmin(cf),3)],["Cs",np.round(Cs,3)],["Cd",np.round(Cd,3)],["CsCd",np.round(CsCd,3)]]
    else:
        CsCd=1  
        CsCd_Out=[[]]  

    return CsCd,CsCd_Out
def LogDecaiment(cf,rho_ar,vm,n1,me,b):
    deltas=0.05 #DecLog_estrutural torre aço
    #deltas=0.03 #DecLog_estrutural torre betão
    deltad=0 #DecLog_devido_a_dispositivos
    cfmin=np.nanmin(cf[1:])
    deltaa=cfmin*rho_ar*b*vm/(2*n1*me) #DecLog_aerodinamico Pag.151
    delta=deltas+deltaa+deltad
    return delta
