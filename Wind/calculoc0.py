import numpy as np
from Utilities.Utilities import find_closest_value_position

def calculoc0(Pais,zmin,h_vector,Tipo,Alt_colina,Lu,Ld,Xtopo,Ac,A500,A1000):
    if Pais in ["Portugal","Portugal-RSA"]:
        c0=calculoc0_Portugal(Pais,zmin,h_vector,Tipo,Alt_colina,Lu,Ld,Xtopo,Ac,A500,A1000)
        return c0
    if Pais in ["France"]:
        c0=calculoc0_France(Pais,zmin,h_vector,Tipo,Alt_colina,Lu,Ld,Xtopo,Ac,A500,A1000)
        return c0
    if Pais in ["Spain"]:
        c0=calculoc0_France(Pais,zmin,h_vector,Tipo,Alt_colina,Lu,Ld,Xtopo,Ac,A500,A1000)
        return c0
def calculoc0_Portugal(Pais,zmin,h_vector,Tipo,Alt_colina,Lu,Ld,Xtopo,Ac,A500,A1000):
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

def calculoc0_France(Pais,zmin,h_vector,Tipo,Alt_colina,Lu,Ld,Xtopo,Ac,A500,A1000):
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

def calculoc0_Spain(Pais,zmin,h_vector,Tipo,Alt_colina,Lu,Ld,Xtopo,Ac,A500,A1000):
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
