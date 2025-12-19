import numpy as np 
from Wind.Wind_Ec import LogDecaiment

def CsCd_calc(Zs,h_torre,IV,Tube,z0,zmin,n1,me,cf,rho_ar,vm,Alt,Pais):
    ###Codigo para calcular o efeito dinamico na estrutura (CsCd no Eurocodigo)
    if Pais in ["Portugal"]:
        CsCd,CsCd_Out=CsCd_calc_Portugal(Zs,h_torre,IV,Tube,z0,zmin,n1,me,cf,rho_ar,vm,Alt,Pais)
        return CsCd,CsCd_Out
    elif Pais in ["France"]:
        CsCd,CsCd_Out=CsCd_calc_France(Zs,h_torre,IV,Tube,z0,zmin,n1,me,cf,rho_ar,vm,Alt,Pais)
        return CsCd,CsCd_Out
    elif Pais in ["Spain"]:
        CsCd,CsCd_Out=CsCd_calc_Spain(Zs,h_torre,IV,Tube,z0,zmin,n1,me,cf,rho_ar,vm,Alt,Pais)
        return CsCd,CsCd_Out
    elif Pais in ["Italy"]:
        CsCd,CsCd_Out=CsCd_calc_Spain(Zs,h_torre,IV,Tube,z0,zmin,n1,me,cf,rho_ar,vm,Alt,Pais)
        return CsCd,CsCd_Out
    elif Pais in ["Portugal-RSA"]:
        CsCd=1
        CsCd_Out=[[]]
        return CsCd,CsCd_Out

def CsCd_calc_Portugal(Zs,h_torre,IV,Tube,z0,zmin,n1,me,cf,rho_ar,vm,Alt,Pais):
    ###Tube é um vetor em que 0 é a largura da base, 1 largura a 0.6h e 2 (-1) a largura no topo
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
    return CsCd,CsCd_Out
def CsCd_calc_Italy(Zs,h_torre,IV,Tube,z0,zmin,n1,me,cf,rho_ar,vm,Alt,Pais):
    ###Tube é um vetor em que 0 é a largura da base, 1 largura a 0.6h e 2 (-1) a largura no topo
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
    return CsCd,CsCd_Out
def CsCd_calc_France(Zs,h_torre,IV,Tube,z0,zmin,n1,me,cf,rho_ar,vm,Alt,Pais):
    ###Tube é um vetor em que 0 é a largura da base, 1 largura a 0.6h e 2 (-1) a largura no topo
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
    return CsCd,CsCd_Out
def CsCd_calc_Spain(Zs,h_torre,IV,Tube,z0,zmin,n1,me,cf,rho_ar,vm,Alt,Pais):
    ###Tube é um vetor em que 0 é a largura da base, 1 largura a 0.6h e 2 (-1) a largura no topo
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
    return CsCd,CsCd_Out