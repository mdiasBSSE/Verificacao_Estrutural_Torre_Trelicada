import numpy as np
def Coef_forca(Af,Ac,Acsup,As,theta,Indice_Cheios,torre,Pais):
    ##### Caclculo Cf da estrutura
    if Pais in ["Portugal","Portugal-RSA"]:
        cf0f,cf0c,cf0csup,cfs0,cfs=Coef_forca_Portugal(Af,Ac,Acsup,As,theta,Indice_Cheios,torre)
        return cf0f,cf0c,cf0csup,cfs0,cfs
    if Pais in ["France"]:
        cf0f,cf0c,cf0csup,cfs0,cfs=Coef_forca_France(Af,Ac,Acsup,As,theta,Indice_Cheios,torre)
        return cf0f,cf0c,cf0csup,cfs0,cfs
    if Pais in ["Spain"]:
        cf0f,cf0c,cf0csup,cfs0,cfs=Coef_forca_Spain(Af,Ac,Acsup,As,theta,Indice_Cheios,torre)
        return cf0f,cf0c,cf0csup,cfs0,cfs
def Coef_forca_Portugal(Af,Ac,Acsup,As,theta,Indice_Cheios,torre):
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
def Coef_forca_France(Af,Ac,Acsup,As,theta,Indice_Cheios,torre):
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
def Coef_forca_Spain(Af,Ac,Acsup,As,theta,Indice_Cheios,torre):
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