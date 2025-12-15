import numpy as np
def Iv_calc(z,z0,zmin,Alt,Pais,c0):
    if Pais in ["Portugal","Portugal-RSA"]:
        Iv=Iv_calc_Portugal(z,z0,zmin,Alt,Pais,c0)
        return Iv   
    if Pais in ["France"]:
        Iv=Iv_calc_France(z,z0,zmin,Alt,Pais,c0)
        return Iv
    if Pais in ["Spain"]:
        Iv=Iv_calc_Spain(z,z0,zmin,Alt,Pais,c0)
        return Iv 
def Iv_calc_Portugal(z,z0,zmin,Alt,Pais,c0):
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
def Iv_calc_France(z,z0,zmin,Alt,Pais,c0):
    #Confirmar valores para portugal
    #Props do terreno
    #vb,z0,zmin=Zone(zona,tipo,Pais)
    z02=0.05
    
    Iv=np.zeros(len(z))
    ki=1
    Iv=np.zeros(len(z))
    mask=z>0
    Iv[mask]=ki/(c0[mask]*np.log(z[mask]/z0))
    Iv[np.where(z<zmin)]= ki/(c0[np.where(z<zmin)]*np.log(zmin/z0))
    return Iv  
def Iv_calc_Spain(z,z0,zmin,Alt,Pais,c0):
    #Confirmar valores para portugal
    #Props do terreno
    #vb,z0,zmin=Zone(zona,tipo,Pais)
    z02=0.05
    Iv=np.zeros(len(z))
    ki=1
    Iv=np.zeros(len(z))
    mask=z>0
    Iv[mask]=ki/(c0[mask]*np.log(z[mask]/z0))
    Iv[np.where(z<zmin)]= ki/(c0[np.where(z<zmin)]*np.log(zmin/z0))
    return Iv  