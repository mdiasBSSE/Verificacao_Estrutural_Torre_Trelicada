import numpy as np
def Velocidade_media(z,vb,z0,zmin,Pais,c0):
    #Props do terreno
    #Confirmar valores Portugal
    if Pais in ["Portugal","Portugal-RSA"]:
        vm=Velocidade_media_Portugal(z,vb,z0,zmin,Pais,c0)
        return vm
    elif Pais in ["France"]:
        vm=Velocidade_media_France(z,vb,z0,zmin,Pais,c0)
        return vm
    elif Pais in ["Spain"]:
        vm=Velocidade_media_Spain(z,vb,z0,zmin,Pais,c0)
        return vm

def Velocidade_media_Portugal(z,vb,z0,zmin,Pais,c0):
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
def Velocidade_media_France(z,vb,z0,zmin,Pais,c0):
    #Props do terreno
    #vb,z0,zmin=Zone(zona,tipo,Pais)
    z02=0.05
    Cr=np.zeros(len(z))
    kr=0.19*(z0/z02)**0.07
    Cr=np.zeros(len(z))
    mask=z>0
    Cr[mask]=kr*np.log(z[mask]/z0)
    Cr[np.where(z<zmin)]=kr*np.log(zmin/z0)
    vm=Cr*c0*vb
    return vm
def Velocidade_media_Spain(z,vb,z0,zmin,Pais,c0):
    #Props do terreno
    #vb,z0,zmin=Zone(zona,tipo,Pais)
    z02=0.05
    Cr=np.zeros(len(z))
    kr=0.19*(z0/z02)**0.07
    Cr=np.zeros(len(z))
    mask=z>0
    Cr[mask]=kr*np.log(z[mask]/z0)
    Cr[np.where(z<zmin)]=kr*np.log(zmin/z0)
    vm=Cr*c0*vb
    return vm