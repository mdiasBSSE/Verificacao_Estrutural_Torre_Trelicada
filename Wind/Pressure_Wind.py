import numpy as np
from Wind.Velocidade_Media import Velocidade_media
from Wind.Iv_Calc import Iv_calc
def Pressure_wind(h_vector,vb,z0,zmin,Alt,Pais,rho_ar,c0):
    ###Código para o calculo da pressão do vento
    if Pais in ["Portugal","Portugal-RSA"]:
        Wind_Pressure=Pressure_wind_Portugal(h_vector,vb,z0,zmin,Alt,Pais,rho_ar,c0)
        return Wind_Pressure
    elif Pais in ["France"]:
        Wind_Pressure=Pressure_wind_France(h_vector,vb,z0,zmin,Alt,Pais,rho_ar,c0)
        return Wind_Pressure
    elif Pais in ["Spain"]:
        Wind_Pressure=Pressure_wind_France(h_vector,vb,z0,zmin,Alt,Pais,rho_ar,c0)
        return Wind_Pressure

def Pressure_wind_France(h_vector,vb,z0,zmin,Alt,Pais,rho_ar,c0):
    #Portugal_Eurocode
    Wind_Pressure=np.zeros(len(h_vector))
    vm=Velocidade_media(h_vector,vb,z0,zmin,Pais,c0)
    IV=Iv_calc(h_vector, z0,zmin,Alt,Pais,c0)
    Wind_Pressure=(1+7*IV)*0.5*rho_ar*vm**2
    return Wind_Pressure
def Pressure_wind_Portugal(h_vector,vb,z0,zmin,Alt,Pais,rho_ar,c0):
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
def Pressure_wind_Spain(h_vector,vb,z0,zmin,Alt,Pais,rho_ar,c0):
    #Portugal_Eurocode
    Wind_Pressure=np.zeros(len(h_vector))
    vm=Velocidade_media(h_vector,vb,z0,zmin,Pais,c0)
    IV=Iv_calc(h_vector, z0,zmin,Alt,Pais,c0)
    Wind_Pressure=(1+7*IV)*0.5*rho_ar*vm**2
    return Wind_Pressure