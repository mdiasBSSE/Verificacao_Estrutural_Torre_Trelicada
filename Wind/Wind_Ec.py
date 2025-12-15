import numpy as np


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


def LogDecaiment(cf,rho_ar,vm,n1,me,b):
    deltas=0.05 #DecLog_estrutural torre aço
    #deltas=0.03 #DecLog_estrutural torre betão
    deltad=0 #DecLog_devido_a_dispositivos
    cfmin=np.nanmin(cf[1:])
    deltaa=cfmin*rho_ar*b*vm/(2*n1*me) #DecLog_aerodinamico Pag.151
    delta=deltas+deltaa+deltad
    return delta
