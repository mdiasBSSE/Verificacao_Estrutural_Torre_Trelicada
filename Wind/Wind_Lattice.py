import numpy as np 
from Utilities.Utilities_ops import no_mais_proximo, distancia_entre_nos,elementos_do_no

def Elements_Areas(Pais,nos_por_classe,nos_vento,Troco,TrussType,Dim,z_vento,ElementosID,ExpVento,Comprimento_Barra,divisor,Reynolds,esp_gelo):
    ####Calcula Area cantoneiras/ Tubos, etc.
    if Pais in ["Portugal","Portugal-RSA"]:
        nos_vento_troco,Largura,Largura_total,Largura_media,At,Af,Ac,Acsup,As,Indice_Cheios=Elements_Areas_Portugal(nos_por_classe,nos_vento,Troco,TrussType,Dim,z_vento,ElementosID,ExpVento,Comprimento_Barra,divisor,Reynolds,esp_gelo)
        return nos_vento_troco,Largura,Largura_total,Largura_media,At,Af,Ac,Acsup,As,Indice_Cheios
    elif Pais in ["France"]:
        nos_vento_troco,Largura,Largura_total,Largura_media,At,Af,Ac,Acsup,As,Indice_Cheios=Elements_Areas_France(nos_por_classe,nos_vento,Troco,TrussType,Dim,z_vento,ElementosID,ExpVento,Comprimento_Barra,divisor,Reynolds,esp_gelo)
        return nos_vento_troco,Largura,Largura_total,Largura_media,At,Af,Ac,Acsup,As,Indice_Cheios
    elif Pais in ["Spain"]:
        nos_vento_troco,Largura,Largura_total,Largura_media,At,Af,Ac,Acsup,As,Indice_Cheios=Elements_Areas_Spain(nos_por_classe,nos_vento,Troco,TrussType,Dim,z_vento,ElementosID,ExpVento,Comprimento_Barra,divisor,Reynolds,esp_gelo)
        return nos_vento_troco,Largura,Largura_total,Largura_media,At,Af,Ac,Acsup,As,Indice_Cheios


def Elements_Areas_France(nos_por_classe,nos_vento,Troco,TrussType,Dim,z_vento,ElementosID,ExpVento,Comprimento_Barra,divisor,Reynolds,esp_gelo):
    nos_vento_troco=np.zeros((2,2,len(nos_por_classe)+1))
    Barras_vento_troco=np.zeros((2,2,len(nos_por_classe)+1))
    Largura=np.zeros((2,len(nos_por_classe)+1))
    At=np.zeros(len(nos_por_classe)+1)
    Af=np.zeros(len(nos_por_classe)+1)
    Ac=np.zeros(len(nos_por_classe)+1)
    Acsup=np.zeros(len(nos_por_classe)+1)
    Largura_media=np.zeros(len(nos_por_classe)+1)
    Largura_total=np.zeros((2,len(nos_por_classe)+1))

    for i in range(nos_vento.shape[0]-1):
        nos_vento_troco[0,0,i+1]=nos_vento[i,0] ###No topo 0
        nos_vento_troco[0,1,i+1]=nos_vento[i,1] ###No topo 1
        nos_vento_troco[1,0,i+1],dist0=no_mais_proximo(int(nos_vento_troco[0,0,i+1]),nos_vento[i+1,:]) ###No base 0
        nos_vento_troco[1,1,i+1],dist1=no_mais_proximo(int(nos_vento_troco[0,1,i+1]),nos_vento[i+1,:]) ###No base 1

        Largura[0,i+1]=distancia_entre_nos(int(nos_vento_troco[0,0,i+1]),int(nos_vento_troco[0,1,i+1])) ### Largura entre eixos topo 
        Largura[1,i+1]=distancia_entre_nos(int(nos_vento_troco[1,0,i+1]),int(nos_vento_troco[1,1,i+1])) ### Largura entre eixos base
        for j in range(2):
            for k in range(2):
                Lista_barra=np.array(elementos_do_no(nos_vento_troco[k,j,i+1]))
                Mascara=(Troco[Lista_barra] == i + 1) & (TrussType[Lista_barra] == "Leg") ###Mascara para filtrar montante inicial e final do troço
                Barras_vento_troco[k,j,i+1]=Lista_barra[Mascara]
        Largura_total[0,i+1]=Largura[0,i+1]+ Dim[int(Barras_vento_troco[0,0,i+1])]/2+ Dim[int(Barras_vento_troco[0,1,i+1])]/2+2*esp_gelo
        Largura_total[1,i+1]=Largura[1,i+1]+ Dim[int(Barras_vento_troco[1,0,i+1])]/2+ Dim[int(Barras_vento_troco[1,1,i+1])]/2+2*esp_gelo
        Largura_media[i+1]=np.mean(Largura_total[:,i+1])
        At[i+1]=Largura_media[i+1]*(z_vento[i]-z_vento[i+1])
        #### Asssumindo que os 3/4 lados da  estrutura é igual####

        Mascara = ((Troco == i + 1) &np.isin(TrussType, ["Horizontal Bar", "Diagonal Bar"]) &(ExpVento == "Flat"))###Mascara para filtrar barras flat da seccção
        idx = ElementosID[Mascara[1:]].astype(int)
        Af[i+1] =Af[i+1]+ np.sum((Dim[idx]+2*esp_gelo) * Comprimento_Barra[idx])/divisor
        Mascara = ((Troco == i + 1) &np.isin(TrussType, ["Leg"]) &(ExpVento == "Flat"))###Mascara para filtrar barras flat da seccção
        idx = ElementosID[Mascara[1:]].astype(int)
        Af[i+1] =Af[i+1]+ np.sum((Dim[idx]+2*esp_gelo)* Comprimento_Barra[idx])*2/divisor
        Mascara= ((Troco == i + 1) &np.isin(TrussType, ["Horizontal Bar", "Diagonal Bar"]) &(ExpVento == "Circ") & (Reynolds<=4*10**5)) ###Mascara para filtrar barras horizotais e diagonais circulares em regime subcritico da seccção
        idx = ElementosID[Mascara[1:]].astype(int)
        Ac[i+1]=Ac[i+1] + np.sum((Dim[idx]+2*esp_gelo) * Comprimento_Barra[idx])/divisor
        Mascara= ((Troco == i + 1) &np.isin(TrussType, ["Leg"]) &(ExpVento == "Circ") & (Reynolds<=4*10**5)) ###Mascara para filtrar montantes circulares em regime subcritico da seccção
        idx = ElementosID[Mascara[1:]].astype(int)
        Ac[i+1]=Ac[i+1] + np.sum((Dim[idx]+2*esp_gelo) * Comprimento_Barra[idx])*2/divisor
        Mascara= ((Troco == i + 1) &np.isin(TrussType, ["Horizontal Bar", "Diagonal Bar"]) &(ExpVento == "Circ") & (Reynolds>4*10**5)) ###Mascara para filtrar barras horizotais e diagonais circulares em regime subcritico da seccção
        idx = ElementosID[Mascara[1:]].astype(int)
        Acsup[i+1]=np.sum((Dim[idx]+2*esp_gelo) * Comprimento_Barra[idx])/divisor
        Mascara= ((Troco == i + 1) &np.isin(TrussType, ["Leg"]) &(ExpVento == "Circ") & (Reynolds>4*10**5)) ###Mascara para filtrar montantes circulares em regime subcritico da seccção
        idx = ElementosID[Mascara[1:]].astype(int)
        Acsup[i+1]=Acsup[i+1] + np.sum((Dim[idx]+2*esp_gelo) * Comprimento_Barra[idx])*2/divisor

    As=Af+Ac+Acsup
    Indice_Cheios=np.zeros(len(As))
    Indice_Cheios[1:]=As[1:]/At[1:]
    return nos_vento_troco,Largura,Largura_total,Largura_media,At,Af,Ac,Acsup,As,Indice_Cheios

def Elements_Areas_Portugal(nos_por_classe,nos_vento,Troco,TrussType,Dim,z_vento,ElementosID,ExpVento,Comprimento_Barra,divisor,Reynolds,esp_gelo):
    nos_vento_troco=np.zeros((2,2,len(nos_por_classe)+1))
    Barras_vento_troco=np.zeros((2,2,len(nos_por_classe)+1))
    Largura=np.zeros((2,len(nos_por_classe)+1))
    At=np.zeros(len(nos_por_classe)+1)
    Af=np.zeros(len(nos_por_classe)+1)
    Ac=np.zeros(len(nos_por_classe)+1)
    Acsup=np.zeros(len(nos_por_classe)+1)
    Largura_media=np.zeros(len(nos_por_classe)+1)
    Largura_total=np.zeros((2,len(nos_por_classe)+1))

    for i in range(nos_vento.shape[0]-1):
        nos_vento_troco[0,0,i+1]=nos_vento[i,0] ###No topo 0
        nos_vento_troco[0,1,i+1]=nos_vento[i,1] ###No topo 1
        nos_vento_troco[1,0,i+1],dist0=no_mais_proximo(int(nos_vento_troco[0,0,i+1]),nos_vento[i+1,:]) ###No base 0
        nos_vento_troco[1,1,i+1],dist1=no_mais_proximo(int(nos_vento_troco[0,1,i+1]),nos_vento[i+1,:]) ###No base 1

        Largura[0,i+1]=distancia_entre_nos(int(nos_vento_troco[0,0,i+1]),int(nos_vento_troco[0,1,i+1])) ### Largura entre eixos topo 
        Largura[1,i+1]=distancia_entre_nos(int(nos_vento_troco[1,0,i+1]),int(nos_vento_troco[1,1,i+1])) ### Largura entre eixos base
        for j in range(2):
            for k in range(2):
                Lista_barra=np.array(elementos_do_no(nos_vento_troco[k,j,i+1]))
                Mascara=(Troco[Lista_barra] == i + 1) & (TrussType[Lista_barra] == "Leg") ###Mascara para filtrar montante inicial e final do troço
                Barras_vento_troco[k,j,i+1]=Lista_barra[Mascara]
        Largura_total[0,i+1]=Largura[0,i+1]+ Dim[int(Barras_vento_troco[0,0,i+1])]/2+ Dim[int(Barras_vento_troco[0,1,i+1])]/2+2*esp_gelo
        Largura_total[1,i+1]=Largura[1,i+1]+ Dim[int(Barras_vento_troco[1,0,i+1])]/2+ Dim[int(Barras_vento_troco[1,1,i+1])]/2+2*esp_gelo
        Largura_media[i+1]=np.mean(Largura_total[:,i+1])
        At[i+1]=Largura_media[i+1]*(z_vento[i]-z_vento[i+1])
        #### Asssumindo que os 3/4 lados da  estrutura é igual####

        Mascara = ((Troco == i + 1) &np.isin(TrussType, ["Horizontal Bar", "Diagonal Bar"]) &(ExpVento == "Flat"))###Mascara para filtrar barras flat da seccção
        idx = ElementosID[Mascara[1:]].astype(int)
        Af[i+1] =Af[i+1]+ np.sum((Dim[idx]+2*esp_gelo) * Comprimento_Barra[idx])/divisor
        Mascara = ((Troco == i + 1) &np.isin(TrussType, ["Leg"]) &(ExpVento == "Flat"))###Mascara para filtrar barras flat da seccção
        idx = ElementosID[Mascara[1:]].astype(int)
        Af[i+1] =Af[i+1]+ np.sum((Dim[idx]+2*esp_gelo)* Comprimento_Barra[idx])*2/divisor
        Mascara= ((Troco == i + 1) &np.isin(TrussType, ["Horizontal Bar", "Diagonal Bar"]) &(ExpVento == "Circ") & (Reynolds<=4*10**5)) ###Mascara para filtrar barras horizotais e diagonais circulares em regime subcritico da seccção
        idx = ElementosID[Mascara[1:]].astype(int)
        Ac[i+1]=Ac[i+1] + np.sum((Dim[idx]+2*esp_gelo) * Comprimento_Barra[idx])/divisor
        Mascara= ((Troco == i + 1) &np.isin(TrussType, ["Leg"]) &(ExpVento == "Circ") & (Reynolds<=4*10**5)) ###Mascara para filtrar montantes circulares em regime subcritico da seccção
        idx = ElementosID[Mascara[1:]].astype(int)
        Ac[i+1]=Ac[i+1] + np.sum((Dim[idx]+2*esp_gelo) * Comprimento_Barra[idx])*2/divisor
        Mascara= ((Troco == i + 1) &np.isin(TrussType, ["Horizontal Bar", "Diagonal Bar"]) &(ExpVento == "Circ") & (Reynolds>4*10**5)) ###Mascara para filtrar barras horizotais e diagonais circulares em regime subcritico da seccção
        idx = ElementosID[Mascara[1:]].astype(int)
        Acsup[i+1]=np.sum((Dim[idx]+2*esp_gelo) * Comprimento_Barra[idx])/divisor
        Mascara= ((Troco == i + 1) &np.isin(TrussType, ["Leg"]) &(ExpVento == "Circ") & (Reynolds>4*10**5)) ###Mascara para filtrar montantes circulares em regime subcritico da seccção
        idx = ElementosID[Mascara[1:]].astype(int)
        Acsup[i+1]=Acsup[i+1] + np.sum((Dim[idx]+2*esp_gelo) * Comprimento_Barra[idx])*2/divisor

    As=Af+Ac+Acsup
    Indice_Cheios=np.zeros(len(As))
    Indice_Cheios[1:]=As[1:]/At[1:]
    return nos_vento_troco,Largura,Largura_total,Largura_media,At,Af,Ac,Acsup,As,Indice_Cheios

def Elements_Areas_Spain(nos_por_classe,nos_vento,Troco,TrussType,Dim,z_vento,ElementosID,ExpVento,Comprimento_Barra,divisor,Reynolds,esp_gelo):
    nos_vento_troco=np.zeros((2,2,len(nos_por_classe)+1))
    Barras_vento_troco=np.zeros((2,2,len(nos_por_classe)+1))
    Largura=np.zeros((2,len(nos_por_classe)+1))
    At=np.zeros(len(nos_por_classe)+1)
    Af=np.zeros(len(nos_por_classe)+1)
    Ac=np.zeros(len(nos_por_classe)+1)
    Acsup=np.zeros(len(nos_por_classe)+1)
    Largura_media=np.zeros(len(nos_por_classe)+1)
    Largura_total=np.zeros((2,len(nos_por_classe)+1))

    for i in range(nos_vento.shape[0]-1):
        nos_vento_troco[0,0,i+1]=nos_vento[i,0] ###No topo 0
        nos_vento_troco[0,1,i+1]=nos_vento[i,1] ###No topo 1
        nos_vento_troco[1,0,i+1],dist0=no_mais_proximo(int(nos_vento_troco[0,0,i+1]),nos_vento[i+1,:]) ###No base 0
        nos_vento_troco[1,1,i+1],dist1=no_mais_proximo(int(nos_vento_troco[0,1,i+1]),nos_vento[i+1,:]) ###No base 1

        Largura[0,i+1]=distancia_entre_nos(int(nos_vento_troco[0,0,i+1]),int(nos_vento_troco[0,1,i+1])) ### Largura entre eixos topo 
        Largura[1,i+1]=distancia_entre_nos(int(nos_vento_troco[1,0,i+1]),int(nos_vento_troco[1,1,i+1])) ### Largura entre eixos base
        for j in range(2):
            for k in range(2):
                Lista_barra=np.array(elementos_do_no(nos_vento_troco[k,j,i+1]))
                Mascara=(Troco[Lista_barra] == i + 1) & (TrussType[Lista_barra] == "Leg") ###Mascara para filtrar montante inicial e final do troço
                Barras_vento_troco[k,j,i+1]=Lista_barra[Mascara]
        Largura_total[0,i+1]=Largura[0,i+1]+ Dim[int(Barras_vento_troco[0,0,i+1])]/2+ Dim[int(Barras_vento_troco[0,1,i+1])]/2+2*esp_gelo
        Largura_total[1,i+1]=Largura[1,i+1]+ Dim[int(Barras_vento_troco[1,0,i+1])]/2+ Dim[int(Barras_vento_troco[1,1,i+1])]/2+2*esp_gelo
        Largura_media[i+1]=np.mean(Largura_total[:,i+1])
        At[i+1]=Largura_media[i+1]*(z_vento[i]-z_vento[i+1])
        #### Asssumindo que os 3/4 lados da  estrutura é igual####

        Mascara = ((Troco == i + 1) &np.isin(TrussType, ["Horizontal Bar", "Diagonal Bar"]) &(ExpVento == "Flat"))###Mascara para filtrar barras flat da seccção
        idx = ElementosID[Mascara[1:]].astype(int)
        Af[i+1] =Af[i+1]+ np.sum((Dim[idx]+2*esp_gelo) * Comprimento_Barra[idx])/divisor
        Mascara = ((Troco == i + 1) &np.isin(TrussType, ["Leg"]) &(ExpVento == "Flat"))###Mascara para filtrar barras flat da seccção
        idx = ElementosID[Mascara[1:]].astype(int)
        Af[i+1] =Af[i+1]+ np.sum((Dim[idx]+2*esp_gelo)* Comprimento_Barra[idx])*2/divisor
        Mascara= ((Troco == i + 1) &np.isin(TrussType, ["Horizontal Bar", "Diagonal Bar"]) &(ExpVento == "Circ") & (Reynolds<=4*10**5)) ###Mascara para filtrar barras horizotais e diagonais circulares em regime subcritico da seccção
        idx = ElementosID[Mascara[1:]].astype(int)
        Ac[i+1]=Ac[i+1] + np.sum((Dim[idx]+2*esp_gelo) * Comprimento_Barra[idx])/divisor
        Mascara= ((Troco == i + 1) &np.isin(TrussType, ["Leg"]) &(ExpVento == "Circ") & (Reynolds<=4*10**5)) ###Mascara para filtrar montantes circulares em regime subcritico da seccção
        idx = ElementosID[Mascara[1:]].astype(int)
        Ac[i+1]=Ac[i+1] + np.sum((Dim[idx]+2*esp_gelo) * Comprimento_Barra[idx])*2/divisor
        Mascara= ((Troco == i + 1) &np.isin(TrussType, ["Horizontal Bar", "Diagonal Bar"]) &(ExpVento == "Circ") & (Reynolds>4*10**5)) ###Mascara para filtrar barras horizotais e diagonais circulares em regime subcritico da seccção
        idx = ElementosID[Mascara[1:]].astype(int)
        Acsup[i+1]=np.sum((Dim[idx]+2*esp_gelo) * Comprimento_Barra[idx])/divisor
        Mascara= ((Troco == i + 1) &np.isin(TrussType, ["Leg"]) &(ExpVento == "Circ") & (Reynolds>4*10**5)) ###Mascara para filtrar montantes circulares em regime subcritico da seccção
        idx = ElementosID[Mascara[1:]].astype(int)
        Acsup[i+1]=Acsup[i+1] + np.sum((Dim[idx]+2*esp_gelo) * Comprimento_Barra[idx])*2/divisor

    As=Af+Ac+Acsup
    Indice_Cheios=np.zeros(len(As))
    Indice_Cheios[1:]=As[1:]/At[1:]
    return nos_vento_troco,Largura,Largura_total,Largura_media,At,Af,Ac,Acsup,As,Indice_Cheios
