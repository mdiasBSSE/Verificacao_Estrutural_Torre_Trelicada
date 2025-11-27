import numpy as np
from openseespy import opensees as ops
from Utilities.Utilities import agrupar_indices
def Node_Filter_NCN(Nodes_Matrix):
    Node_Type=np.zeros(len(Nodes_Matrix)+1).astype(str)
    coords = {}
    coords_filtrados = {}
    for i in range(len(Nodes_Matrix)):
        tag = int(Nodes_Matrix[i,0])
        x, y, z = float(Nodes_Matrix[i,2])*10**-3,float(Nodes_Matrix[i,3])*10**-3,float(Nodes_Matrix[i,4])*10**-3
        ops.node(tag, x, y, z)
        coords[tag] = (x, y, z)
        Node_Type[tag]=Nodes_Matrix[i,5]
        if Node_Type[tag]!="Non-Connectable Node":
            coords_filtrados[tag]=(x, y, z)
        
    nos_grupo, z_grupo=agrupar_indices(coords_filtrados)

    nos_grupo=np.array(nos_grupo)
    z_grupo=np.sort(np.array(z_grupo))[::-1]

    return nos_grupo,z_grupo,Node_Type,coords


def Node_Filter_Montante_TrussType(Elements_Matrix):
    nos_montante=[]
    nEle=Elements_Matrix.shape[0]
    TrussType=np.zeros(nEle+1).astype(str)
    for i in range(len(Elements_Matrix)):
        eleTag = int(Elements_Matrix[i,0])
        ni = int(Elements_Matrix[i,2])
        nj = int(Elements_Matrix[i,3]) 
        TrussType[eleTag]=Elements_Matrix[i,4]
        if TrussType[eleTag]=="Manual Bar":
            nome_barra = str(Elements_Matrix[i,1])  
            if nome_barra[0].upper()=='E':
                TrussType[eleTag]="External Manual Bar"
            elif nome_barra[0].upper()=='I':
                TrussType[eleTag]="Internal Manual Bar"
        if TrussType[eleTag]=="Leg":
            nos_montante.append(ni)
            nos_montante.append(nj)
    nos_montante=np.array(nos_montante)
    nos_montante=np.unique(nos_montante)

    return nos_montante,TrussType

def Node_Constraint_func(Nodes_Matrix,Node_Type):
    i_Node_Type=Node_Type=="Tower Foundation"
    Node_Constraint=Nodes_Matrix[i_Node_Type[1:],0]
    if len(Node_Constraint)==3:
        torre="Triangular"
    else:
        torre="Quadrangular"
    for i in range(len(Node_Constraint)):
        ops.fix(int(Node_Constraint[i]), 1, 1, 1,0,0,0)  # fixa translação x,y,z
    return Node_Constraint,torre

def Node_Wind(Nodes_Matrix,Node_Type,torre):
    ###Carregar nós vento
    ###Assumindo que os nós vem organizados da base para o topo
    i_NodeType = np.logical_or(Node_Type == "Tower Foundation",
                            Node_Type == "Connection Shaft")
    nos_vento=Nodes_Matrix[i_NodeType[1:],0].astype(int)
    z_vento=Nodes_Matrix[i_NodeType[1:],4].astype(float)*10**-3
    nos_vento=np.array(nos_vento)
    z_vento=np.array(z_vento)
    if torre=="Triangular":
        nos_vento = nos_vento.reshape(-1, 3)
        z_vento = z_vento.reshape(-1, 3)
        divisor=3
    elif torre=="Quadrangular":
        nos_vento = nos_vento.reshape(-1, 4)  
        z_vento = z_vento.reshape(-1, 4)
        divisor=4
    z_vento=np.max(z_vento, axis=1)

    return nos_vento,z_vento,divisor

def Node_Class_Troco(z_vento,nEle):
    # Obter todos os nós do modelo
    todos_nos = ops.getNodeTags()

    # Coordenadas Z dos nós
    coord_z = {n: ops.nodeCoord(n)[2] for n in todos_nos}

    # Dicionário para guardar nós por classe
    nos_por_classe = {i: [] for i in range(len(z_vento)-1)}

    # Classificação
    for n, z in coord_z.items():
        for i in range(len(z_vento)-1):
            z_min = z_vento[i]
            z_max = z_vento[i+1]
            if (z >= 0.9*z_min) and (z <= z_max or (i == len(z_vento)-1 and z == z_max)):
                nos_por_classe[i].append(n)
    # Obter todos os elementos
    Elementos = ops.getEleTags()

    # Dicionário elemento → [nó inicial, nó final]
    nos_por_elemento = {ele: ops.eleNodes(ele) for ele in Elementos}
    Troco=np.zeros(nEle+1)
    # Mostrar resultados
    for ele, nos in nos_por_elemento.items():
        for i in range(len(nos_por_classe)):
            if nos[0] in nos_por_classe[i] and nos[1] in nos_por_classe[i]:
                Troco[ele]=i+1
                break
    
    Troco=np.max(Troco)+1-Troco
    
    return nos_por_classe,Troco,todos_nos, nos_por_elemento