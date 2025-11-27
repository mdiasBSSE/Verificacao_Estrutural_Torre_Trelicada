import numpy as np
from openseespy import opensees as ops


def distancia_entre_nos(no1, no2):
    """
    Retorna a distância entre dois nós no modelo OpenSeesPy.
    """
    coord1 = np.array(ops.nodeCoord(int(no1)))
    coord2 = np.array(ops.nodeCoord(int(no2)))
    return np.linalg.norm(coord2 - coord1)

def comprimentos_elementos_array():
    """
    Retorna dois arrays:
    - eleTags: IDs dos elementos
    - comprimentos: comprimentos correspondentes
    Funciona mesmo que o elemento tenha mais de dois nós.
    """
    eleTags = np.array(ops.getEleTags())
    comprimentos = np.zeros(len(eleTags)+1)

    for idx, ele in enumerate(eleTags):
        nos = ops.eleNodes(int(ele))
        coord1 = np.array(ops.nodeCoord(nos[0]))
        coord2 = np.array(ops.nodeCoord(nos[1]))
        comprimentos[idx+1] = np.linalg.norm(coord2 - coord1)
    
    return comprimentos

def no_mais_proximo(no_ref, lista_nos):
    coord_ref = np.array(ops.nodeCoord(int(no_ref)))
    distancias = []
    for n in lista_nos:
        coord_n = np.array(ops.nodeCoord(int(n)))
        dist = np.linalg.norm(coord_n - coord_ref)
        distancias.append((n, dist))
    # Ordenar pela distância
    distancias.sort(key=lambda x: x[1])
    return distancias[0]  # retorna (no, distância)

def elementos_do_no(no):
    """
    Retorna uma lista com os elementos que estão conectados ao nó dado.
    """
    elementos = ops.getEleTags()  # Todos os elementos do modelo
    elementos_no = []
    
    for ele in elementos:
        nos_ele = ops.eleNodes(ele)  # tupla com os nós do elemento
        if no in nos_ele:
            elementos_no.append(ele)
    
    return elementos_no


def elementos_do_no_filt(noTag,TrussType, tipo):
    """
    Retorna os elementos do modelo OpenSeesPy conectados ao nó dado
    e que sejam do tipo especificado (ex.: 'Leg').
    """
    elems = []
    for ele in ops.getEleTags():
        if noTag in ops.eleNodes(ele):
            if TrussType[ele] == tipo:
                elems.append(ele)
    return elems

def elementos_entre(z_min, z_max):
    """
    Retorna todos os elementos cujos nós de extremidade
    estão completamente entre z_min e z_max.
    """
    elems = []
    for eleTag in ops.getEleTags():
        nodes = ops.eleNodes(eleTag)
        coords = [ops.nodeCoord(nd)[2] for nd in nodes]  # pega z de cada nó
        if min(coords) >= z_min and max(coords) <= z_max:
            elems.append(eleTag)
    return np.array(elems)

def elementos_acima(z_min):
    """
    Retorna lista de tuplos (eleTag, perc_acima)
    onde perc_acima é a percentagem do comprimento do elemento
    que está acima de z_min.
    """
    elems = []
    for eleTag in ops.getEleTags():
        nodes = ops.eleNodes(eleTag)
        coords = [ops.nodeCoord(nd) for nd in nodes]
        z = [c[2] for c in coords]

        z1, z2 = z
        L = np.linalg.norm(np.array(coords[1]) - np.array(coords[0]))  # comprimento total

        if max(z1, z2) <= z_min:
            # totalmente abaixo
            continue
        elif min(z1, z2) >= z_min:
            # totalmente acima
            continue
        else:
            # cruza o plano em z_min → calcular fração
            dz_total = abs(z2 - z1)
            dz_acima = max(z2, z1) - z_min
            perc = dz_acima / dz_total

        elems.append((eleTag, perc))

    return np.array(elems)
def elementos_abaixo(z_max):
    """
    Retorna lista de tuplos (eleTag, perc_abaixo)
    onde perc_abaixo é a percentagem do comprimento do elemento
    que está abaixo de z_max.
    """
    elems = []
    for eleTag in ops.getEleTags():
        nodes = ops.eleNodes(eleTag)
        coords = [ops.nodeCoord(nd) for nd in nodes]
        z = [c[2] for c in coords]

        z1, z2 = z
        L = np.linalg.norm(np.array(coords[1]) - np.array(coords[0]))  # comprimento total

        if min(z1, z2) >= z_max:
            # totalmente acima
            continue
        elif max(z1, z2) <= z_max:
            # totalmente abaixo
            continue
        else:
            # cruza o plano em z_max → calcular fração
            dz_total = abs(z2 - z1)
            dz_abaixo = z_max - min(z2, z1)
            perc = dz_abaixo / dz_total

        elems.append((eleTag, perc))

    return np.array(elems)


def get_nodes_same_elevation(node_ref, tol=1e-6):
    """
    Retorna todos os nós com a mesma cota (z) do nó de referência.
    
    Parameters
    ----------
    node_ref : int
        Tag do nó de referência
    tol : float, opcional
        Tolerância para comparação de coordenadas (padrão: 1e-6)
    
    Returns
    -------
    nodes_same_z : list
        Lista com os nós na mesma cota
    """
    # Coordenadas do nó de referência
    x_ref, y_ref, z_ref = ops.nodeCoord(node_ref)
    
    # Lista de todos os nós do modelo
    all_nodes = ops.getNodeTags()
    
    # Filtrar por mesma cota
    nodes_same_z = []
    for n in all_nodes:
        x, y, z = ops.nodeCoord(n)
        if abs(z - z_ref) <= tol:
            nodes_same_z.append(n)
    nodes_same_z=np.array(nodes_same_z)
    return nodes_same_z