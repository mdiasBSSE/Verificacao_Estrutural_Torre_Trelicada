import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
def angulo_entre_vetores(u, v, graus=True):
    """
    Calcula o ângulo entre dois vetores u e v.
    
    Parameters
    ----------
    u, v : array_like
        Vetores em qualquer dimensão (2D, 3D, ...).
    graus : bool
        Se True retorna o ângulo em graus, senão em radianos.
    
    Returns
    -------
    float
        Ângulo entre u e v.
    """
    u = np.array(u)
    v = np.array(v)
    
    # produto escalar
    dot = np.dot(u, v)
    
    # normas
    norm_u = np.linalg.norm(u)
    norm_v = np.linalg.norm(v)
    
    # cos(θ)
    cos_theta = dot / (norm_u * norm_v)
    
    # corrigir possíveis erros numéricos (valor ligeiramente fora de [-1,1])
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    
    theta = np.arccos(cos_theta)
    return np.degrees(theta) if graus else theta



def loaddistribution(y,y1,y2,load):
    ratio=(y-y1)/(y2-y1)
    load1=load*(1-ratio)
    load2=load*ratio
    return load1,load2

def encontra_vizinhos(vetor, valor):
    vetor = np.array(vetor)
    idx = np.searchsorted(vetor, valor)

    if idx == 0 or idx == len(vetor):
        return None
    # ajustar índice para o vetor original (decrescente)
    i_maior = len(vetor) - idx
    i_menor = i_maior - 1
    return np.array([i_menor, i_maior])

def find_closest_value_position(vector, reference_value):
    difference= abs(vector-reference_value)
    closest_position=np.argmin(difference)   
    return closest_position





def integral(x, y):
    return np.trapezoid(y, x)     
def agrupar_indices(d, tol=0.1):
    """
    Agrupa IDs de um dicionário {id: (x,y,z)} pela 3ª coluna (z) com tolerância,
    ignorando IDs ausentes (NaN). Ordena os grupos do z mais alto para o mais baixo.
    
    Parâmetros
    ----------
    d : dict
        {id: (x, y, z)}
    tol : float
        tolerância para agrupar por altura
    
    Retorna
    -------
    grupos : list of lists
        Cada lista contém IDs agrupados, do z mais alto para o mais baixo.
        IDs ausentes (NaN) não aparecem em nenhum grupo.
    """
    # Determinar intervalo completo de IDs
    id_min = min(d.keys())
    id_max = max(d.keys())
    n_ids = id_max - id_min + 1

    # Criar matriz com NaN para IDs ausentes
    mat = np.full((n_ids, 3), np.nan)
    ids_list = np.arange(id_min, id_max+1)
    for k, v in d.items():
        mat[k - id_min, :] = v

    # Seleciona apenas linhas válidas (não NaN) para o clustering
    valid_idx = ~np.isnan(mat[:, 2])
    valid_ids = ids_list[valid_idx]
    valid_coords = mat[valid_idx, :]

    # Clustering hierárquico pela 3ª coluna
    z = valid_coords[:, 2].reshape(-1, 1)
    Z = linkage(z, method='single', metric='euclidean')
    labels = fcluster(Z, t=tol, criterion='distance')

    grupos = []
    medias_z = []
    for lbl in np.unique(labels):
        idx = np.where(labels == lbl)[0]
        grupos.append(list(valid_ids[idx]))
        medias_z.append(valid_coords[idx, 2].mean())

    # Ordenar grupos do z mais alto para o mais baixo
    grupos_sorted = [g for _, g in sorted(zip(medias_z, grupos), key=lambda x: -x[0])]

    return grupos_sorted,medias_z