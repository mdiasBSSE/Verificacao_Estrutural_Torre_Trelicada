import numpy as np
from openseespy import opensees as ops
from Output.Output_Lattice_Reduced import Buckling_Output
from Output.Output_Lattice_Expanded import Buckling_Exp_Output
from Utilities.Utilities_ops import elementos_do_no_filt
def norm(v):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    return v / n if n != 0 else v

def rotation_matrix_from_vectors(a, b):
    a = norm(a)
    b = norm(b)
    v = np.cross(a, b)
    c = np.dot(a, b)
    
    if np.isclose(c, 1):
        return np.eye(3)
    if np.isclose(c, -1):
        # 180º rotation → pick any orthogonal vector
        axis = norm(np.array([1,0,0]) if abs(a[0]) < 0.9 else np.array([0,1,0]))
        return -np.eye(3) + 2*np.outer(axis, axis)

    s = np.linalg.norm(v)
    K = np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

    return np.eye(3) + K + K @ K * ((1 - c) / (s**2))
def Buckling_Lenght(Ele, Length, TrussType):
    import numpy as np
    nos_elemento = ops.eleNodes(Ele)
    n1, n2 = nos_elemento

    # ==========================
    # 1) Direção principal da barra
    # ==========================
    c1 = np.array(ops.nodeCoord(n1), dtype=float)
    c2 = np.array(ops.nodeCoord(n2), dtype=float)
    v_main = c2 - c1
    L = np.linalg.norm(v_main)
    v_main_u = v_main / L

    # ==========================
    # 2) Construção do referencial local
    #     R_global_to_local @ v_main = (0, 0, 1)
    # ==========================
    R = rotation_matrix_from_vectors(v_main_u, np.array([0, 0, 1]))

    elementos = ops.getEleTags()

    Blockedy  = np.empty(2, dtype='U64')
    Blockedv  = np.empty(2, dtype='U64')
    Bracing_A = np.empty(2, dtype='U64')

    for idx, no in enumerate(nos_elemento):

        Dir_global = []

        # ==========================
        # 3) Obter direções das barras ligadas ao nó
        # ==========================
        for ele in elementos:
            nA, nB = ops.eleNodes(ele)
            #if no in (nA, nB) and TrussType[ele] != "Leg":
            if no in (nA, nB):
                cA = np.array(ops.nodeCoord(nA), dtype=float)
                cB = np.array(ops.nodeCoord(nB), dtype=float)
                v = cB - cA

                if np.linalg.norm(v) > 0:
                    Dir_global.append(v / np.linalg.norm(v))

        Dir_global = np.array(Dir_global)

        # ==========================
        # 4) Transformar para coordenadas locais
        #    v_local = R * v_global
        # ==========================
        # print(Ele)
        # print(R)
        # print(Dir_global)
        Dir_local = (R @ Dir_global.T).T
        

        # ==========================
        # 5) CLASSIFICAÇÕES NO REFERENCIAL LOCAL
        # ==========================
        has_x = np.any(np.abs(Dir_local[:,0]) > 1e-4)
        has_y = np.any(np.abs(Dir_local[:,1]) > 1e-4)
        has_z = np.any(np.abs(Dir_local[:,2]) > 1e-4)

        # ---- Bloqueio local no eixo Y ----
        if (has_x and has_y) or ops.nodeCoord(no)[2] == 0:
            Blockedy[idx] = "True"
        else:
            Blockedy[idx] = "False"

        # ---- Simetria local ----
        if has_x != has_y:
            Bracing_A[idx] = "Unsymmetrical"
        else:
            Bracing_A[idx] = "Symmetrical"

        # ---- Bloqueio após rotação adicional (45° dentro do local) ----
        theta = np.radians(45)
        Rz = np.array([
            [np.cos(theta), -np.sin(theta), 0],
            [np.sin(theta),  np.cos(theta), 0],
            [0, 0, 1]
        ])

        Dir_v = (Rz @ Dir_local.T).T
        
        has_x2 = np.any(np.abs(Dir_v[:,0]) > 1e-4)
        has_y2 = np.any(np.abs(Dir_v[:,1]) > 1e-4)

        if (has_x2 and has_y2) or ops.nodeCoord(no)[2] == 0:
            Blockedv[idx] = "True"
        else:
            Blockedv[idx] = "False"

    # ==========================
    # 6) RESULTADOS FINAIS
    # ==========================
    L_buckling_y = Length if np.all(Blockedy == "True") else Length * 2
    L_buckling_v = Length if np.all(Blockedv == "True") else Length * 2
    Bracing = "Symmetrical" if np.all(Bracing_A == "Symmetrical") else "Unsymmetrical"
    return L_buckling_y, L_buckling_v, Bracing

# def Buckling_Lenght_old(Ele,Lenght,TrussType):
#     elementos = ops.getEleTags()  # Todos os elementos do modelo
#     elementos_no = []
#     Dir_no = []
#     Blockedy=np.empty(2).astype(str)
#     Bracing=np.empty(2).astype(str)
#     Blockedv=np.empty(2).astype(str)
#     nos = ops.eleNodes(Ele) 
#     i_buckling=0
#     for no in nos:
#         elementos_no = []
#         Dir_no = []
#         for ele in elementos:
#             nos_ele = ops.eleNodes(ele)  # tupla com os nós do elemento
#             if no in nos_ele:
#                 if TrussType[ele]!="Leg":
#                     elementos_no.append(ele)
#                     coords1 = np.array(ops.nodeCoord(nos_ele[0]))
#                     coords2 = np.array(ops.nodeCoord(nos_ele[1]))
#                     Dir=(coords2-coords1)/np.linalg.norm(coords2-coords1)
#                     Dir_no.append(Dir)
#         Dir_no=np.array(Dir_no)
#         cols_nonzero = np.any(Dir_no != 0, axis=0)

        
#         if cols_nonzero[0] and cols_nonzero[1] or ops.nodeCoord(no)[2]==0:
#             Blockedy[i_buckling]="True"
#         else:
#             Blockedy[i_buckling]="False"
        
#         if cols_nonzero[0]!=cols_nonzero[1]:
#             Bracing[i_buckling]="Unsymmetrical"
#         else:
#             Bracing[i_buckling]="Symmetrical"
        
#         for j in range(Dir_no.shape[0]):
#             theta_rot=np.radians(45)
#             TransRot=np.array([[np.cos(theta_rot),-np.sin(theta_rot),0],[np.sin(theta_rot),np.cos(theta_rot),0],[0,0,1]])
#             Dir_no[j,:]=TransRot@Dir_no[j,:]
#         cols_nonzero = np.any(Dir_no != 0, axis=0)
#         if cols_nonzero[0] and cols_nonzero[1]or ops.nodeCoord(no)[2]==0:
#             Blockedv[i_buckling]="True"
#         else:
#             Blockedv[i_buckling]="False"
#         i_buckling=i_buckling+1
#     if np.all(Blockedy=="True"):
#         L_buckling_y=Lenght
#     else:
#         L_buckling_y=Lenght*2
#     if np.all(Blockedv=="True"):
#         L_buckling_v=Lenght
#     else:
#         L_buckling_v=Lenght*2
#     if np.all(Bracing=="Symmetrical"):
#         Bracing="Symmetrical"
#     else:
#         Bracing="Unsymmetrical"

    
    
    
#     return L_buckling_y,L_buckling_v,Bracing


def Imp_factor_buckling_curve(Buckling_Curve):
    if Buckling_Curve=="a0":
        alpha_Enc=0.13
    elif Buckling_Curve=="a":
        alpha_Enc=0.21
    elif Buckling_Curve=="b":
        alpha_Enc=0.34   
    elif Buckling_Curve=="c":
        alpha_Enc=0.49
    elif Buckling_Curve=="d":
        alpha_Enc=0.76
    return alpha_Enc
# def Block_direction_old (Ele,nos_montante):
#     elementos = ops.getEleTags()  # Todos os elementos do modelo
#     elementos_no = []
#     Dir_no = []
#     Block_vert_down=False
#     Block_vert_up=False
#     Block_inside=False
#     nos = ops.eleNodes(Ele)
#     if nos[0] in nos_montante and nos[1] not in nos_montante: 
#         j=1
#         Flag_Enc=False
#     elif nos[1] in nos_montante and nos[0] not in nos_montante: 
#         j=0
#         Flag_Enc=False
#     elif nos[0] in nos_montante and nos[1] in nos_montante: #Se os dois nós tiverem nos montantes
#         Block_vert_up=True
#         Block_vert_down=True
#         Block_inside=True
#         Flag_Enc=True #Se ambos os nos tiverem no montante significa que está bloqueado
#         return Block_vert_up,Block_vert_down,Block_inside,Flag_Enc
#     else:
#         j=0 #Se nenhum estiver nos montantes é uma interna manual
#         Flag_Enc=False
#     elementos_no = []
#     Dir_no = []
#     for ele in elementos:
#         nos_ele = ops.eleNodes(ele)  # tupla com os nós do elemento
#         if nos[j] in nos_ele: #Analisa apenas o no que não está no montante
#             elementos_no.append(ele)
#             coords1 = np.array(ops.nodeCoord(int(nos[j])))
#             if nos[j]==nos_ele[0]:
#                 coords2 = np.array(ops.nodeCoord(nos_ele[1]))
#             else:
#                 coords2 = np.array(ops.nodeCoord(nos_ele[0]))
#             Dir=(coords2-coords1)/np.linalg.norm(coords2-coords1) ##Vetor a sair do nó do elemento em analise
            
#             Dir_no.append(Dir)
#     Dir_no=np.array(Dir_no)
#     cols_vert_up = np.any(Dir_no > 0, axis=0)
#     if cols_vert_up[2]: ###Se tiver algum bloqueio na vertical (Cima)
#         Block_vert_up=True
    
#     cols_vert_down=  np.any(Dir_no < 0, axis=0)
#     if cols_vert_down[2]: ###Se tiver algum bloqueio na vertical (Baixo)
#         Block_vert_down=True
    
#     Block_inside = np.any((Dir_no[:,2]==0) & (Dir_no[:,0]!=0) & (Dir_no[:,1]!=0)) ###Block_inside é True se houver um vetor que tenha z=0, x!=0 e y!=0
#     return Block_vert_up,Block_vert_down,Block_inside,Flag_Enc

def Block_direction (Ele,nos_montante,TrussType,Troco):
    elementos = ops.getEleTags()  # Todos os elementos do modelo
    elementos_no = []
    Dir_no = []
    Block_vert_down=False
    Block_vert_up=False
    Block_inside=False
    nos = ops.eleNodes(Ele)
    if nos[0] in nos_montante and nos[1] not in nos_montante: 
        j=1 #Nó fora do montante
        k=0 #Nó do montante 
        Flag_Enc=False
    elif nos[1] in nos_montante and nos[0] not in nos_montante: 
        j=0
        k=1
        Flag_Enc=False
    elif nos[0] in nos_montante and nos[1] in nos_montante: #Se os dois nós tiverem nos montantes
        Block_vert_up=True
        Block_vert_down=True
        Block_inside=True
        Flag_Enc=True #Se ambos os nos tiverem no montante significa que está bloqueado
        return Block_vert_up,Block_vert_down,Block_inside,Flag_Enc
    else:
        j=0 #Se nenhum estiver nos montantes é uma interna manual e portanto o comprimento de encurvadura é L
        k=0
        Block_vert_up=True
        Block_vert_down=True
        Block_inside=True
        Flag_Enc=True #Se ambos os nos tiverem no montante significa que está bloqueado
        return Block_vert_up,Block_vert_down,Block_inside,Flag_Enc
    elementos_no = []
    Dir_no = []
    if TrussType[Ele]=="Diagonal": #Trabalhar na seleção do montante
        if Troco[Ele]==Troco[np.max(elementos_do_no_filt(nos[k],TrussType,"Leg"))]:
            ELE_Montante=np.max(elementos_do_no_filt(nos[k],TrussType,"Leg"))
        else:
            ELE_Montante=np.min(elementos_do_no_filt(nos[k],TrussType,"Leg"))
    else:
        ELE_Montante=np.min(elementos_do_no_filt(nos[k],TrussType,"Leg"))
    nos_ele_montante = ops.eleNodes(int(ELE_Montante))
    coords1_montante = np.array(ops.nodeCoord(int(nos[k])))
    if nos[k]==nos_ele_montante[0]:
        coords2_montante = np.array(ops.nodeCoord(nos_ele_montante[1]))
    else:
        coords2_montante = np.array(ops.nodeCoord(nos_ele_montante[0]))
    Dir_montante=(coords2_montante-coords1_montante)/np.linalg.norm(coords2_montante-coords1_montante)
    coords2_Ele = np.array(ops.nodeCoord(int(nos[j])))
    if nos[j]==nos[0]:
        coords1_Ele = np.array(ops.nodeCoord(nos[1]))
    else:
        coords1_Ele = np.array(ops.nodeCoord(nos[0]))
    Dir_ELE=(coords2_Ele-coords1_Ele)/np.linalg.norm(coords2_Ele-coords1_Ele)
    Dir_Z_Local=np.cross(Dir_montante, Dir_ELE)
    Dir_Z_Local=Dir_Z_Local/np.linalg.norm(Dir_Z_Local)
    Dir_X_Local=Dir_ELE
    Dir_Y_Local=np.cross(Dir_X_Local, Dir_Z_Local)
    # print("Dir_Montante: ",Dir_montante)
    # print("X: ", Dir_X_Local)
    # print("Y: ", Dir_Y_Local)
    # print("Z: ", Dir_Z_Local)
    R = np.vstack((Dir_X_Local, Dir_Y_Local, Dir_Z_Local)) 
    for ele in elementos:
        nos_ele = ops.eleNodes(ele)  # tupla com os nós do elemento
        if nos[j] in nos_ele: #Analisa apenas o no que não está no montante
            elementos_no.append(ele)
            coords1 = np.array(ops.nodeCoord(int(nos[j])))
            if nos[j]==nos_ele[0]:
                coords2 = np.array(ops.nodeCoord(nos_ele[1]))
            else:
                coords2 = np.array(ops.nodeCoord(nos_ele[0]))
            Dir=(coords2-coords1)/np.linalg.norm(coords2-coords1) ##Vetor a sair do nó do elemento em analise
            Dir_Local=R @ Dir
            Dir_no.append(Dir_Local)
    
    Dir_no=np.array(Dir_no)
    cols_vert_up = np.any(Dir_no > 0, axis=0)
    if cols_vert_up[1]: ###Se tiver algum bloqueio na vertical (Cima) Y positivo
        Block_vert_up=True
    
    cols_vert_down=  np.any(Dir_no < 0, axis=0)
    if cols_vert_down[1]: ###Se tiver algum bloqueio na vertical (Baixo) Y negativo
        Block_vert_down=True
    Block_inside = np.any(np.abs(Dir_no[:, 2]) > 10**-6)###Block_inside é True se houver um vetor que tenha z!=0
    return Block_vert_up,Block_vert_down,Block_inside,Flag_Enc

def Horizontal_Buckling(Ele,nos_montante,Type_Con,N_bolt_con,Exp_vento,TrussType,Troco):
    Block_vert_up,Block_vert_down,Block_inside,Flag_Enc=Block_direction(Ele,nos_montante,TrussType,Troco)
    eta=1.0
    if Flag_Enc: #Se os dois nós tiverem nos montantes
        Multy=1
        Multv=1
        if Type_Con=="No" and N_bolt_con==1 and Exp_vento=="Flat":
            eta=0.8
        else:
            eta=1.0
    else: #Só entra aqui se não for uma barra de montante a montante
        if Block_inside and Block_vert_up and Block_vert_down: #Barras em X completo com travamento interno
            Multy=1
            Multv=1
            if Type_Con=="No" and N_bolt_con==1 and Exp_vento=="Flat":
                eta=0.9
            else:
                eta=1.0
        elif Block_vert_up and Block_vert_down:  #Barras em X completo sem travamento interno
            Multy=2
            Multv=1
            if Type_Con=="No" and N_bolt_con==1 and Exp_vento=="Flat":
                eta=0.9
            else:
                eta=1.0     
        elif (Block_vert_down and Block_inside) or (Block_vert_up and Block_inside) : #Barras em meio X  com travamento interno
            Multy=1
            Multv=1
            if Type_Con=="No" and N_bolt_con==1 and Exp_vento=="Flat":
                eta=0.8
            else:
                eta=1.0   
        elif (Block_vert_down) or (Block_vert_up) or (Block_inside):
            Multy=2
            Multv=1
            if Type_Con=="No" and N_bolt_con==1 and Exp_vento=="Flat":
                eta=0.8
            else:
                eta=1.0  
    return Multy,Multv,eta
def Diagonal_Buckling(Ele,nos_montante,Type_Con,N_bolt_con,Exp_vento,TrussType,Troco):
    Block_vert_up,Block_vert_down,Block_inside,Flag_Enc=Block_direction(Ele,nos_montante,TrussType,Troco)
    
    eta=1.0
    if Flag_Enc:
        Multy=1
        Multv=1
        if Type_Con=="No" and N_bolt_con==1 and Exp_vento=="Flat":
            eta=0.8
        else:
            eta=1.0
    else:
        if Block_inside and Block_vert_up and Block_vert_down:
            Multy=1
            Multv=1
            
            if Type_Con=="No" and N_bolt_con==1 and Exp_vento=="Flat":
                eta=0.9
            else:
                eta=1.0
        elif Block_vert_up and Block_vert_down:
            Multy=1
            Multv=1
            if Type_Con=="No" and N_bolt_con==1 and Exp_vento=="Flat":
                eta=0.9
            else:
                eta=1.0     
        elif (Block_vert_down and Block_inside) or (Block_vert_up and Block_inside) :
            Multy=1
            Multv=1
            if Type_Con=="No" and N_bolt_con==1 and Exp_vento=="Flat":
                eta=0.8
            else:
                eta=1.0   
        elif (Block_vert_down) or (Block_vert_up) or (Block_inside):
            Multy=1
            Multv=1
            if Type_Con=="No" and N_bolt_con==1 and Exp_vento=="Flat":
                eta=0.8
            else:
                eta=1.0  

    return Multy,Multv,eta


def Buckling_Function(Pais,nEle,Inertiav,Inertiau,Inertiay,Inertiaz,Area,Comprimento_Barra,Sigma,Rigidez,nos_vento,Shape,Travamento,Top_bars,Middle_bars,ExpVento,Dim,Esp,TrussType,nos_montante,Conection_Ele,N_Bolt_Ele,GammaM1,F_Compressao,Size_Profile,Elements_Matrix,NTracRd,Troco):
    if Pais in ["Portugal","Portugal-RSA"]:
        Ratio_Enc,Enc_Out_csv,Enc_Out_Exp=Buckling_Function_Portugal(nEle,Inertiav,Inertiau,Inertiay,Inertiaz,Area,Comprimento_Barra,Sigma,Rigidez,nos_vento,Shape,Travamento,Top_bars,Middle_bars,ExpVento,Dim,Esp,TrussType,nos_montante,Conection_Ele,N_Bolt_Ele,GammaM1,F_Compressao,Size_Profile,Elements_Matrix,NTracRd,Troco)
        return Ratio_Enc,Enc_Out_csv,Enc_Out_Exp
    elif Pais in ["France"]:
        Ratio_Enc,Enc_Out_csv,Enc_Out_Exp=Buckling_Function_France(nEle,Inertiav,Inertiau,Inertiay,Inertiaz,Area,Comprimento_Barra,Sigma,Rigidez,nos_vento,Shape,Travamento,Top_bars,Middle_bars,ExpVento,Dim,Esp,TrussType,nos_montante,Conection_Ele,N_Bolt_Ele,GammaM1,F_Compressao,Size_Profile,Elements_Matrix,NTracRd,Troco)
        return Ratio_Enc,Enc_Out_csv,Enc_Out_Exp        
    elif Pais in ["Spain"]:
        Ratio_Enc,Enc_Out_csv,Enc_Out_Exp=Buckling_Function_France(nEle,Inertiav,Inertiau,Inertiay,Inertiaz,Area,Comprimento_Barra,Sigma,Rigidez,nos_vento,Shape,Travamento,Top_bars,Middle_bars,ExpVento,Dim,Esp,TrussType,nos_montante,Conection_Ele,N_Bolt_Ele,GammaM1,F_Compressao,Size_Profile,Elements_Matrix,NTracRd,Troco)
        return Ratio_Enc,Enc_Out_csv,Enc_Out_Exp   
def Buckling_Function_Portugal(nEle,Inertiav,Inertiau,Inertiay,Inertiaz,Area,Comprimento_Barra,Sigma,Rigidez,nos_vento,Shape,Travamento,Top_bars,Middle_bars,ExpVento,Dim,Esp,TrussType,nos_montante,Conection_Ele,N_Bolt_Ele,GammaM1,F_Compressao,Size_Profile,Elements_Matrix,NTracRd,Troco):
    Buckling_Curve=np.zeros(nEle+1).astype(str)
    Class_Enc=np.zeros(nEle+1)
    Check_Enc=np.zeros(nEle+1).astype(str)
    Comprimento_Enc=np.zeros(nEle+1)
    Lambda_Lim=np.zeros(nEle+1)

    Lambda_v=np.zeros(nEle+1)
    rho_Enc=np.ones(nEle+1)
    Lambda_u=np.zeros(nEle+1)
    Lambda_y=np.zeros(nEle+1)
    Lambda_Barra_p=np.zeros(nEle+1)
    Lambda_z=np.zeros(nEle+1)
    Lambda_Ratio=np.zeros(nEle+1) 
    Lambda_Ratio_v=np.zeros(nEle+1) 
    Lambda_Ratio_u=np.zeros(nEle+1) 
    Lambda_Ratio_y=np.zeros(nEle+1) 
    Lambda_Ratio_z=np.zeros(nEle+1) 
    Lambda_Enc_eff=np.zeros(nEle+1) 
    Lambda_Enc_eff_v=np.zeros(nEle+1) 
    Lambda_Enc_eff_u=np.zeros(nEle+1) 
    Lambda_Enc_eff_y=np.zeros(nEle+1)
    Lambda_Enc_eff_z=np.zeros(nEle+1) 
    Lim_Class=np.zeros(nEle+1) 
    Ratio_Class=np.zeros(nEle+1) 
    Lambda_Check=np.zeros(nEle+1).astype(str)
    alpha_Enc=np.zeros(nEle+1).astype(str)
    k_Enc=np.ones(nEle+1)
    Red=np.ones(nEle+1)
    Bracing=np.ones(nEle+1).astype(str)
    eta_Enc=np.ones(nEle+1)
    Aeff=np.ones(nEle+1)
    ksigma=np.ones(nEle+1)
    k_v=np.ones(nEle+1)
    k_u=np.ones(nEle+1)
    k_y=np.ones(nEle+1)
    k_z=np.ones(nEle+1)
    iv=np.ones(nEle+1)
    iu=np.ones(nEle+1)
    iy=np.ones(nEle+1)
    iz=np.ones(nEle+1)
    iv[1:]=np.sqrt(Inertiav[1:]/Area[1:])
    iu[1:]=np.sqrt(Inertiau[1:]/Area[1:])
    iy[1:]=np.sqrt(Inertiay[1:]/Area[1:])
    iz[1:]=np.sqrt(Inertiaz[1:]/Area[1:])
    Lambda=Comprimento_Barra/np.min(np.vstack([iv, iu, iz, iy]), axis=0)
    epsilon=np.ones(nEle+1)
    Lambda_1=np.ones(nEle+1)
    epsilon[1:]=np.sqrt(235/(Sigma[1:]*10**-6))
    Lambda_1[1:]=np.pi*np.sqrt(Rigidez[1:]/Sigma[1:])
    #####Possivelemente dá para vetorizar isto #####
    nos_vento_inv=nos_vento[::-1,:]
    for i in range(1,len(Class_Enc)):
        if (Shape[i]=="X-Bracing") or (Shape[i]=="Panel Bracing") or ((Shape[i]=="Diamond Bracing") and (Top_bars[i]=="yes") and (Middle_bars[i]=="yes")) or ((Shape[i]=="Diamond Bracing 2") and (Top_bars[i]=="yes") and (Middle_bars[i]=="yes")) or ((Shape[i]=="Alternating Diagonal") and (Top_bars[i]=="yes") and (Middle_bars[i]=="yes")):
            Travamento[i]="Case b"
        elif ((Shape[i]=="Diamond Bracing") and ((Top_bars[i]!="yes") or (Middle_bars[i]!="yes"))) or ((Shape[i]=="Diamond Bracing 2") and ((Top_bars[i]!="yes") or (Middle_bars[i]!="yes"))):
            Travamento[i]="Case a"
        elif ((Shape[i]=="Alternating Diagonal") and ((Top_bars[i]!="yes") or (Middle_bars[i]!="yes"))):
            Travamento[i]="Case d"
        if ExpVento[i]=="Circ":
            Ratio_Class[i]=Dim[i]/Esp[i]
            Lim_Class[i]=90*epsilon[i]**2
            if Dim[i]/Esp[i]<=90*epsilon[i]**2:
                Class_Enc[i]=3
            else:
                Class_Enc[i]=4
            if TrussType[i]!="Leg":
                
                if TrussType[i] in ["Horizontal Bar", "External Manual Bar"]:
                    Multy,Multv,eta_Enc[i]=Horizontal_Buckling(int(i),nos_montante,Conection_Ele[i],N_Bolt_Ele[i],ExpVento[i],TrussType,Troco)
                    Comprimento_Enc[i]=np.max([Comprimento_Barra[i]*Multv,Comprimento_Barra[i]*Multy])
                
                elif TrussType[i]=="Internal Manual Bar" and ops.eleNodes(int(i))[0] in nos_montante and ops.eleNodes(int(i))[1] in nos_montante:
                    Comprimento_Enc[i]=Comprimento_Barra[i]*1
                else:
                    Multy,Multv,eta_Enc[i]=Diagonal_Buckling(int(i),nos_montante,Conection_Ele[i],N_Bolt_Ele[i],ExpVento[i],TrussType,Troco)
                    Comprimento_Enc[i]=np.max([Comprimento_Barra[i]*Multv,Comprimento_Barra[i]*Multy])
                if Conection_Ele[i]=="Yes":
                    k_Enc[i]=0.7
                elif Conection_Ele[i]=="No":
                    k_Enc[i]=0.95
                Lambda_y[i]=Comprimento_Enc[i]/iy[i]

                Lambda[i]=np.max([Lambda_v[i],Lambda_u[i],Lambda_y[i],Lambda_z[i]])
                Lambda_Lim[i]=180
                Lambda_Ratio[i]=Lambda[i]/Lambda_1[i]
                # if Lambda[i]<=Lambda_Lim[i]:
                #     Lambda_Check[i]="OK"
                # else:
                #     Lambda_Check[i]="KO"
            elif TrussType[i]=="Leg":
                L_enc_y,L_enc_v,Bracing[i]=Buckling_Lenght(int(i),Comprimento_Barra[i],TrussType)
                k_Enc[i]=1
                Comprimento_Enc[i]=L_enc_y

                Lambda_y[i]=Comprimento_Enc[i]/iy[i]

                Lambda[i]=np.max([Lambda_v[i],Lambda_u[i],Lambda_y[i],Lambda_z[i]])
                Lambda_Lim[i]=120
                Lambda_Ratio[i]=Lambda[i]/Lambda_1[i]
                # if Lambda[i]<=Lambda_Lim[i]:
                #     Lambda_Check[i]="OK"
                # else:
                #     Lambda_Check[i]="KO"
            Buckling_Curve[i]="c"
            alpha_Enc[i]=Imp_factor_buckling_curve(Buckling_Curve[i])
            Red[i]=1
            Aeff[i]=Area[i]
            Lambda_Enc_eff[i]=k_Enc[i]*Lambda_Ratio[i]
        elif ExpVento[i]=="Flat":
            Ratio_Class[i]=(Dim[i]-2*Esp[i])/Esp[i]
            Lim_Class[i]=15*epsilon[i]
            if ((Dim[i]-2*Esp[i])/Esp[i]<=15*epsilon[i]):
                Class_Enc[i]=3
            else:
                Class_Enc[i]=4   
            ksigma[i]=0.43

            if TrussType[i]!="Leg":
                if TrussType[i] in ["Horizontal Bar", "External Manual Bar"]:
                    Multy,Multv,eta_Enc[i]=Horizontal_Buckling(int(i),nos_montante,Conection_Ele[i],N_Bolt_Ele[i],ExpVento[i],TrussType,Troco)
                    L2_Enc=Comprimento_Barra[i]*Multv
                    L1_Enc=Comprimento_Barra[i]*Multy
                elif TrussType[i]=="Internal Manual Bar" and ops.eleNodes(int(i))[0] in nos_montante and ops.eleNodes(int(i))[1] in nos_montante:
                    L2_Enc=Comprimento_Barra[i]*0.5
                    L1_Enc=Comprimento_Barra[i]*1
                    if Conection_Ele[i]=="No" and N_Bolt_Ele[i]==1:
                        k_Enc[i]=0.9
                    else:
                        k_Enc[i]=1
                else:
                    Multy,Multv,eta_Enc[i]=Diagonal_Buckling(int(i),nos_montante,Conection_Ele[i],N_Bolt_Ele[i],ExpVento[i],TrussType,Troco)
                    L2_Enc=Comprimento_Barra[i]*Multv
                    L1_Enc=Comprimento_Barra[i]*Multy
                Lambda_Barra_p[i]=((Dim[i]-2*Esp[i])/Esp[i])/(28.4*epsilon[i]*np.sqrt(ksigma[i]))
                if Class_Enc[i]==4:
                    if Lambda_Barra_p[i]<=0.748:
                        rho_Enc[i]=1
                    else:
                        rho_Enc[i]=(Lambda_Barra_p[i]-0.188)/(Lambda_Barra_p[i]**2)
                        if rho_Enc[i]>=1:
                            rho_Enc[i]=1
                    Lambda_v[i]=L2_Enc/iv[i]*np.sqrt(rho_Enc[i]) ### 1993-1-1 altera o lambda em casos de classe 4
                    Lambda_u[i]=L2_Enc/iu[i]*np.sqrt(rho_Enc[i])
                    Lambda_y[i]=L1_Enc/iy[i]*np.sqrt(rho_Enc[i])
                    Lambda_z[i]=L1_Enc/iz[i]*np.sqrt(rho_Enc[i])
                    Lambda[i]=np.max([Lambda_v[i],Lambda_u[i],Lambda_y[i],Lambda_z[i]])
                    i_Enc_crit=np.argmax([Lambda_v[i],Lambda_u[i],Lambda_y[i],Lambda_z[i]])      
                else:
                    Lambda_v[i]=L2_Enc/iv[i]
                    Lambda_u[i]=L2_Enc/iu[i]
                    Lambda_y[i]=L1_Enc/iy[i]
                    Lambda_z[i]=L1_Enc/iz[i]
                    Lambda[i]=np.max([Lambda_v[i],Lambda_u[i],Lambda_y[i],Lambda_z[i]])
                    i_Enc_crit=np.argmax([Lambda_v[i],Lambda_u[i],Lambda_y[i],Lambda_z[i]])             
                Lambda_Lim[i]=180
                Lambda_Ratio[i]=Lambda[i]/Lambda_1[i]
                Lambda_Ratio_v[i]=Lambda_v[i]/Lambda_1[i]
                Lambda_Ratio_u[i]=Lambda_u[i]/Lambda_1[i]
                Lambda_Ratio_y[i]=Lambda_y[i]/Lambda_1[i]
                Lambda_Ratio_z[i]=Lambda_z[i]/Lambda_1[i]

                if (Conection_Ele[i]=="No") & (N_Bolt_Ele[i]==1) & (eta_Enc[i]==0.8):
                    k_v[i]=0.7+0.35/Lambda_Ratio_v[i]
                    k_y[i]=0.7+0.58/Lambda_Ratio_y[i]
                    k_z[i]=0.7+0.58/Lambda_Ratio_z[i]
                    k_u[i]=0.0001
                else:
                    k_v[i]=0.7+0.35/Lambda_Ratio_v[i]
                    k_y[i]=0.7+0.4/Lambda_Ratio_y[i]
                    k_z[i]=0.7+0.4/Lambda_Ratio_z[i]
                    k_u[i]=0.0001
                Lambda_Enc_eff_v[i]=k_v[i]*Lambda_Ratio_v[i]
                Lambda_Enc_eff_u[i]=k_u[i]*Lambda_Ratio_u[i]
                Lambda_Enc_eff_y[i]=k_y[i]*Lambda_Ratio_y[i]
                Lambda_Enc_eff_z[i]=k_z[i]*Lambda_Ratio_z[i]
                
                Lambda_Enc_eff[i]=np.max([Lambda_Enc_eff_v[i],Lambda_Enc_eff_u[i],Lambda_Enc_eff_y[i],Lambda_Enc_eff_z[i]])
                i_Enc_crit=np.argmax([Lambda_Enc_eff_v[i],Lambda_Enc_eff_u[i],Lambda_Enc_eff_y[i],Lambda_Enc_eff_z[i]]) 
                k_Enc_Aux=np.array([k_v[i],k_u[i],k_y[i],k_z[i]])
                k_Enc[i]=k_Enc_Aux[i_Enc_crit]
                if i_Enc_crit in [0,1]:
                    Comprimento_Enc[i]=L2_Enc
                else:
                    Comprimento_Enc[i]=L1_Enc
            elif TrussType[i]=="Leg":
                L_enc_y,L_enc_v,Bracing[i]=Buckling_Lenght(int(i),Comprimento_Barra[i],TrussType)
                Lambda_Barra_p[i]=((Dim[i]-2*Esp[i])/Esp[i])/(28.4*epsilon[i]*np.sqrt(ksigma[i]))
                if Class_Enc[i]==4:
                    if Lambda_Barra_p[i]<=0.748:
                        rho_Enc[i]=1
                    else:
                        rho_Enc[i]=(Lambda_Barra_p[i]-0.188)/(Lambda_Barra_p[i]**2)
                        if rho_Enc[i]>=1:
                            rho_Enc[i]=1
                
                if Bracing[i] in ["Symmetrical"]:
                    Comprimento_Enc[i]=L_enc_v
                    if Class_Enc[i]==4:
                        Lambda_v[i]=Comprimento_Enc[i]/iv[i]*np.sqrt(rho_Enc[i]) 
                    else:
                        Lambda_v[i]=Comprimento_Enc[i]/iv[i]
                    Lambda_Ratio_v[i]=Lambda_v[i]/Lambda_1[i]

                    k_v[i]=0.8+Lambda_Ratio_v[i]/10
                    if k_v[i]<=0.9:
                        k_v[i]=0.9
                    elif k_v[i]>=1.0:
                        k_v[i]=1.0
                    Lambda[i]=Lambda_v[i]
                    Lambda_Ratio[i]=Lambda_Ratio_v[i]
                    k_Enc[i]=k_v[i]
                    Lambda_Enc_eff[i]=k_Enc[i]*Lambda_Ratio[i]

                elif Bracing[i] in ["Unsymmetrical"]:    
                    
                    L2_Enc=L_enc_v
                    L1_Enc=L_enc_y
                    if Class_Enc[i]==4:
                        Lambda_v[i]= L2_Enc/iv[i]*np.sqrt(rho_Enc[i]) ### 1993-1-1 altera o lambda em casos de classe 4
                        Lambda_y[i]=L1_Enc/iy[i]*np.sqrt(rho_Enc[i])    
                    else:   
                        Lambda_v[i]=L2_Enc/iv[i]
                        Lambda_y[i]=L1_Enc/iy[i]
                    Lambda_Ratio_v[i]=Lambda_v[i]/Lambda_1[i]
                    Lambda_Ratio_y[i]=Lambda_y[i]/Lambda_1[i]
                    k_v[i]=1.2*(0.8+Lambda_Ratio_v[i]/10)
                    if k_v[i]<=1.08:
                        k_v[i]=1.08
                    elif k_v[i]>=1.20:
                        k_v[i]=1.20
                    k_y[i]=1.2*(0.8+Lambda_Ratio_y[i]/10)
                    if k_y[i]<=1.08:
                        k_y[i]=1.08
                    elif k_y[i]>=1.20:
                        k_y[i]=1.20                    
                    Lambda_Enc_eff_v[i]=k_v[i]*Lambda_Ratio_v[i]
                    Lambda_Enc_eff_y[i]=k_y[i]*Lambda_Ratio_y[i]
                    if Lambda_Enc_eff_v[i]>=Lambda_Enc_eff_y[i]:
                        Lambda[i]=Lambda_v[i]
                        Lambda_Ratio[i]=Lambda_Ratio_v[i]
                        k_Enc[i]=k_v[i]
                        Lambda_Enc_eff[i]=k_Enc[i]*Lambda_Ratio[i]
                        Comprimento_Enc[i]=L2_Enc
                    else:
                        Lambda[i]=Lambda_y[i]
                        Lambda_Ratio[i]=Lambda_Ratio_y[i]
                        k_Enc[i]=k_y[i]
                        Lambda_Enc_eff[i]=k_Enc[i]*Lambda_Ratio[i]
                        Comprimento_Enc[i]=L1_Enc  

                Lambda_Lim[i]=120
                # if Lambda[i]<=Lambda_Lim[i]:
                #     Lambda_Check[i]="OK"
                # else:
                #     Lambda_Check[i]="KO"            
                    
            Aeff[i]=rho_Enc[i]*Area[i]
            Buckling_Curve[i]="b"
            alpha_Enc[i]=Imp_factor_buckling_curve(Buckling_Curve[i])
            Red[i]=1
            if Lambda_Enc_eff[i]*Lambda_1[i]<=Lambda_Lim[i]:
                Lambda_Check[i]="OK"
            else:
                Lambda_Check[i]="KO"
    alpha_Enc=alpha_Enc.astype(float)
    Phi_Enc=0.5*(1+alpha_Enc*(Lambda_Enc_eff-0.2)+Lambda_Enc_eff**2)
    Chi_Enc=1/(Phi_Enc+np.sqrt(Phi_Enc**2-Lambda_Enc_eff**2))
    Chi_Enc[np.where(Chi_Enc>1)]=1

    NbRd=Chi_Enc*Aeff*Sigma*Red*eta_Enc/GammaM1
    Ratio_Enc=np.zeros(len(F_Compressao))
    Ratio_Enc[1:]=np.abs(F_Compressao[1:])/NbRd[1:]
    Check_Enc[np.where(Ratio_Enc<=1)]="OK"
    Check_Enc[np.where(Ratio_Enc>1)]="KO"
    Truss_Out=np.array(["Leg","Diagonal Bar","Horizontal Bar","External Manual Bar","Internal Manual Bar"])
    Enc_Out_csv=Buckling_Output(Class_Enc,nos_vento,TrussType,Truss_Out,Size_Profile,Elements_Matrix,Lambda_Enc_eff,Lambda_1,NTracRd,NbRd,Troco,Ratio_Enc,F_Compressao)
    Enc_Out_Exp=Buckling_Exp_Output(Class_Enc,nos_vento,TrussType,Truss_Out,Size_Profile,Elements_Matrix,Lambda_Enc_eff,Lambda_1,NTracRd,NbRd,Troco,Ratio_Enc,F_Compressao,Comprimento_Enc,Ratio_Class,Lim_Class,Lambda_Lim,k_Enc,alpha_Enc,Phi_Enc,Chi_Enc,Aeff,Red,eta_Enc,Lambda_Check)
    return Ratio_Enc,Enc_Out_csv,Enc_Out_Exp

def Buckling_Function_France(nEle,Inertiav,Inertiau,Inertiay,Inertiaz,Area,Comprimento_Barra,Sigma,Rigidez,nos_vento,Shape,Travamento,Top_bars,Middle_bars,ExpVento,Dim,Esp,TrussType,nos_montante,Conection_Ele,N_Bolt_Ele,GammaM1,F_Compressao,Size_Profile,Elements_Matrix,NTracRd,Troco):
    Buckling_Curve=np.zeros(nEle+1).astype(str)
    Class_Enc=np.zeros(nEle+1)
    Check_Enc=np.zeros(nEle+1).astype(str)
    Comprimento_Enc=np.zeros(nEle+1)
    Lambda_Lim=np.zeros(nEle+1)

    Lambda_v=np.zeros(nEle+1)
    rho_Enc=np.ones(nEle+1)
    Lambda_u=np.zeros(nEle+1)
    Lambda_y=np.zeros(nEle+1)
    Lambda_Barra_p=np.zeros(nEle+1)
    Lambda_z=np.zeros(nEle+1)
    Lambda_Ratio=np.zeros(nEle+1) 
    Lambda_Ratio_v=np.zeros(nEle+1) 
    Lambda_Ratio_u=np.zeros(nEle+1) 
    Lambda_Ratio_y=np.zeros(nEle+1) 
    Lambda_Ratio_z=np.zeros(nEle+1) 
    Lambda_Enc_eff=np.zeros(nEle+1) 
    Lambda_Enc_eff_v=np.zeros(nEle+1) 
    Lambda_Enc_eff_u=np.zeros(nEle+1) 
    Lambda_Enc_eff_y=np.zeros(nEle+1)
    Lambda_Enc_eff_z=np.zeros(nEle+1) 
    Lim_Class=np.zeros(nEle+1) 
    Ratio_Class=np.zeros(nEle+1) 
    Lambda_Check=np.zeros(nEle+1).astype(str)
    alpha_Enc=np.zeros(nEle+1).astype(str)
    k_Enc=np.ones(nEle+1)
    Red=np.ones(nEle+1)
    Bracing=np.ones(nEle+1).astype(str)
    eta_Enc=np.ones(nEle+1)
    Aeff=np.ones(nEle+1)
    ksigma=np.ones(nEle+1)
    k_v=np.ones(nEle+1)
    k_u=np.ones(nEle+1)
    k_y=np.ones(nEle+1)
    k_z=np.ones(nEle+1)
    iv=np.ones(nEle+1)
    iu=np.ones(nEle+1)
    iy=np.ones(nEle+1)
    iz=np.ones(nEle+1)
    iv[1:]=np.sqrt(Inertiav[1:]/Area[1:])
    iu[1:]=np.sqrt(Inertiau[1:]/Area[1:])
    iy[1:]=np.sqrt(Inertiay[1:]/Area[1:])
    iz[1:]=np.sqrt(Inertiaz[1:]/Area[1:])
    Lambda=Comprimento_Barra/np.min(np.vstack([iv, iu, iz, iy]), axis=0)
    epsilon=np.ones(nEle+1)
    Lambda_1=np.ones(nEle+1)
    epsilon[1:]=np.sqrt(235/(Sigma[1:]*10**-6))
    Lambda_1[1:]=np.pi*np.sqrt(Rigidez[1:]/Sigma[1:])
    #####Possivelemente dá para vetorizar isto #####
    nos_vento_inv=nos_vento[::-1,:]
    for i in range(1,len(Class_Enc)):
        if (Shape[i]=="X-Bracing") or (Shape[i]=="Panel Bracing") or ((Shape[i]=="Diamond Bracing") and (Top_bars[i]=="yes") and (Middle_bars[i]=="yes")) or ((Shape[i]=="Diamond Bracing 2") and (Top_bars[i]=="yes") and (Middle_bars[i]=="yes")) or ((Shape[i]=="Alternating Diagonal") and (Top_bars[i]=="yes") and (Middle_bars[i]=="yes")):
            Travamento[i]="Case b"
        elif ((Shape[i]=="Diamond Bracing") and ((Top_bars[i]!="yes") or (Middle_bars[i]!="yes"))) or ((Shape[i]=="Diamond Bracing 2") and ((Top_bars[i]!="yes") or (Middle_bars[i]!="yes"))):
            Travamento[i]="Case a"
        elif ((Shape[i]=="Alternating Diagonal") and ((Top_bars[i]!="yes") or (Middle_bars[i]!="yes"))):
            Travamento[i]="Case d"
        if ExpVento[i]=="Circ":
            Ratio_Class[i]=Dim[i]/Esp[i]
            Lim_Class[i]=90*epsilon[i]**2
            if Dim[i]/Esp[i]<=90*epsilon[i]**2:
                Class_Enc[i]=3
            else:
                Class_Enc[i]=4
            if TrussType[i]!="Leg":
                
                if TrussType[i] in ["Horizontal Bar", "External Manual Bar"]:
                    Multy,Multv,eta_Enc[i]=Horizontal_Buckling(int(i),nos_montante,Conection_Ele[i],N_Bolt_Ele[i],ExpVento[i],TrussType,Troco)
                    Comprimento_Enc[i]=np.max([Comprimento_Barra[i]*Multv,Comprimento_Barra[i]*Multy])
                
                elif TrussType[i]=="Internal Manual Bar" and ops.eleNodes(int(i))[0] in nos_montante and ops.eleNodes(int(i))[1] in nos_montante:
                    Comprimento_Enc[i]=Comprimento_Barra[i]*1
                else:
                    Multy,Multv,eta_Enc[i]=Diagonal_Buckling(int(i),nos_montante,Conection_Ele[i],N_Bolt_Ele[i],ExpVento[i],TrussType,Troco)
                    Comprimento_Enc[i]=np.max([Comprimento_Barra[i]*Multv,Comprimento_Barra[i]*Multy])
                if Conection_Ele[i]=="Yes":
                    k_Enc[i]=0.7
                elif Conection_Ele[i]=="No":
                    k_Enc[i]=0.95
                Lambda_y[i]=Comprimento_Enc[i]/iy[i]

                Lambda[i]=np.max([Lambda_v[i],Lambda_u[i],Lambda_y[i],Lambda_z[i]])
                Lambda_Lim[i]=180
                Lambda_Ratio[i]=Lambda[i]/Lambda_1[i]
                # if Lambda[i]<=Lambda_Lim[i]:
                #     Lambda_Check[i]="OK"
                # else:
                #     Lambda_Check[i]="KO"
            elif TrussType[i]=="Leg":
                L_enc_y,L_enc_v,Bracing[i]=Buckling_Lenght(int(i),Comprimento_Barra[i],TrussType)
                k_Enc[i]=1
                Comprimento_Enc[i]=L_enc_y

                Lambda_y[i]=Comprimento_Enc[i]/iy[i]

                Lambda[i]=np.max([Lambda_v[i],Lambda_u[i],Lambda_y[i],Lambda_z[i]])
                Lambda_Lim[i]=120
                Lambda_Ratio[i]=Lambda[i]/Lambda_1[i]
                # if Lambda[i]<=Lambda_Lim[i]:
                #     Lambda_Check[i]="OK"
                # else:
                #     Lambda_Check[i]="KO"
            Buckling_Curve[i]="c"
            alpha_Enc[i]=Imp_factor_buckling_curve(Buckling_Curve[i])
            Red[i]=1
            Aeff[i]=Area[i]
            Lambda_Enc_eff[i]=k_Enc[i]*Lambda_Ratio[i]
        elif ExpVento[i]=="Flat":
            Ratio_Class[i]=(Dim[i]-2*Esp[i])/Esp[i]
            Lim_Class[i]=15*epsilon[i]
            if ((Dim[i]-2*Esp[i])/Esp[i]<=15*epsilon[i]):
                Class_Enc[i]=3
            else:
                Class_Enc[i]=4   
            ksigma[i]=0.43

            if TrussType[i]!="Leg":
                if TrussType[i] in ["Horizontal Bar", "External Manual Bar"]:
                    Multy,Multv,eta_Enc[i]=Horizontal_Buckling(int(i),nos_montante,Conection_Ele[i],N_Bolt_Ele[i],ExpVento[i],TrussType,Troco)
                    L2_Enc=Comprimento_Barra[i]*Multv
                    L1_Enc=Comprimento_Barra[i]*Multy
                elif TrussType[i]=="Internal Manual Bar" and ops.eleNodes(int(i))[0] in nos_montante and ops.eleNodes(int(i))[1] in nos_montante:
                    L2_Enc=Comprimento_Barra[i]*0.5
                    L1_Enc=Comprimento_Barra[i]*1
                    if Conection_Ele[i]=="No" and N_Bolt_Ele[i]==1:
                        k_Enc[i]=0.9
                    else:
                        k_Enc[i]=1
                else:
                    Multy,Multv,eta_Enc[i]=Diagonal_Buckling(int(i),nos_montante,Conection_Ele[i],N_Bolt_Ele[i],ExpVento[i],TrussType,Troco)
                    L2_Enc=Comprimento_Barra[i]*Multv
                    L1_Enc=Comprimento_Barra[i]*Multy
                Lambda_Barra_p[i]=((Dim[i]-2*Esp[i])/Esp[i])/(28.4*epsilon[i]*np.sqrt(ksigma[i]))
                if Class_Enc[i]==4:
                    if Lambda_Barra_p[i]<=0.748:
                        rho_Enc[i]=1
                    else:
                        rho_Enc[i]=(Lambda_Barra_p[i]-0.188)/(Lambda_Barra_p[i]**2)
                        if rho_Enc[i]>=1:
                            rho_Enc[i]=1
                    Lambda_v[i]=L2_Enc/iv[i]*np.sqrt(rho_Enc[i]) ### 1993-1-1 altera o lambda em casos de classe 4
                    Lambda_u[i]=L2_Enc/iu[i]*np.sqrt(rho_Enc[i])
                    Lambda_y[i]=L1_Enc/iy[i]*np.sqrt(rho_Enc[i])
                    Lambda_z[i]=L1_Enc/iz[i]*np.sqrt(rho_Enc[i])
                    Lambda[i]=np.max([Lambda_v[i],Lambda_u[i],Lambda_y[i],Lambda_z[i]])
                    i_Enc_crit=np.argmax([Lambda_v[i],Lambda_u[i],Lambda_y[i],Lambda_z[i]])      
                else:
                    Lambda_v[i]=L2_Enc/iv[i]
                    Lambda_u[i]=L2_Enc/iu[i]
                    Lambda_y[i]=L1_Enc/iy[i]
                    Lambda_z[i]=L1_Enc/iz[i]
                    Lambda[i]=np.max([Lambda_v[i],Lambda_u[i],Lambda_y[i],Lambda_z[i]])
                    i_Enc_crit=np.argmax([Lambda_v[i],Lambda_u[i],Lambda_y[i],Lambda_z[i]])             
                Lambda_Lim[i]=180
                Lambda_Ratio[i]=Lambda[i]/Lambda_1[i]
                Lambda_Ratio_v[i]=Lambda_v[i]/Lambda_1[i]
                Lambda_Ratio_u[i]=Lambda_u[i]/Lambda_1[i]
                Lambda_Ratio_y[i]=Lambda_y[i]/Lambda_1[i]
                Lambda_Ratio_z[i]=Lambda_z[i]/Lambda_1[i]

                if (Conection_Ele[i]=="No") & (N_Bolt_Ele[i]==1) & (eta_Enc[i]==0.8):
                    k_v[i]=0.7+0.35/Lambda_Ratio_v[i]
                    k_y[i]=0.7+0.58/Lambda_Ratio_y[i]
                    k_z[i]=0.7+0.58/Lambda_Ratio_z[i]
                    k_u[i]=0.0001
                else:
                    k_v[i]=0.7+0.35/Lambda_Ratio_v[i]
                    k_y[i]=0.7+0.4/Lambda_Ratio_y[i]
                    k_z[i]=0.7+0.4/Lambda_Ratio_z[i]
                    k_u[i]=0.0001
                Lambda_Enc_eff_v[i]=k_v[i]*Lambda_Ratio_v[i]
                Lambda_Enc_eff_u[i]=k_u[i]*Lambda_Ratio_u[i]
                Lambda_Enc_eff_y[i]=k_y[i]*Lambda_Ratio_y[i]
                Lambda_Enc_eff_z[i]=k_z[i]*Lambda_Ratio_z[i]
                
                Lambda_Enc_eff[i]=np.max([Lambda_Enc_eff_v[i],Lambda_Enc_eff_u[i],Lambda_Enc_eff_y[i],Lambda_Enc_eff_z[i]])
                i_Enc_crit=np.argmax([Lambda_Enc_eff_v[i],Lambda_Enc_eff_u[i],Lambda_Enc_eff_y[i],Lambda_Enc_eff_z[i]]) 
                k_Enc_Aux=np.array([k_v[i],k_u[i],k_y[i],k_z[i]])
                k_Enc[i]=k_Enc_Aux[i_Enc_crit]
                if i_Enc_crit in [0,1]:
                    Comprimento_Enc[i]=L2_Enc
                else:
                    Comprimento_Enc[i]=L1_Enc
            elif TrussType[i]=="Leg":
                L_enc_y,L_enc_v,Bracing[i]=Buckling_Lenght(int(i),Comprimento_Barra[i],TrussType)
                Lambda_Barra_p[i]=((Dim[i]-2*Esp[i])/Esp[i])/(28.4*epsilon[i]*np.sqrt(ksigma[i]))
                if Class_Enc[i]==4:
                    if Lambda_Barra_p[i]<=0.748:
                        rho_Enc[i]=1
                    else:
                        rho_Enc[i]=(Lambda_Barra_p[i]-0.188)/(Lambda_Barra_p[i]**2)
                        if rho_Enc[i]>=1:
                            rho_Enc[i]=1
                
                if Bracing[i] in ["Symmetrical"]:
                    Comprimento_Enc[i]=L_enc_v
                    if Class_Enc[i]==4:
                        Lambda_v[i]=Comprimento_Enc[i]/iv[i]*np.sqrt(rho_Enc[i]) 
                    else:
                        Lambda_v[i]=Comprimento_Enc[i]/iv[i]
                    Lambda_Ratio_v[i]=Lambda_v[i]/Lambda_1[i]

                    k_v[i]=0.8+Lambda_Ratio_v[i]/10
                    if k_v[i]<=0.9:
                        k_v[i]=0.9
                    elif k_v[i]>=1.0:
                        k_v[i]=1.0
                    Lambda[i]=Lambda_v[i]
                    Lambda_Ratio[i]=Lambda_Ratio_v[i]
                    k_Enc[i]=k_v[i]
                    Lambda_Enc_eff[i]=k_Enc[i]*Lambda_Ratio[i]

                elif Bracing[i] in ["Unsymmetrical"]:    
                    
                    L2_Enc=L_enc_v
                    L1_Enc=L_enc_y
                    if Class_Enc[i]==4:
                        Lambda_v[i]= L2_Enc/iv[i]*np.sqrt(rho_Enc[i]) ### 1993-1-1 altera o lambda em casos de classe 4
                        Lambda_y[i]=L1_Enc/iy[i]*np.sqrt(rho_Enc[i])    
                    else:   
                        Lambda_v[i]=L2_Enc/iv[i]
                        Lambda_y[i]=L1_Enc/iy[i]
                    Lambda_Ratio_v[i]=Lambda_v[i]/Lambda_1[i]
                    Lambda_Ratio_y[i]=Lambda_y[i]/Lambda_1[i]
                    k_v[i]=1.2*(0.8+Lambda_Ratio_v[i]/10)
                    if k_v[i]<=1.08:
                        k_v[i]=1.08
                    elif k_v[i]>=1.20:
                        k_v[i]=1.20
                    k_y[i]=1.2*(0.8+Lambda_Ratio_y[i]/10)
                    if k_y[i]<=1.08:
                        k_y[i]=1.08
                    elif k_y[i]>=1.20:
                        k_y[i]=1.20                    
                    Lambda_Enc_eff_v[i]=k_v[i]*Lambda_Ratio_v[i]
                    Lambda_Enc_eff_y[i]=k_y[i]*Lambda_Ratio_y[i]
                    if Lambda_Enc_eff_v[i]>=Lambda_Enc_eff_y[i]:
                        Lambda[i]=Lambda_v[i]
                        Lambda_Ratio[i]=Lambda_Ratio_v[i]
                        k_Enc[i]=k_v[i]
                        Lambda_Enc_eff[i]=k_Enc[i]*Lambda_Ratio[i]
                        Comprimento_Enc[i]=L2_Enc
                    else:
                        Lambda[i]=Lambda_y[i]
                        Lambda_Ratio[i]=Lambda_Ratio_y[i]
                        k_Enc[i]=k_y[i]
                        Lambda_Enc_eff[i]=k_Enc[i]*Lambda_Ratio[i]
                        Comprimento_Enc[i]=L1_Enc  

                Lambda_Lim[i]=120
                # if Lambda[i]<=Lambda_Lim[i]:
                #     Lambda_Check[i]="OK"
                # else:
                #     Lambda_Check[i]="KO"            
                    
            Aeff[i]=rho_Enc[i]*Area[i]
            Buckling_Curve[i]="b"
            alpha_Enc[i]=Imp_factor_buckling_curve(Buckling_Curve[i])
            Red[i]=1
            if Lambda_Enc_eff[i]*Lambda_1[i]<=Lambda_Lim[i]:
                Lambda_Check[i]="OK"
            else:
                Lambda_Check[i]="KO"
    alpha_Enc=alpha_Enc.astype(float)
    Phi_Enc=0.5*(1+alpha_Enc*(Lambda_Enc_eff-0.2)+Lambda_Enc_eff**2)
    Chi_Enc=1/(Phi_Enc+np.sqrt(Phi_Enc**2-Lambda_Enc_eff**2))
    Chi_Enc[np.where(Chi_Enc>1)]=1

    NbRd=Chi_Enc*Aeff*Sigma*Red*eta_Enc/GammaM1
    Ratio_Enc=np.zeros(len(F_Compressao))
    Ratio_Enc[1:]=np.abs(F_Compressao[1:])/NbRd[1:]
    Check_Enc[np.where(Ratio_Enc<=1)]="OK"
    Check_Enc[np.where(Ratio_Enc>1)]="KO"
    Truss_Out=np.array(["Leg","Diagonal Bar","Horizontal Bar","External Manual Bar","Internal Manual Bar"])
    Enc_Out_csv=Buckling_Output(Class_Enc,nos_vento,TrussType,Truss_Out,Size_Profile,Elements_Matrix,Lambda_Enc_eff,Lambda_1,NTracRd,NbRd,Troco,Ratio_Enc,F_Compressao)
    Enc_Out_Exp=Buckling_Exp_Output(Class_Enc,nos_vento,TrussType,Truss_Out,Size_Profile,Elements_Matrix,Lambda_Enc_eff,Lambda_1,NTracRd,NbRd,Troco,Ratio_Enc,F_Compressao,Comprimento_Enc,Ratio_Class,Lim_Class,Lambda_Lim,k_Enc,alpha_Enc,Phi_Enc,Chi_Enc,Aeff,Red,eta_Enc,Lambda_Check)
    return Ratio_Enc,Enc_Out_csv,Enc_Out_Exp

def Buckling_Function_Spain(nEle,Inertiav,Inertiau,Inertiay,Inertiaz,Area,Comprimento_Barra,Sigma,Rigidez,nos_vento,Shape,Travamento,Top_bars,Middle_bars,ExpVento,Dim,Esp,TrussType,nos_montante,Conection_Ele,N_Bolt_Ele,GammaM1,F_Compressao,Size_Profile,Elements_Matrix,NTracRd,Troco):
    Buckling_Curve=np.zeros(nEle+1).astype(str)
    Class_Enc=np.zeros(nEle+1)
    Check_Enc=np.zeros(nEle+1).astype(str)
    Comprimento_Enc=np.zeros(nEle+1)
    Lambda_Lim=np.zeros(nEle+1)

    Lambda_v=np.zeros(nEle+1)
    rho_Enc=np.ones(nEle+1)
    Lambda_u=np.zeros(nEle+1)
    Lambda_y=np.zeros(nEle+1)
    Lambda_Barra_p=np.zeros(nEle+1)
    Lambda_z=np.zeros(nEle+1)
    Lambda_Ratio=np.zeros(nEle+1) 
    Lambda_Ratio_v=np.zeros(nEle+1) 
    Lambda_Ratio_u=np.zeros(nEle+1) 
    Lambda_Ratio_y=np.zeros(nEle+1) 
    Lambda_Ratio_z=np.zeros(nEle+1) 
    Lambda_Enc_eff=np.zeros(nEle+1) 
    Lambda_Enc_eff_v=np.zeros(nEle+1) 
    Lambda_Enc_eff_u=np.zeros(nEle+1) 
    Lambda_Enc_eff_y=np.zeros(nEle+1)
    Lambda_Enc_eff_z=np.zeros(nEle+1) 
    Lim_Class=np.zeros(nEle+1) 
    Ratio_Class=np.zeros(nEle+1) 
    Lambda_Check=np.zeros(nEle+1).astype(str)
    alpha_Enc=np.zeros(nEle+1).astype(str)
    k_Enc=np.ones(nEle+1)
    Red=np.ones(nEle+1)
    Bracing=np.ones(nEle+1).astype(str)
    eta_Enc=np.ones(nEle+1)
    Aeff=np.ones(nEle+1)
    ksigma=np.ones(nEle+1)
    k_v=np.ones(nEle+1)
    k_u=np.ones(nEle+1)
    k_y=np.ones(nEle+1)
    k_z=np.ones(nEle+1)
    iv=np.ones(nEle+1)
    iu=np.ones(nEle+1)
    iy=np.ones(nEle+1)
    iz=np.ones(nEle+1)
    iv[1:]=np.sqrt(Inertiav[1:]/Area[1:])
    iu[1:]=np.sqrt(Inertiau[1:]/Area[1:])
    iy[1:]=np.sqrt(Inertiay[1:]/Area[1:])
    iz[1:]=np.sqrt(Inertiaz[1:]/Area[1:])
    Lambda=Comprimento_Barra/np.min(np.vstack([iv, iu, iz, iy]), axis=0)
    epsilon=np.ones(nEle+1)
    Lambda_1=np.ones(nEle+1)
    epsilon[1:]=np.sqrt(235/(Sigma[1:]*10**-6))
    Lambda_1[1:]=np.pi*np.sqrt(Rigidez[1:]/Sigma[1:])
    #####Possivelemente dá para vetorizar isto #####
    nos_vento_inv=nos_vento[::-1,:]
    for i in range(1,len(Class_Enc)):
        if (Shape[i]=="X-Bracing") or (Shape[i]=="Panel Bracing") or ((Shape[i]=="Diamond Bracing") and (Top_bars[i]=="yes") and (Middle_bars[i]=="yes")) or ((Shape[i]=="Diamond Bracing 2") and (Top_bars[i]=="yes") and (Middle_bars[i]=="yes")) or ((Shape[i]=="Alternating Diagonal") and (Top_bars[i]=="yes") and (Middle_bars[i]=="yes")):
            Travamento[i]="Case b"
        elif ((Shape[i]=="Diamond Bracing") and ((Top_bars[i]!="yes") or (Middle_bars[i]!="yes"))) or ((Shape[i]=="Diamond Bracing 2") and ((Top_bars[i]!="yes") or (Middle_bars[i]!="yes"))):
            Travamento[i]="Case a"
        elif ((Shape[i]=="Alternating Diagonal") and ((Top_bars[i]!="yes") or (Middle_bars[i]!="yes"))):
            Travamento[i]="Case d"
        if ExpVento[i]=="Circ":
            Ratio_Class[i]=Dim[i]/Esp[i]
            Lim_Class[i]=90*epsilon[i]**2
            if Dim[i]/Esp[i]<=90*epsilon[i]**2:
                Class_Enc[i]=3
            else:
                Class_Enc[i]=4
            if TrussType[i]!="Leg":
                
                if TrussType[i] in ["Horizontal Bar", "External Manual Bar"]:
                    Multy,Multv,eta_Enc[i]=Horizontal_Buckling(int(i),nos_montante,Conection_Ele[i],N_Bolt_Ele[i],ExpVento[i],TrussType,Troco)
                    Comprimento_Enc[i]=np.max([Comprimento_Barra[i]*Multv,Comprimento_Barra[i]*Multy])
                
                elif TrussType[i]=="Internal Manual Bar" and ops.eleNodes(int(i))[0] in nos_montante and ops.eleNodes(int(i))[1] in nos_montante:
                    Comprimento_Enc[i]=Comprimento_Barra[i]*1
                else:
                    Multy,Multv,eta_Enc[i]=Diagonal_Buckling(int(i),nos_montante,Conection_Ele[i],N_Bolt_Ele[i],ExpVento[i],TrussType,Troco)
                    Comprimento_Enc[i]=np.max([Comprimento_Barra[i]*Multv,Comprimento_Barra[i]*Multy])
                if Conection_Ele[i]=="Yes":
                    k_Enc[i]=0.7
                elif Conection_Ele[i]=="No":
                    k_Enc[i]=0.95
                Lambda_y[i]=Comprimento_Enc[i]/iy[i]

                Lambda[i]=np.max([Lambda_v[i],Lambda_u[i],Lambda_y[i],Lambda_z[i]])
                Lambda_Lim[i]=180
                Lambda_Ratio[i]=Lambda[i]/Lambda_1[i]
                # if Lambda[i]<=Lambda_Lim[i]:
                #     Lambda_Check[i]="OK"
                # else:
                #     Lambda_Check[i]="KO"
            elif TrussType[i]=="Leg":
                L_enc_y,L_enc_v,Bracing[i]=Buckling_Lenght(int(i),Comprimento_Barra[i],TrussType)
                k_Enc[i]=1
                Comprimento_Enc[i]=L_enc_y

                Lambda_y[i]=Comprimento_Enc[i]/iy[i]

                Lambda[i]=np.max([Lambda_v[i],Lambda_u[i],Lambda_y[i],Lambda_z[i]])
                Lambda_Lim[i]=120
                Lambda_Ratio[i]=Lambda[i]/Lambda_1[i]
                # if Lambda[i]<=Lambda_Lim[i]:
                #     Lambda_Check[i]="OK"
                # else:
                #     Lambda_Check[i]="KO"
            Buckling_Curve[i]="c"
            alpha_Enc[i]=Imp_factor_buckling_curve(Buckling_Curve[i])
            Red[i]=1
            Aeff[i]=Area[i]
            Lambda_Enc_eff[i]=k_Enc[i]*Lambda_Ratio[i]
        elif ExpVento[i]=="Flat":
            Ratio_Class[i]=(Dim[i]-2*Esp[i])/Esp[i]
            Lim_Class[i]=15*epsilon[i]
            if ((Dim[i]-2*Esp[i])/Esp[i]<=15*epsilon[i]):
                Class_Enc[i]=3
            else:
                Class_Enc[i]=4   
            ksigma[i]=0.43

            if TrussType[i]!="Leg":
                if TrussType[i] in ["Horizontal Bar", "External Manual Bar"]:
                    Multy,Multv,eta_Enc[i]=Horizontal_Buckling(int(i),nos_montante,Conection_Ele[i],N_Bolt_Ele[i],ExpVento[i],TrussType,Troco)
                    L2_Enc=Comprimento_Barra[i]*Multv
                    L1_Enc=Comprimento_Barra[i]*Multy
                elif TrussType[i]=="Internal Manual Bar" and ops.eleNodes(int(i))[0] in nos_montante and ops.eleNodes(int(i))[1] in nos_montante:
                    L2_Enc=Comprimento_Barra[i]*0.5
                    L1_Enc=Comprimento_Barra[i]*1
                    if Conection_Ele[i]=="No" and N_Bolt_Ele[i]==1:
                        k_Enc[i]=0.9
                    else:
                        k_Enc[i]=1
                else:
                    Multy,Multv,eta_Enc[i]=Diagonal_Buckling(int(i),nos_montante,Conection_Ele[i],N_Bolt_Ele[i],ExpVento[i],TrussType,Troco)
                    L2_Enc=Comprimento_Barra[i]*Multv
                    L1_Enc=Comprimento_Barra[i]*Multy
                Lambda_Barra_p[i]=((Dim[i]-2*Esp[i])/Esp[i])/(28.4*epsilon[i]*np.sqrt(ksigma[i]))
                if Class_Enc[i]==4:
                    if Lambda_Barra_p[i]<=0.748:
                        rho_Enc[i]=1
                    else:
                        rho_Enc[i]=(Lambda_Barra_p[i]-0.188)/(Lambda_Barra_p[i]**2)
                        if rho_Enc[i]>=1:
                            rho_Enc[i]=1
                    Lambda_v[i]=L2_Enc/iv[i]*np.sqrt(rho_Enc[i]) ### 1993-1-1 altera o lambda em casos de classe 4
                    Lambda_u[i]=L2_Enc/iu[i]*np.sqrt(rho_Enc[i])
                    Lambda_y[i]=L1_Enc/iy[i]*np.sqrt(rho_Enc[i])
                    Lambda_z[i]=L1_Enc/iz[i]*np.sqrt(rho_Enc[i])
                    Lambda[i]=np.max([Lambda_v[i],Lambda_u[i],Lambda_y[i],Lambda_z[i]])
                    i_Enc_crit=np.argmax([Lambda_v[i],Lambda_u[i],Lambda_y[i],Lambda_z[i]])      
                else:
                    Lambda_v[i]=L2_Enc/iv[i]
                    Lambda_u[i]=L2_Enc/iu[i]
                    Lambda_y[i]=L1_Enc/iy[i]
                    Lambda_z[i]=L1_Enc/iz[i]
                    Lambda[i]=np.max([Lambda_v[i],Lambda_u[i],Lambda_y[i],Lambda_z[i]])
                    i_Enc_crit=np.argmax([Lambda_v[i],Lambda_u[i],Lambda_y[i],Lambda_z[i]])             
                Lambda_Lim[i]=180
                Lambda_Ratio[i]=Lambda[i]/Lambda_1[i]
                Lambda_Ratio_v[i]=Lambda_v[i]/Lambda_1[i]
                Lambda_Ratio_u[i]=Lambda_u[i]/Lambda_1[i]
                Lambda_Ratio_y[i]=Lambda_y[i]/Lambda_1[i]
                Lambda_Ratio_z[i]=Lambda_z[i]/Lambda_1[i]

                if (Conection_Ele[i]=="No") & (N_Bolt_Ele[i]==1) & (eta_Enc[i]==0.8):
                    k_v[i]=0.7+0.35/Lambda_Ratio_v[i]
                    k_y[i]=0.7+0.58/Lambda_Ratio_y[i]
                    k_z[i]=0.7+0.58/Lambda_Ratio_z[i]
                    k_u[i]=0.0001
                else:
                    k_v[i]=0.7+0.35/Lambda_Ratio_v[i]
                    k_y[i]=0.7+0.4/Lambda_Ratio_y[i]
                    k_z[i]=0.7+0.4/Lambda_Ratio_z[i]
                    k_u[i]=0.0001
                Lambda_Enc_eff_v[i]=k_v[i]*Lambda_Ratio_v[i]
                Lambda_Enc_eff_u[i]=k_u[i]*Lambda_Ratio_u[i]
                Lambda_Enc_eff_y[i]=k_y[i]*Lambda_Ratio_y[i]
                Lambda_Enc_eff_z[i]=k_z[i]*Lambda_Ratio_z[i]
                
                Lambda_Enc_eff[i]=np.max([Lambda_Enc_eff_v[i],Lambda_Enc_eff_u[i],Lambda_Enc_eff_y[i],Lambda_Enc_eff_z[i]])
                i_Enc_crit=np.argmax([Lambda_Enc_eff_v[i],Lambda_Enc_eff_u[i],Lambda_Enc_eff_y[i],Lambda_Enc_eff_z[i]]) 
                k_Enc_Aux=np.array([k_v[i],k_u[i],k_y[i],k_z[i]])
                k_Enc[i]=k_Enc_Aux[i_Enc_crit]
                if i_Enc_crit in [0,1]:
                    Comprimento_Enc[i]=L2_Enc
                else:
                    Comprimento_Enc[i]=L1_Enc
            elif TrussType[i]=="Leg":
                L_enc_y,L_enc_v,Bracing[i]=Buckling_Lenght(int(i),Comprimento_Barra[i],TrussType)
                Lambda_Barra_p[i]=((Dim[i]-2*Esp[i])/Esp[i])/(28.4*epsilon[i]*np.sqrt(ksigma[i]))
                if Class_Enc[i]==4:
                    if Lambda_Barra_p[i]<=0.748:
                        rho_Enc[i]=1
                    else:
                        rho_Enc[i]=(Lambda_Barra_p[i]-0.188)/(Lambda_Barra_p[i]**2)
                        if rho_Enc[i]>=1:
                            rho_Enc[i]=1
                
                if Bracing[i] in ["Symmetrical"]:
                    Comprimento_Enc[i]=L_enc_v
                    if Class_Enc[i]==4:
                        Lambda_v[i]=Comprimento_Enc[i]/iv[i]*np.sqrt(rho_Enc[i]) 
                    else:
                        Lambda_v[i]=Comprimento_Enc[i]/iv[i]
                    Lambda_Ratio_v[i]=Lambda_v[i]/Lambda_1[i]

                    k_v[i]=0.8+Lambda_Ratio_v[i]/10
                    if k_v[i]<=0.9:
                        k_v[i]=0.9
                    elif k_v[i]>=1.0:
                        k_v[i]=1.0
                    Lambda[i]=Lambda_v[i]
                    Lambda_Ratio[i]=Lambda_Ratio_v[i]
                    k_Enc[i]=k_v[i]
                    Lambda_Enc_eff[i]=k_Enc[i]*Lambda_Ratio[i]

                elif Bracing[i] in ["Unsymmetrical"]:    
                    
                    L2_Enc=L_enc_v
                    L1_Enc=L_enc_y
                    if Class_Enc[i]==4:
                        Lambda_v[i]= L2_Enc/iv[i]*np.sqrt(rho_Enc[i]) ### 1993-1-1 altera o lambda em casos de classe 4
                        Lambda_y[i]=L1_Enc/iy[i]*np.sqrt(rho_Enc[i])    
                    else:   
                        Lambda_v[i]=L2_Enc/iv[i]
                        Lambda_y[i]=L1_Enc/iy[i]
                    Lambda_Ratio_v[i]=Lambda_v[i]/Lambda_1[i]
                    Lambda_Ratio_y[i]=Lambda_y[i]/Lambda_1[i]
                    k_v[i]=1.2*(0.8+Lambda_Ratio_v[i]/10)
                    if k_v[i]<=1.08:
                        k_v[i]=1.08
                    elif k_v[i]>=1.20:
                        k_v[i]=1.20
                    k_y[i]=1.2*(0.8+Lambda_Ratio_y[i]/10)
                    if k_y[i]<=1.08:
                        k_y[i]=1.08
                    elif k_y[i]>=1.20:
                        k_y[i]=1.20                    
                    Lambda_Enc_eff_v[i]=k_v[i]*Lambda_Ratio_v[i]
                    Lambda_Enc_eff_y[i]=k_y[i]*Lambda_Ratio_y[i]
                    if Lambda_Enc_eff_v[i]>=Lambda_Enc_eff_y[i]:
                        Lambda[i]=Lambda_v[i]
                        Lambda_Ratio[i]=Lambda_Ratio_v[i]
                        k_Enc[i]=k_v[i]
                        Lambda_Enc_eff[i]=k_Enc[i]*Lambda_Ratio[i]
                        Comprimento_Enc[i]=L2_Enc
                    else:
                        Lambda[i]=Lambda_y[i]
                        Lambda_Ratio[i]=Lambda_Ratio_y[i]
                        k_Enc[i]=k_y[i]
                        Lambda_Enc_eff[i]=k_Enc[i]*Lambda_Ratio[i]
                        Comprimento_Enc[i]=L1_Enc  

                Lambda_Lim[i]=120
                # if Lambda[i]<=Lambda_Lim[i]:
                #     Lambda_Check[i]="OK"
                # else:
                #     Lambda_Check[i]="KO"            
                    
            Aeff[i]=rho_Enc[i]*Area[i]
            Buckling_Curve[i]="b"
            alpha_Enc[i]=Imp_factor_buckling_curve(Buckling_Curve[i])
            Red[i]=1
            if Lambda_Enc_eff[i]*Lambda_1[i]<=Lambda_Lim[i]:
                Lambda_Check[i]="OK"
            else:
                Lambda_Check[i]="KO"
    alpha_Enc=alpha_Enc.astype(float)
    Phi_Enc=0.5*(1+alpha_Enc*(Lambda_Enc_eff-0.2)+Lambda_Enc_eff**2)
    Chi_Enc=1/(Phi_Enc+np.sqrt(Phi_Enc**2-Lambda_Enc_eff**2))
    Chi_Enc[np.where(Chi_Enc>1)]=1

    NbRd=Chi_Enc*Aeff*Sigma*Red*eta_Enc/GammaM1
    Ratio_Enc=np.zeros(len(F_Compressao))
    Ratio_Enc[1:]=np.abs(F_Compressao[1:])/NbRd[1:]
    Check_Enc[np.where(Ratio_Enc<=1)]="OK"
    Check_Enc[np.where(Ratio_Enc>1)]="KO"
    Truss_Out=np.array(["Leg","Diagonal Bar","Horizontal Bar","External Manual Bar","Internal Manual Bar"])
    Enc_Out_csv=Buckling_Output(Class_Enc,nos_vento,TrussType,Truss_Out,Size_Profile,Elements_Matrix,Lambda_Enc_eff,Lambda_1,NTracRd,NbRd,Troco,Ratio_Enc,F_Compressao)
    Enc_Out_Exp=Buckling_Exp_Output(Class_Enc,nos_vento,TrussType,Truss_Out,Size_Profile,Elements_Matrix,Lambda_Enc_eff,Lambda_1,NTracRd,NbRd,Troco,Ratio_Enc,F_Compressao,Comprimento_Enc,Ratio_Class,Lim_Class,Lambda_Lim,k_Enc,alpha_Enc,Phi_Enc,Chi_Enc,Aeff,Red,eta_Enc,Lambda_Check)
    return Ratio_Enc,Enc_Out_csv,Enc_Out_Exp

