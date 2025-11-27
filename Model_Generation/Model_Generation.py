import numpy as np
from openseespy import opensees as ops
from Model_Generation.Node_Filters import Node_Filter_NCN,Node_Filter_Montante_TrussType,Node_Constraint_func
from Sections_Materials.Section_Properties import Section,get_fu_fy
from Utilities.Utilities import angulo_entre_vetores
from Utilities.Utilities_ops import comprimentos_elementos_array
def Model_Generation(Nodes_Matrix, Elements_Matrix):  
    ##Geração do modelo  
    Node_Type=np.zeros(len(Nodes_Matrix)+1).astype(str)

    nos_grupo, z_grupo,Node_Type,coords=Node_Filter_NCN(Nodes_Matrix)


    nEle=Elements_Matrix.shape[0]
    zmed=np.zeros(nEle+1)
    zmax=np.zeros(nEle+1)
    Dim=np.zeros(nEle+1)
    Esp=np.zeros(nEle+1)
    Area=np.zeros(nEle+1)
    Inertiay=np.zeros(nEle+1)
    Sigma=np.zeros(nEle+1)
    Shape=np.zeros(nEle+1).astype(str)
    Top_bars=np.zeros(nEle+1).astype(str)
    Middle_bars=np.zeros(nEle+1).astype(str)
    Sigma_u=np.zeros(nEle+1)
    Conection_Ele=np.zeros(nEle+1).astype(str)
    N_Bolt_Ele=np.zeros(nEle+1)
    Rigidez=np.zeros(nEle+1)
    Inertiaz=np.zeros(nEle+1)
    munit=np.zeros(nEle+1)
    
    ExpVento=np.zeros(nEle+1).astype(str)
    TrussType=np.zeros(nEle+1).astype(str)
    Type_Profile=np.zeros(nEle+1).astype(str)
    Size_Profile=np.zeros(nEle+1).astype(str)

    Inertiav=np.zeros(nEle+1)
    Inertiau=np.zeros(nEle+1)
    ops.geomTransf("Linear", 1, 0, 1, 0)
    ops.geomTransf("Linear", 2, 1, 0, 0)

    nos_montante,TrussType= Node_Filter_Montante_TrussType(Elements_Matrix)

    for i in range(len(Elements_Matrix)):
        eleTag = int(Elements_Matrix[i,0])
        ni = int(Elements_Matrix[i,2])
        nj = int(Elements_Matrix[i,3])
        Type_Profile[eleTag]=Elements_Matrix[i,5]
        Size_Profile[eleTag]=Elements_Matrix[i,6]
        Dim[eleTag],Esp[eleTag],Area[eleTag],Inertiay[eleTag],Inertiaz[eleTag],Inertiav[eleTag],Inertiau[eleTag],munit[eleTag],ExpVento[eleTag]=Section(Type_Profile[eleTag],Size_Profile[eleTag],7850)
        E = 210*10**9
        Sigma_u[eleTag],Sigma[eleTag]=get_fu_fy(Elements_Matrix[i,7])
        Conection_Ele[eleTag]=Elements_Matrix[i,8]
        N_Bolt_Ele[eleTag]=Elements_Matrix[i,9]
        Shape[eleTag]=Elements_Matrix[i,10]
        Top_bars[eleTag]=Elements_Matrix[i,11]
        Middle_bars[eleTag]=Elements_Matrix[i,12]
        matTag = 1000 + eleTag
        ops.uniaxialMaterial("Elastic", matTag, E)
        #ops.element("truss", eleTag, ni, nj,Area[eleTag] , matTag, "-rho", munit[eleTag])
        #G=E/(2*(1+0.3))
        G=81*10**9
        J=10**-9
        Rigidez[eleTag]=E
        coords_i = np.array(ops.nodeCoord(ni))
        coords_j = np.array(ops.nodeCoord(nj))
        Vetor_Bar=coords_j-coords_i
        Prod=np.abs(np.dot(Vetor_Bar, [0,1,0]))

        Angle_Bar=angulo_entre_vetores(Vetor_Bar,[0,1,0])
        if Angle_Bar==0 or Angle_Bar==180:
            Trans=2
        else:
            Trans=1
        if TrussType[eleTag]=="Leg":
            ops.element("elasticBeamColumn", eleTag, ni, nj, Area[eleTag], E, G, J, Inertiay[eleTag], Inertiaz[eleTag], Trans,"-mass", munit[eleTag])
        elif (TrussType[eleTag]!="Leg" and int(ni) in nos_montante and int(nj) in nos_montante) or (TrussType[eleTag]=="Internal Manual Bar") or (TrussType[eleTag]=="External Manual Bar"):
            ops.element("elasticBeamColumn", eleTag, ni, nj, Area[eleTag], E, G, J, Inertiay[eleTag], Inertiaz[eleTag],Trans,"-mass", munit[eleTag],"-release", 1, 4, 5,6,"-release", 2, 4, 5,6)
            #print("Release ",eleTag , ": ", ni," e ", nj)
        elif TrussType[eleTag]!="Leg" and int(ni) in nos_montante:
            ops.element("elasticBeamColumn", eleTag, ni, nj, Area[eleTag], E, G, J, Inertiay[eleTag], Inertiaz[eleTag], Trans,"-mass", munit[eleTag],"-release", 1, 4, 5,6)
            #print("Release ",eleTag , ": ", ni)
        elif TrussType[eleTag]!="Leg" and int(nj) in nos_montante:
            ops.element("elasticBeamColumn", eleTag, ni, nj, Area[eleTag], E, G, J, Inertiay[eleTag], Inertiaz[eleTag], Trans,"-mass", munit[eleTag],"-release", 2, 4, 5,6)
            #print("Release ",eleTag , ": ", nj)
        
        zmax[eleTag]=max(ops.nodeCoord(ni)[2],ops.nodeCoord(nj)[2])
        zmed[eleTag]=np.round((ops.nodeCoord(ni)[2]+ops.nodeCoord(nj)[2])/2,5)

    Comprimento_Barra=comprimentos_elementos_array()
    ElementosID = np.array(ops.getEleTags(), dtype=int)  
    Sigma_u=Sigma_u*10**6
    Sigma=Sigma*10**6

    Node_Constraint,torre= Node_Constraint_func(Nodes_Matrix,Node_Type)

    return Type_Profile, Size_Profile, Dim, Esp, Area, Inertiay,Inertiaz,Inertiav,Inertiau,munit,ExpVento,E,Sigma_u,Sigma, Conection_Ele,N_Bolt_Ele,Shape,Top_bars,Middle_bars,zmax,zmed,Comprimento_Barra,ElementosID ,nos_grupo, z_grupo,nEle,Node_Type,torre,TrussType, nos_montante,Rigidez,coords