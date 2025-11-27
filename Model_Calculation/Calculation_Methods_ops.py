from openseespy import opensees as ops
import numpy as np

def Create_Load_Case(i,Nome):
    ops.timeSeries("Constant", i)
    ops.pattern("Plain", i, i)
    return Nome

def Calc_Load_Case(i,coords,nos_por_elemento,todos_nos, nEle):
    # === 5. Análise estática ===
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-6, 10)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1)
    ops.analysis("Static")
    ops.analyze(1)
    ops.reactions()
    F_Func=np.zeros((nEle+1,12))
    F_axial_Func=np.zeros((nEle+1))
    Desloc_Func=np.zeros((len(todos_nos)+2,6))
    Reactions_Func=np.zeros((len(todos_nos)+2,6))
    Desloc_aux_Func=np.zeros((len(todos_nos)+2))
    for tag in coords.keys():
        u = ops.nodeDisp(tag)
        Desloc_Func[tag,:] = u
        Reactions_Func[tag,:]= ops.nodeReaction(tag)
        Desloc_aux_Func[tag]=np.sqrt(u[0]**2+u[1]**2)
    for ele, nos in nos_por_elemento.items():
        F_Func[ele,:] = ops.eleForce(ele)
        F_axial_Func[ele]=ops.eleResponse(int(ele), 'basicForce')[0]
    
    ops.remove('loadPattern',i)
    ops.wipeAnalysis()
    
    return Desloc_Func,Reactions_Func,Desloc_aux_Func,F_Func,F_axial_Func