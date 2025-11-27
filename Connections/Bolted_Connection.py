import numpy as np
from Sections_Materials.Section_Properties import get_bolt_props,get_fu_fy,nut_across_flats
from Output.Output_Lattice_Reduced import Connection_Bracing_Output,Connection_Shaft_Output,Flanged_Shaft_Output
from Output.Output_Lattice_Expanded import Connection_Bracing_Exp_Output,Connection_Shaft_Exp_Output,Flanged_Shaft_Exp_Output
from Utilities.Utilities_ops import elementos_do_no_filt
def Bolted_Connection_Bracing(Dig_Hor_Matrix,F_Tracao,F_Compressao,TrussType,Dim,GammaM2,ExpVento,Esp,Area,Leg_Type,Sigma_u,Sigma,GammaM0,Troco,nos_vento):    
    ID_Lig_Bolted=len(Dig_Hor_Matrix)
    alpha_v_Bolt=np.zeros(ID_Lig_Bolted+1)

    As_bolt=np.zeros(ID_Lig_Bolted+1)
    A_bolt=np.zeros(ID_Lig_Bolted+1)
    d_0=np.zeros(ID_Lig_Bolted+1)
    d_Bolt=np.zeros(ID_Lig_Bolted+1)
    FvRd=np.zeros(ID_Lig_Bolted+1)
    Ratio_bolt_lig=np.zeros(ID_Lig_Bolted+1)
    Ratio_Esmag_Lig=np.zeros(ID_Lig_Bolted+1)
    Ratio_Block_Lig=np.zeros(ID_Lig_Bolted+1)
    Fele= np.maximum(np.abs(F_Tracao), np.abs(F_Compressao))
    alpha_d_end=np.zeros(ID_Lig_Bolted+1)
    alpha_d_inner=np.zeros(ID_Lig_Bolted+1)
    k1_edge=np.zeros(ID_Lig_Bolted+1)
    k1_inner=np.zeros(ID_Lig_Bolted+1)
    alpha_b_end=np.zeros(ID_Lig_Bolted+1)
    alpha_b_inner=np.zeros(ID_Lig_Bolted+1)
    Corte_Crit=np.zeros(ID_Lig_Bolted+1).astype(str)
    FbRd=np.zeros(ID_Lig_Bolted+1)
    FbRd_Edge_Inner=np.zeros(ID_Lig_Bolted+1)
    FbRd_Edge_End=np.zeros(ID_Lig_Bolted+1)
    Bolt_Steel=np.zeros(ID_Lig_Bolted+1)
    Veff1Rd=np.zeros(ID_Lig_Bolted+1)
    Lig_Type=np.zeros(ID_Lig_Bolted+1).astype(str)
    Veff2Rd=np.zeros(ID_Lig_Bolted+1)
    Anv_block=np.zeros(ID_Lig_Bolted+1)
    beta_Cant=np.zeros(ID_Lig_Bolted+1)
    t_gusset=np.zeros(ID_Lig_Bolted+1)
    Ant_block=np.zeros(ID_Lig_Bolted+1)
    fub=np.zeros(ID_Lig_Bolted+1)
    Bolt_Steel=np.zeros(ID_Lig_Bolted+1).astype(str)
    fyb=np.zeros(ID_Lig_Bolted+1)
    NuRd=np.zeros(ID_Lig_Bolted+1)
    Anet=np.zeros(ID_Lig_Bolted+1)
    fu_Lig,fy_Lig=np.zeros(ID_Lig_Bolted+1),np.zeros(ID_Lig_Bolted+1)
    N_Bolt_Lig=np.zeros(ID_Lig_Bolted+1)
    Ele_Lig=np.zeros(ID_Lig_Bolted+1)
    Steel_Lig=np.zeros(ID_Lig_Bolted+1).astype(str)
    Corte_Lig=np.ones(ID_Lig_Bolted+1)
    Bolt=np.ones(ID_Lig_Bolted+1).astype(str)
    p1_gusset=np.zeros(ID_Lig_Bolted+1)
    e2_gusset=np.zeros(ID_Lig_Bolted+1)
    Bolt_Lig=np.zeros(ID_Lig_Bolted+1).astype(str)
    Calc_Bolt=np.zeros(ID_Lig_Bolted+1).astype(str)
    Calc_Esmag=np.zeros(ID_Lig_Bolted+1).astype(str)
    Calc_Block=np.zeros(ID_Lig_Bolted+1).astype(str)
    e1_gusset=np.zeros(ID_Lig_Bolted+1)
    p2_gusset=np.zeros(ID_Lig_Bolted+1)
    Check_Bolt=np.zeros(ID_Lig_Bolted+1).astype(str)
    Check_Esmag=np.zeros(ID_Lig_Bolted+1).astype(str)
    Check_Block=np.zeros(ID_Lig_Bolted+1).astype(str)
    
    Check_Lig=np.zeros(ID_Lig_Bolted+1).astype(str)

    Calc_Nu=np.zeros(ID_Lig_Bolted+1).astype(str)
    Check_DisMin=np.zeros(ID_Lig_Bolted+1).astype(str)
    TrussOut=np.zeros(ID_Lig_Bolted+1).astype(str)
    Ratio_Lig=np.zeros(ID_Lig_Bolted+1)
    DimOut=np.zeros(ID_Lig_Bolted+1)
    EspOut=np.zeros(ID_Lig_Bolted+1)
    FLig_Rd=np.zeros(ID_Lig_Bolted+1)
    N_Bolt_Lig[1:]=Dig_Hor_Matrix[:,3]
    Ele_Lig[1:]=Dig_Hor_Matrix[:,1]
    Lig_Type[1:]=Dig_Hor_Matrix[:,2]
    Bolt_Steel[1:]=Dig_Hor_Matrix[:,5]
    p1_gusset[1:]=Dig_Hor_Matrix[:,9].astype(float)*10**-3
    e1_gusset[1:]=Dig_Hor_Matrix[:,7].astype(float)*10**-3
    t_gusset[1:]=Dig_Hor_Matrix[:,11].astype(float)*10**-3
    e2_gusset[1:]=Dig_Hor_Matrix[:,8].astype(float)*10**-3
    p2_gusset[1:]=Dig_Hor_Matrix[:,10].astype(float)*10**-3
    Bolt[1:]=Dig_Hor_Matrix[:,4].astype(str)
    Steel_Lig[1:]=Dig_Hor_Matrix[:,6]




    for i in range(1,ID_Lig_Bolted+1):
        Calc_Bolt[i]="Sim"
        Calc_Esmag[i]="Sim"
        Calc_Block[i]="Sim"
        Calc_Nu[i]="Sim"
        TrussOut[i]=TrussType[int(Ele_Lig[i])]
        DimOut[i]=Dim[int(Ele_Lig[i])]
        EspOut[i]=Esp[int(Ele_Lig[i])]
        if Lig_Type[i]=="No":
            if e1_gusset[i]==0 or e2_gusset[i]==0:
                Calc_Bolt[i]="Nao"
                Calc_Esmag[i]="Nao"
                Calc_Block[i]="Nao"
                Calc_Nu[i]="Nao"
            
            if Bolt[i]=="0" or Bolt_Steel[i]=="0" or N_Bolt_Lig[i]=="0":
                Calc_Bolt[i]="Nao"
                Calc_Esmag[i]="Nao"
                Calc_Block[i]="Nao"
                Calc_Nu[i]="Nao"
            if Bolt[i]!="0" and Bolt_Steel[i]!="0" and N_Bolt_Lig[i]!="0":
                Calc_Bolt[i]="Sim"
            if Bolt[i]!="0" and N_Bolt_Lig[i]!="0":
                Calc_Bolt[i]="Sim"
                Bolt_Steel[i]!="8.8"

            if Calc_Bolt[i]=="Sim":
                ###### Verificação parafuso ######
                fub[i],fyb[i]=get_fu_fy(Bolt_Steel[i])
                fub[i]=fub[i]*10**6
                fyb[i]=fyb[i]*10**6
                alpha_v_Bolt[i]=0.6
                d_Bolt[i], d_0[i], As_bolt[i], A_bolt[i]=get_bolt_props(Bolt[i])
                FvRd[i]=Corte_Lig[i]*alpha_v_Bolt[i]*fub[i]*As_bolt[i]*N_Bolt_Lig[i]/GammaM2
                Ratio_bolt_lig[i]=Fele[int(Ele_Lig[i])]/FvRd[i]

            if Calc_Esmag[i]=="Sim":       
                ##### Validação distancias minimas
                if e1_gusset[i]!=0 and e1_gusset[i]<=1.2*d_0[i]:
                    Check_DisMin[i]="KO"
                elif e2_gusset[i]!=0 and e2_gusset[i]<=1.2*d_0[i]:
                    Check_DisMin[i]="KO"
                elif p1_gusset[i]!=0 and p1_gusset[i]<=2.2*d_0[i]:
                    Check_DisMin[i]="KO"
                elif p2_gusset[i]!=0 and p2_gusset[i]<=2.4*d_0[i]:
                    Check_DisMin[i]="KO"
                else:
                    Check_DisMin[i]="OK"
                ##### Verificação corte #######
                if ExpVento[int(Ele_Lig[i])]=="Flat":
                    if t_gusset[i]==0:
                        t_gusset[i]=Esp[int(Ele_Lig[i])]
                    else:
                        t_gusset[i]=np.min([Esp[int(Ele_Lig[i])],t_gusset[i]])

                if Leg_Type=="Flat" and ExpVento[int(Ele_Lig[i])]=="Flat":
                    ####Cantoneira-Cantoneira
                    fu_Lig[i],fy_Lig[i]=Sigma_u[int(Ele_Lig[i])],Sigma[int(Ele_Lig[i])]
                else:
                    if Steel_Lig[i]!="0":
                        fu_Lig[i],fy_Lig[i]=get_fu_fy(Steel_Lig[i])
                        fu_Lig[i]=fu_Lig[i]*10**6
                        fy_Lig[i]=fy_Lig[i]*10**6
                    else:
                        fu_Lig[i],fy_Lig[i]=Sigma_u[int(Ele_Lig[i])],Sigma[int(Ele_Lig[i])]
                if ((p1_gusset[i]!=0) and (p2_gusset[i]!=0)):

                    alpha_d_end[i]=e1_gusset[i]/(3*d_0[i])
                    alpha_d_inner[i]=p1_gusset[i]/(3*d_0[i])-1/4
                    k1_edge[i]=np.min([2.8*e2_gusset[i]/d_0[i]-1.7,1.4*p2_gusset[i]/d_0[i]-1.7,2.5])
                    k1_inner[i]=np.min([1.4*p2_gusset[i]/d_0[i]-1.7,2.5])     
                    alpha_b_end[i]=np.min([alpha_d_end[i],fub[i]/fu_Lig[i],1.0])
                    alpha_b_inner[i]=np.min([alpha_d_inner[i],fub[i]/fu_Lig[i],1.0])
                    FbRd_Edge_Inner[i]=alpha_b_inner[i]* k1_edge[i]*fu_Lig[i]*d_Bolt[i]*t_gusset[i]/GammaM2
                    FbRd_Edge_End[i]=alpha_b_end[i]* k1_edge[i]*fu_Lig[i]*d_Bolt[i]*t_gusset[i]/GammaM2            
                    FbRd[i]=FbRd_Edge_Inner[i]*(N_Bolt_Lig[i]-2)+FbRd_Edge_End[i]*(2)
                elif ((p2_gusset[i]==0) and (p1_gusset[i]!=0)):
                    alpha_d_end[i]=e1_gusset[i]/(3*d_0[i])
                    alpha_d_inner[i]=p1_gusset[i]/(3*d_0[i])-1/4
                    k1_edge[i]=np.min([2.8*e2_gusset[i]/d_0[i]-1.7,2.5])
                    k1_inner[i]=np.min([2.5])     
                    alpha_b_end[i]=np.min([alpha_d_end[i],fub[i]/fu_Lig[i],1.0])
                    alpha_b_inner[i]=np.min([alpha_d_inner[i],fub[i]/fu_Lig[i],1.0])
                    FbRd_Edge_Inner[i]=alpha_b_inner[i]* k1_edge[i]*fu_Lig[i]*d_Bolt[i]*t_gusset[i]/GammaM2
                    FbRd_Edge_End[i]=alpha_b_end[i]* k1_edge[i]*fu_Lig[i]*d_Bolt[i]*t_gusset[i]/GammaM2            
                    FbRd[i]=FbRd_Edge_Inner[i]*(N_Bolt_Lig[i]-1)+FbRd_Edge_End[i]*(1)
                elif ((p1_gusset[i]==0) and (p2_gusset[i]==0)):
                    alpha_d_end[i]=e1_gusset[i]/(3*d_0[i])
                    k1_edge[i]=np.min([2.8*e2_gusset[i]/d_0[i]-1.7,2.5]) 
                    alpha_b_end[i]=np.min([alpha_d_end[i],fub[i]/fu_Lig[i],1.0])
                    FbRd_Edge_Inner[i]=0
                    FbRd_Edge_End[i]=alpha_b_end[i]* k1_edge[i]*fu_Lig[i]*d_Bolt[i]*t_gusset[i]/GammaM2  
                    FbRd[i]=FbRd_Edge_End[i]
                elif ((p2_gusset[i]!=0) and (p1_gusset[i]==0)):
                    alpha_d_end[i]=e1_gusset[i]/(3*d_0[i])
                    k1_edge[i]=np.min([2.8*e2_gusset[i]/d_0[i]-1.7,1.4*p2_gusset[i]/d_0[i]-1.7,2.5])    
                    alpha_b_end[i]=np.min([alpha_d_end[i],fub[i]/fu_Lig[i],1.0])
                    FbRd_Edge_Inner[i]=0
                    FbRd_Edge_End[i]=alpha_b_end[i]* k1_edge[i]*fu_Lig[i]*d_Bolt[i]*t_gusset[i]/GammaM2 
                    FbRd[i]=FbRd_Edge_End[i]*N_Bolt_Lig[i]
                Ratio_Esmag_Lig[i]=Fele[int(Ele_Lig[i])]/FbRd[i]
            if Calc_Block[i]=="Sim":                        
                ##### Verificação Bloco #######
                if p2_gusset[i]!=0:
                    # if p1_gusset[i]==0:
                    #     Anv_block[i]=t_gusset[i]*(e1_gusset[i]-0.5*d_0[i])
                    #     Ant_block[i]=t_gusset[i]*(e2_gusset[i]+(N_Bolt_Lig[i]-1)*p2_gusset[i]-((N_Bolt_Lig[i]-0.5)*d_0[i]))
                    # elif p2_gusset[i]==0:
                    #     Anv_block[i]=t_gusset[i]*(e1_gusset[i]+(N_Bolt_Lig[i]-1)*p1_gusset[i]-((N_Bolt_Lig[i]-0.5)*d_0[i]))
                    #     Ant_block[i]=t_gusset[i]*(e2_gusset[i]-0.5*d_0[i])
                    # elif (p1_gusset[i]!=0) and (p1_gusset[i]!=0):
                    #     Anv_block[i]=t_gusset[i]*(e1_gusset[i]+(N_Bolt_Lig[i]/2-1)*p1_gusset[i]-((N_Bolt_Lig[i]/2-0.5)*d_0[i]))            
                    #     Ant_block[i]=t_gusset[i]*(e2_gusset[i]+(N_Bolt_Lig[i]/2-1)*p2_gusset[i]-((N_Bolt_Lig[i]/2-0.5)*d_0[i]))
                    Anv_block[i]= t_gusset[i]*(e1_gusset[i]+ (N_Bolt_Lig[i]-2)/2*p1_gusset[i] - (N_Bolt_Lig[i]/2-0.5)*d_0[i])*2
                    Ant_block[i]= t_gusset[i]*(p2_gusset[i]-d_0[i])
                else:
                    Anv_block[i]=0            
                    Ant_block[i]=0          
                Veff1Rd[i]=fu_Lig[i]*Ant_block[i]/GammaM2+(1/np.sqrt(3))*fy_Lig[i]*Anv_block[i]/GammaM0
                if ExpVento[int(Ele_Lig[i])]=="Flat":
                    Veff2Rd[i]=0.5*fu_Lig[i]*Ant_block[i]/GammaM2+(1/np.sqrt(3))*fy_Lig[i]*Anv_block[i]/GammaM0
                else:
                    Veff2Rd[i]=0
            if Calc_Nu[i]=="Sim":
                if (ExpVento[int(Ele_Lig[i])]=="Flat") and (p2_gusset[i]==0):
                    if N_Bolt_Lig[i]==1:
                        NuRd[i]=2.0*(e2_gusset[i]-0.5*d_0[i])*t_gusset[i]*fu_Lig[i]/GammaM2
                    else:
                        if N_Bolt_Lig[i]==2:
                            if p1_gusset[i]<=2.5*d_0[i]:
                                beta_Cant[i]=0.4
                            elif p1_gusset[i]>=5*d_0[i]:
                                beta_Cant[i]=0.7
                            else:
                                beta_Cant[i]=0.4+(0.7-0.4)/(5*d_0[i]-2.5*d_0[i])*(p1_gusset[i]-2.5*d_0[i])
                        else:
                            if p1_gusset[i]<=2.5*d_0[i]:
                                beta_Cant[i]=0.5
                            elif p1_gusset[i]>=5*d_0[i]:
                                beta_Cant[i]=0.7
                            else:
                                beta_Cant[i]=0.5+(0.7-0.5)/(5*d_0[i]-2.5*d_0[i])*(p1_gusset[i]-2.5*d_0[i])
                        #Anet[i]=2*Dim[int(Ele_Lig[i])]*Esp[int(Ele_Lig[i])]-Esp[int(Ele_Lig[i])]**2-Esp[int(Ele_Lig[i])]*d_0[i]
                        Anet[i]=Area[int(Ele_Lig[i])]-Esp[int(Ele_Lig[i])]*d_0[i]
                        NuRd[i]=beta_Cant[i]*Anet[i]*fu_Lig[i]/GammaM2
                else:
                    NuRd[i]=0
            valores = [Veff1Rd[i],Veff2Rd[i],NuRd[i]]
            valores_sem_zero = [v for v in valores if v != 0]

            if valores_sem_zero:  # garante que não está vazio
                FLig_aux= min(valores_sem_zero)
            if np.max([Veff1Rd[i],Veff2Rd[i],NuRd[i]])==0:
                Ratio_Block_Lig[i]=0
            else:
                Ratio_Block_Lig[i]=Fele[int(Ele_Lig[i])]/FLig_aux

            Ratio_Lig[i]=np.nanmax([Ratio_bolt_lig[i],Ratio_Esmag_Lig[i],Ratio_Block_Lig[i]])
            if Ratio_Lig[i]==0:
                Check_Lig[i]="NA"
            elif Ratio_Lig[i]<1:
                Check_Lig[i]="OK"
            else:
                Check_Lig[i]="KO"

            valores = [FvRd[i], FbRd[i],Veff1Rd[i],Veff2Rd[i],NuRd[i]]
            valores_sem_zero = [v for v in valores if v != 0]

            if valores_sem_zero:  # garante que não está vazio
                FLig_Rd[i] = min(valores_sem_zero)
            else:
                FLig_Rd[i]  = 0  # ou None, dependendo do que quiseres
        else:
            Check_Lig[i]="Welded"
    # i_Bolt=Ratio_bolt_lig==0
    # i_Esmag=Ratio_Esmag_Lig==0
    # i_Block=Ratio_Block_Lig==0
    # Check_Bolt[i_Bolt]="NA"
    # Check_Esmag[i_Esmag]="NA"
    # Check_Block[i_Block]="NA"

    # i_Bolt=Ratio_bolt_lig<=1
    # i_Esmag=Ratio_Esmag_Lig<=1
    # i_Block=Ratio_Block_Lig<=1

    # Check_Bolt[i_Bolt]="OK"
    # Check_Esmag[i_Esmag]="OK"
    # Check_Block[i_Block]="OK"

    # i_Bolt=Ratio_bolt_lig>1
    # i_Esmag=Ratio_Esmag_Lig>1
    # i_Block=Ratio_Block_Lig>1

    # Check_Bolt[i_Bolt]="KO"
    # Check_Esmag[i_Esmag]="KO"
    # Check_Block[i_Block]="KO"
    Truss_Out=np.array(["Leg","Diagonal Bar","Horizontal Bar","External Manual Bar","Internal Manual Bar"])
    Lig_Out_csv=Connection_Bracing_Output(TrussType,Troco,nos_vento, Truss_Out,Ratio_Lig,t_gusset,e1_gusset,e2_gusset,p1_gusset,p2_gusset,fy_Lig,N_Bolt_Lig,Bolt_Steel,Bolt,Fele,FLig_Rd,Ele_Lig,Check_Lig,Calc_Block,Calc_Esmag,Calc_Bolt,Calc_Nu)
    Lig_Out_Exp=Connection_Bracing_Exp_Output(TrussType, Troco, nos_vento, Truss_Out, Ratio_Lig, t_gusset, e1_gusset, e2_gusset,
                                                p1_gusset, p2_gusset, fy_Lig, N_Bolt_Lig, Bolt_Steel, Bolt, Fele, FLig_Rd, Ele_Lig,
                                                d_Bolt, d_0, As_bolt, A_bolt, FvRd, Ratio_bolt_lig, alpha_d_end, alpha_d_inner,
                                                k1_edge, k1_inner, alpha_b_end, alpha_b_inner, FbRd_Edge_Inner, FbRd_Edge_End,
                                                FbRd, Anv_block, Ant_block, Veff1Rd, Veff2Rd, beta_Cant, NuRd)
    return Ratio_Lig,Lig_Out_csv,Lig_Out_Exp

def Bolted_Connection_Shaft(Shaft_Bolted_Matrix,TrussType,F_Tracao,F_Compressao,GammaM2,Sigma_u,Sigma,GammaM0,Area,Esp,Troco,nos_vento,divisor):
    ID_Bolted=len(Shaft_Bolted_Matrix)
    alpha_v_Bolt_B=np.zeros(ID_Bolted+1)
    Fele= np.maximum(np.abs(F_Tracao), np.abs(F_Compressao))
    As_bolt_B=np.zeros(ID_Bolted+1)
    A_bolt_B=np.zeros(ID_Bolted+1)
    d_0_B=np.zeros(ID_Bolted+1)
    d_Bolt_B=np.zeros(ID_Bolted+1)
    FvRd_B=np.zeros(ID_Bolted+1)
    Ratio_bolt_B=np.zeros(ID_Bolted+1)
    Ratio_Esmag_B=np.zeros(ID_Bolted+1)
    Ratio_Block_B=np.zeros(ID_Bolted+1)
    Check_Lig_B=np.zeros(ID_Bolted+1).astype(str)
    alpha_d_end_B=np.zeros(ID_Bolted+1)
    alpha_d_inner_B=np.zeros(ID_Bolted+1)
    k1_edge_B=np.zeros(ID_Bolted+1)
    k1_inner_B=np.zeros(ID_Bolted+1)
    alpha_b_end_B=np.zeros(ID_Bolted+1)
    Ratio_Lig_B=np.zeros(ID_Bolted+1)
    FLig_Rd_B=np.zeros(ID_Bolted+1)
    alpha_b_inner_B=np.zeros(ID_Bolted+1)
    Corte_Crit_B=np.zeros(ID_Bolted+1).astype(str)
    FbRd_B=np.zeros(ID_Bolted+1)
    FbRd_Edge_Inner_B=np.zeros(ID_Bolted+1)
    FbRd_Edge_End_B=np.zeros(ID_Bolted+1)
    Bolt_Steel_B=np.zeros(ID_Bolted+1)
    Veff1Rd_B=np.zeros(ID_Bolted+1)
    Lig_Type_B=np.zeros(ID_Bolted+1).astype(str)
    Veff2Rd_B=np.zeros(ID_Bolted+1)
    Anv_block_B=np.zeros(ID_Bolted+1)
    beta_Cant_B=np.zeros(ID_Bolted+1)
    t_B=np.zeros(ID_Bolted+1)
    Ant_block_B=np.zeros(ID_Bolted+1)
    fub_B=np.zeros(ID_Bolted+1)
    Bolt_Steel_B=np.zeros(ID_Bolted+1).astype(str)
    fyb_B=np.zeros(ID_Bolted+1)
    NuRd_B=np.zeros(ID_Bolted+1)
    Anet_B=np.zeros(ID_Bolted+1)
    fu_B,fy_B=np.zeros(ID_Bolted+1),np.zeros(ID_Bolted+1)
    N_Bolt_B=np.zeros(ID_Bolted+1)
    Node_B=np.zeros(ID_Bolted+1)
    Steel_B=np.zeros(ID_Bolted+1).astype(str)
    Corte_B=np.ones(ID_Bolted+1)
    Bolt_B=np.ones(ID_Bolted+1).astype(str)
    p1_B=np.zeros(ID_Bolted+1)
    e2_B=np.zeros(ID_Bolted+1)
    Bolt_B=np.zeros(ID_Bolted+1).astype(str)
    Calc_Bolt_B=np.zeros(ID_Bolted+1).astype(str)
    Calc_Esmag_B=np.zeros(ID_Bolted+1).astype(str)
    Calc_Block_B=np.zeros(ID_Bolted+1).astype(str)
    e1_B=np.zeros(ID_Bolted+1)
    p2_B=np.zeros(ID_Bolted+1)
    Check_Bolt_B=np.zeros(ID_Bolted+1).astype(str)
    Check_Esmag_B=np.zeros(ID_Bolted+1).astype(str)
    Check_Block_B=np.zeros(ID_Bolted+1).astype(str)
    Calc_Nu_B=np.zeros(ID_Bolted+1).astype(str)
    F_Ele_B=np.zeros(ID_Bolted+1)
    Check_DisMin_B=np.zeros(ID_Bolted+1).astype(str)
    if len(Shaft_Bolted_Matrix)==0:
        Shaft_Bolted_Matrix=np.empty((0,12))
    N_Bolt_B[1:]=Shaft_Bolted_Matrix[:,9]
    Node_B[1:]=Shaft_Bolted_Matrix[:,1]
    Bolt_Steel_B[1:]=Shaft_Bolted_Matrix[:,11]
    Corte_B[1:]=Shaft_Bolted_Matrix[:,2].astype(float)
    p1_B[1:]=Shaft_Bolted_Matrix[:,7].astype(float)*10**-3
    e1_B[1:]=Shaft_Bolted_Matrix[:,5].astype(float)*10**-3
    t_B[1:]=Shaft_Bolted_Matrix[:,4].astype(float)*10**-3
    e2_B[1:]=Shaft_Bolted_Matrix[:,6].astype(float)*10**-3
    p2_B[1:]=Shaft_Bolted_Matrix[:,8].astype(float)*10**-3
    Bolt_B[1:]=Shaft_Bolted_Matrix[:,10].astype(str)
    Steel_B[1:]=Shaft_Bolted_Matrix[:,3]



    for i in range(1,ID_Bolted+1):
        ELE_B=elementos_do_no_filt(Node_B[i],TrussType,"Leg")
        ELE_B=np.max(ELE_B)
        F_Ele_B[i]=Fele[int(ELE_B)]
        Calc_Bolt_B[i]="Sim"
        Calc_Esmag_B[i]="Sim"
        Calc_Block_B[i]="Sim"
        Calc_Nu_B[i]="Sim"
        if e1_B[i]==0 or e2_B[i]==0:
            Calc_Bolt_B[i]="Nao"
            Calc_Esmag_B[i]="Nao"
            Calc_Block_B[i]="Nao"
            Calc_Nu_B[i]="Nao"
        
        if Bolt_B[i]=="0" or Bolt_Steel_B[i]=="0" or N_Bolt_B[i]=="0":
            Calc_Bolt_B[i]="Nao"
            Calc_Esmag_B[i]="Nao"
            Calc_Block_B[i]="Nao"
            Calc_Nu_B[i]="Nao"
        if Bolt_B[i]!="0" and Bolt_Steel_B[i]!="0" and N_Bolt_B[i]!="0" and Corte_B[i]!="0" :
            Calc_Bolt_B[i]="Sim"
        if Bolt_B[i]!="0" and N_Bolt_B[i]!="0" and Corte_B[i]!="0":
            Calc_Bolt_B[i]="Sim"
            Bolt_Steel_B[i]!="8.8"

        if Calc_Bolt_B[i]=="Sim":
            ###### Verificação parafuso ######
            fub_B[i],fyb_B[i]=get_fu_fy(Bolt_Steel_B[i])
            fub_B[i]=fub_B[i]*10**6
            fyb_B[i]=fyb_B[i]*10**6
            alpha_v_Bolt_B[i]=0.6
            d_Bolt_B[i], d_0_B[i], As_bolt_B[i], A_bolt_B[i]=get_bolt_props(Bolt_B[i])
            FvRd_B[i]=Corte_B[i]*alpha_v_Bolt_B[i]*fub_B[i]*As_bolt_B[i]*N_Bolt_B[i]/GammaM2
            Ratio_bolt_B[i]=F_Ele_B[i]/FvRd_B[i]

        if Calc_Esmag_B[i]=="Sim":       
            ##### Validação distancias minimas
            if e1_B[i]!=0 and e1_B[i]<=1.2*d_0_B[i]:
                Check_DisMin_B[i]="KO"
            elif e2_B[i]!=0 and e2_B[i]<=1.2*d_0_B[i]:
                Check_DisMin_B[i]="KO"
            elif p1_B[i]!=0 and p1_B[i]<=2.2*d_0_B[i]:
                Check_DisMin_B[i]="KO"
            elif p2_B[i]!=0 and p2_B[i]<=2.4*d_0_B[i]:
                Check_DisMin_B[i]="KO"
            else:
                Check_DisMin_B[i]="OK" 
            ##### Verificação corte #######
            if t_B[i]==0:
                t_B[i]=Esp[int(ELE_B)]
            else:
                t_B[i]=np.min([Esp[int(ELE_B)],t_B[i]])
            #if Leg_Type=="Flat" and ExpVento[int(Ele_Lig[i])]=="Flat":
            #    ####Cantoneira-Cantoneira
            #    fu_Lig[i],fy_Lig[i]=Sigma_u[int(Ele_Lig[i])],Sigma[int(Ele_Lig[i])]
            #else: 
            ########## Montantes não aparafusados um ao outro, 
            if Steel_B[i]!="0":
                fu_B[i],fy_B[i]=get_fu_fy(Steel_B[i])
                fu_B[i]=fu_B[i]*10**6
                fy_B[i]=fy_B[i]*10**6
            else:
                fu_B[i]=Sigma_u[int(ELE_B)]
                fy_B[i]=Sigma[int(ELE_B)]
            if ((p1_B[i]!=0) and (p2_B[i]!=0)):
                alpha_d_end_B[i]=e1_B[i]/(3*d_0_B[i])
                alpha_d_inner_B[i]=p1_B[i]/(3*d_0_B[i])-1/4
                k1_edge_B[i]=np.min([2.8*e2_B[i]/d_0_B[i]-1.7,1.4*p2_B[i]/d_0_B[i]-1.7,2.5])
                k1_inner_B[i]=np.min([1.4*p2_B[i]/d_0_B[i]-1.7,2.5])     
                alpha_b_end_B[i]=np.min([alpha_d_end_B[i],fub_B[i]/fu_B[i],1.0])
                alpha_b_inner_B[i]=np.min([alpha_d_inner_B[i],fub_B[i]/fu_B[i],1.0])
                #Aux=[k1_edge[i]*alpha_b_end[i],k1_edge[i]*alpha_b_inner[i]]
                #iAux=np.argmin([k1_edge[i]*alpha_b_end[i],k1_edge[i]*alpha_b_inner[i]])
                #Case_Corte_Name=np.array(["Edge-End","Inner-End"])
                #Corte_Crit[i]=Case_Corte_Name[iAux]
                FbRd_Edge_Inner_B[i]=Corte_B[i]*alpha_b_inner_B[i]* k1_edge_B[i]*fu_B[i]*d_Bolt_B[i]*t_B[i]/GammaM2
                FbRd_Edge_End_B[i]=Corte_B[i]*alpha_b_end_B[i]* k1_edge_B[i]*fu_B[i]*d_Bolt_B[i]*t_B[i]/GammaM2            
                FbRd_B[i]=FbRd_Edge_Inner_B[i]*(N_Bolt_B[i]-2)+FbRd_Edge_End_B[i]*(2)
            elif ((p2_B[i]==0) and (p1_B[i]!=0)):
                alpha_d_end_B[i]=e1_B[i]/(3*d_0_B[i])
                alpha_d_inner_B[i]=p1_B[i]/(3*d_0_B[i])-1/4
                k1_edge_B[i]=np.min([2.8*e2_B[i]/d_0_B[i]-1.7,2.5])
                k1_inner_B[i]=np.min([2.5])     
                alpha_b_end_B[i]=np.min([alpha_d_end_B[i],fub_B[i]/fu_B[i],1.0])
                alpha_b_inner_B[i]=np.min([alpha_d_inner_B[i],fub_B[i]/fu_B[i],1.0])
                FbRd_Edge_Inner_B[i]=Corte_B[i]*alpha_b_inner_B[i]* k1_edge_B[i]*fu_B[i]*d_Bolt_B[i]*t_B[i]/GammaM2
                FbRd_Edge_End_B[i]=Corte_B[i]*alpha_b_end_B[i]* k1_edge_B[i]*fu_B[i]*d_Bolt_B[i]*t_B[i]/GammaM2            
                FbRd_B[i]=FbRd_Edge_Inner_B[i]*(N_Bolt_B[i]-1)+FbRd_Edge_End_B[i]*(1)
            elif ((p1_B[i]==0) and (p2_B[i]==0)):
                alpha_d_end_B[i]=e1_B[i]/(3*d_0_B[i])
                k1_edge_B[i]=np.min([2.8*e2_B[i]/d_0_B[i]-1.7,2.5]) 
                alpha_b_end_B[i]=np.min([alpha_d_end_B[i],fub_B[i]/fu_B[i],1.0])
                FbRd_Edge_Inner_B[i]=0
                FbRd_Edge_End_B[i]=Corte_B[i]*alpha_b_end_B[i]* k1_edge_B[i]*fu_B[i]*d_Bolt_B[i]*t_B[i]/GammaM2  
                FbRd_B[i]=FbRd_Edge_End_B[i]
            elif ((p2_B[i]!=0) and (p1_B[i]==0)):
                alpha_d_end_B[i]=e1_B[i]/(3*d_0_B[i])
                k1_edge_B[i]=np.min([2.8*e2_B[i]/d_0_B[i]-1.7,1.4*p2_B[i]/d_0_B[i]-1.7,2.5])    
                alpha_b_end_B[i]=np.min([alpha_d_end_B[i],fub_B[i]/fu_B[i],1.0])
                FbRd_Edge_Inner_B[i]=0
                FbRd_Edge_End_B[i]=Corte_B[i]*alpha_b_end_B[i]* k1_edge_B[i]*fu_B[i]*d_Bolt_B[i]*t_B[i]/GammaM2 
                FbRd_B[i]=FbRd_Edge_End_B[i]*N_Bolt_B[i]
            Ratio_Esmag_B[i]=F_Ele_B[i]/FbRd_B[i]
        if Calc_Block_B[i]=="Sim":                        
            ##### Verificação Bloco #######
            if p2_B[i]!=0:
                # if p1_gusset[i]==0:
                #     Anv_block[i]=t_gusset[i]*(e1_gusset[i]-0.5*d_0[i])
                #     Ant_block[i]=t_gusset[i]*(e2_gusset[i]+(N_Bolt_Lig[i]-1)*p2_gusset[i]-((N_Bolt_Lig[i]-0.5)*d_0[i]))
                # elif p2_gusset[i]==0:
                #     Anv_block[i]=t_gusset[i]*(e1_gusset[i]+(N_Bolt_Lig[i]-1)*p1_gusset[i]-((N_Bolt_Lig[i]-0.5)*d_0[i]))
                #     Ant_block[i]=t_gusset[i]*(e2_gusset[i]-0.5*d_0[i])
                # elif (p1_gusset[i]!=0) and (p1_gusset[i]!=0):
                #     Anv_block[i]=t_gusset[i]*(e1_gusset[i]+(N_Bolt_Lig[i]/2-1)*p1_gusset[i]-((N_Bolt_Lig[i]/2-0.5)*d_0[i]))            
                #     Ant_block[i]=t_gusset[i]*(e2_gusset[i]+(N_Bolt_Lig[i]/2-1)*p2_gusset[i]-((N_Bolt_Lig[i]/2-0.5)*d_0[i]))
                Anv_block_B[i]= t_B[i]*(e1_B[i]+ (N_Bolt_B[i]-2)/2*p1_B[i] - (N_Bolt_B[i]/2-0.5)*d_0_B[i])*2
                Ant_block_B[i]= t_B[i]*(p2_B[i]-d_0_B[i])
            else:
                Anv_block_B[i]=0            
                Ant_block_B[i]=0          
            Veff1Rd_B[i]=Corte_B[i]*fu_B[i]*Ant_block_B[i]/GammaM2+(1/np.sqrt(3))*fy_B[i]*Anv_block_B[i]/GammaM0
            Veff2Rd_B[i]=Corte_B[i]*0.5*fu_B[i]*Ant_block_B[i]/GammaM2+(1/np.sqrt(3))*fy_B[i]*Anv_block_B[i]/GammaM0
        if Calc_Nu_B[i]=="Sim":
            if p2_B[i]==0:
                if N_Bolt_B[i]==1:
                    NuRd_B[i]=2.0*(e2_B[i]-0.5*d_0_B[i])*t_B[i]*fu_B[i]/GammaM2
                else:
                    if N_Bolt_B[i]==2:
                        if p1_B[i]<=2.5*d_0_B[i]:
                            beta_Cant_B[i]=0.4
                        elif p1_B[i]>=5*d_0_B[i]:
                            beta_Cant_B[i]=0.7
                        else:
                            beta_Cant_B[i]=0.4+(0.7-0.4)/(5*d_0_B[i]-2.5*d_0_B[i])*(p1_B[i]-2.5*d_0_B[i])
                    else:
                        if p1_B[i]<=2.5*d_0_B[i]:
                            beta_Cant_B[i]=0.5
                        elif p1_B[i]>=5*d_0_B[i]:
                            beta_Cant_B[i]=0.7
                        else:
                            beta_Cant_B[i]=0.5+(0.7-0.5)/(5*d_0_B[i]-2.5*d_0_B[i])*(p1_B[i]-2.5*d_0_B[i])
                    #Anet_B[i]=2*Dim[int(Ele_Lig[i])]*Esp[int(Ele_Lig[i])]-Esp[int(Ele_Lig[i])]**2-Esp[int(Ele_Lig[i])]*d_0_B[i]
                    Anet_B[i]=Area[int(ELE_B)]-Esp[int(ELE_B)]*d_0_B[i]
                    NuRd_B[i]=Corte_B[i]*beta_Cant_B[i]*Anet_B[i]*fu_B[i]/GammaM2
            else:
                NuRd_B[i]=0
        valores = [Veff1Rd_B[i],Veff2Rd_B[i],NuRd_B[i]]
        valores_sem_zero = [v for v in valores if v != 0]

        if valores_sem_zero:  # garante que não está vazio
            FLig_aux_B= min(valores_sem_zero)
        if np.max([Veff1Rd_B[i],Veff2Rd_B[i],NuRd_B[i]])==0:
            Ratio_Block_B[i]=0
        else:
            Ratio_Block_B[i]=F_Ele_B[i]/FLig_aux_B

        Ratio_Lig_B[i]=np.nanmax([Ratio_bolt_B[i],Ratio_Esmag_B[i],Ratio_Block_B[i]])
        if Ratio_Lig_B[i]==0:
            Check_Lig_B[i]="NA"
        elif Ratio_Lig_B[i]<1:
            Check_Lig_B[i]="OK"
        else:
            Check_Lig_B[i]="KO"

        valores = [FvRd_B[i], FbRd_B[i],Veff1Rd_B[i],Veff2Rd_B[i],NuRd_B[i]]
        valores_sem_zero = [v for v in valores if v != 0]

        if valores_sem_zero:  # garante que não está vazio
            FLig_Rd_B[i] = min(valores_sem_zero)
        else:
            FLig_Rd_B[i]  = 0  # ou None, dependendo do que quiseres

    B_Out_csv=Connection_Shaft_Output(TrussType,Troco,Ratio_Lig_B,nos_vento,divisor,t_B,e1_B,e2_B,p1_B,p2_B,Corte_B,fy_B,N_Bolt_B,Bolt_Steel_B,Bolt_B,F_Ele_B,FLig_Rd_B,Check_Lig_B,Calc_Bolt_B,Calc_Esmag_B,Calc_Block_B,Calc_Nu_B)
    B_Out_Exp=Connection_Shaft_Exp_Output(TrussType, Troco, Ratio_Lig_B, nos_vento, divisor,
                                            t_B, e1_B, e2_B, p1_B, p2_B, Corte_B,
                                            fy_B, N_Bolt_B, Bolt_Steel_B, Bolt_B,
                                            F_Ele_B, FLig_Rd_B,
                                            d_Bolt_B, d_0_B, As_bolt_B, A_bolt_B, FvRd_B,
                                            Ratio_bolt_B, alpha_d_end_B, alpha_d_inner_B,
                                            k1_edge_B, k1_inner_B, alpha_b_end_B, alpha_b_inner_B,
                                            FbRd_Edge_Inner_B, FbRd_Edge_End_B, FbRd_B,
                                            Anv_block_B, Ant_block_B, Veff1Rd_B, Veff2Rd_B,
                                            beta_Cant_B, NuRd_B
                                        )
    return Ratio_Lig_B,B_Out_csv,B_Out_Exp

def Flanged_Bolted_Connection(Shaft_Flange_Matrix,TrussType,F_Tracao,Dim,Esp,Sigma_u,Sigma,GammaM2,GammaM0,Area,nos_vento,divisor):
    if Shaft_Flange_Matrix.size == 0:
        Shaft_Flange_Matrix=np.empty((0,15))
    ID_Flanged=len(Shaft_Flange_Matrix)

    Node_F=np.zeros(ID_Flanged+1)
    Ext_Per_F=np.zeros(ID_Flanged+1)
    D_between_Bolt_F=np.zeros(ID_Flanged+1)
    t_F=np.zeros(ID_Flanged+1)
    Steel_F=np.zeros(ID_Flanged+1).astype(str)
    Bolt_F=np.zeros(ID_Flanged+1).astype(str)
    N_Bolt_F=np.zeros(ID_Flanged+1)
    Bolt_Steel_F=np.zeros(ID_Flanged+1).astype(str)
    Stiffeners_F=np.zeros(ID_Flanged+1).astype(str)
    N_Stiff_F=np.zeros(ID_Flanged+1)
    t_Stiff_F=np.zeros(ID_Flanged+1)
    L_Stiff_F=np.zeros(ID_Flanged+1)
    h_Stiff_F=np.zeros(ID_Flanged+1)
    theta_F=np.zeros(ID_Flanged+1)
    D_Ext_F=np.zeros(ID_Flanged+1)
    D_Bolt_Flange_F=np.zeros(ID_Flanged+1)
    D_int_F=np.zeros(ID_Flanged+1)
    e1_F=np.zeros(ID_Flanged+1)
    e2_F=np.zeros(ID_Flanged+1)
    alpha_r0_F=np.zeros(ID_Flanged+1)
    d_Bolt_F=np.zeros(ID_Flanged+1)
    d_0_F=np.zeros(ID_Flanged+1)
    As_bolt_F=np.zeros(ID_Flanged+1)
    A_bolt_F=np.zeros(ID_Flanged+1)
    Lb_F=np.zeros(ID_Flanged+1)
    leff=np.zeros(ID_Flanged+1)
    alpha_r_F=np.zeros(ID_Flanged+1)
    epsilon_F=np.zeros(ID_Flanged+1)
    F_Ele_F=np.zeros(ID_Flanged+1)
    fub_F=np.zeros(ID_Flanged+1)
    fyb_F=np.zeros(ID_Flanged+1)  
    fu_F=np.zeros(ID_Flanged+1)  
    fy_F=np.zeros(ID_Flanged+1)  
    n_F=np.zeros(ID_Flanged+1)  


    elinha_F=np.zeros(ID_Flanged+1)
    m_F=np.zeros(ID_Flanged+1)
    FtRd=np.zeros(ID_Flanged+1)
    Rb_F=np.zeros(ID_Flanged+1)
    R_F=np.zeros(ID_Flanged+1)
    k1_F=np.zeros(ID_Flanged+1)
    mplRd_F=np.zeros(ID_Flanged+1)
    NT1Rd=np.zeros(ID_Flanged+1)
    Re_F=np.zeros(ID_Flanged+1)
    k3_F=np.zeros(ID_Flanged+1)
    NT2Rd=np.zeros(ID_Flanged+1)
    NT3Rd=np.zeros(ID_Flanged+1)
    NT4Rd=np.zeros(ID_Flanged+1)
    NT5Rd=np.zeros(ID_Flanged+1)
    Flag_Data_F=np.zeros(ID_Flanged+1)
    NTRd=np.zeros(ID_Flanged+1)
    aws_F=np.zeros(ID_Flanged+1)
    t_Tube_F=np.zeros(ID_Flanged+1)
    BpRd=np.zeros(ID_Flanged+1)
    dm_F=np.zeros(ID_Flanged+1)
    Ratio_F=np.zeros(ID_Flanged+1)
    Check_F=np.zeros(ID_Flanged+1).astype(str)

    Node_F[1:]=Shaft_Flange_Matrix[:,1].astype(int)
    Ext_Per_F[1:]=Shaft_Flange_Matrix[:,2].astype(float)*10**-3
    D_between_Bolt_F[1:]=Shaft_Flange_Matrix[:,3].astype(float)*10**-3
    t_F[1:]=Shaft_Flange_Matrix[:,4].astype(float)*10**-3
    Steel_F[1:]=Shaft_Flange_Matrix[:,5]
    Bolt_F[1:]=Shaft_Flange_Matrix[:,6]
    N_Bolt_F[1:]=Shaft_Flange_Matrix[:,7].astype(float)
    Bolt_Steel_F[1:]=Shaft_Flange_Matrix[:,8]
    Stiffeners_F[1:]=Shaft_Flange_Matrix[:,9]
    N_Stiff_F[1:]=Shaft_Flange_Matrix[:,10]
    t_Stiff_F[1:]=Shaft_Flange_Matrix[:,11].astype(float)*10**-3
    L_Stiff_F[1:]=Shaft_Flange_Matrix[:,12].astype(float)*10**-3
    h_Stiff_F[1:]=Shaft_Flange_Matrix[:,13].astype(float)*10**-3

    Bride=np.zeros(ID_Flanged+1).astype(str)


    # theta_F[1:]=360/N_Bolt_F[1:]
    # D_Ext_F[1:]=Ext_Per_F[1:]/np.pi
    # D_Bolt_Flange_F[1:]=2*D_between_Bolt_F[1:]/2/np.sin(np.radians(theta_F[1:]/2))

    for i in range(1,ID_Flanged+1):
        ###Assume-se Bride Creuse 
        if Bolt_F[i]!="0" and Ext_Per_F[i]!=0 and D_between_Bolt_F[i]!=0 and t_F[i]!=0 and  Steel_F[i]!="0" and Bolt_F[i]!="0" and N_Bolt_F[i]!=0 and Bolt_Steel_F[i]!="0" and Stiffeners_F[i]!="0":
            theta_F[i]=360/N_Bolt_F[i]
            D_Ext_F[i]=Ext_Per_F[i]/np.pi
            D_Bolt_Flange_F[i]=2*D_between_Bolt_F[i]/2/np.sin(np.radians(theta_F[i]/2))
            ELE_F=elementos_do_no_filt(Node_F[i],TrussType,"Leg")
            # if len(ELE_F)==1:
            #     ELE_F=np.array([ELE_F[0],ELE_F[0]])
            # print(ELE_F)
            if Stiffeners_F[i]=="no":
                ELE_F_Baixo=np.min(ELE_F)
                ELE_F=np.max(ELE_F)
                F_Ele_F[i]=F_Tracao[int(ELE_F_Baixo)]
                D_int_F[i]=Dim[int(ELE_F)]-Esp[int(ELE_F)]
                e2_F[i]=(D_Ext_F[i]-D_Bolt_Flange_F[i])/2
                e1_F[i]=(D_Bolt_Flange_F[i]-D_int_F[i])/2
                t_Tube_F[i]=Esp[int(ELE_F)]
                if e1_F[i] < 0 or e2_F[i] < 0:
                    Ratio_F[i] = 0
                    Flag_Data_F[i]=1
                    Check_F[i]="NA"
                    continue
                alpha_r0_F[i]=(e2_F[i]/e1_F[i]+1)/(e2_F[i]/e1_F[i])**3
                d_Bolt_F[i], d_0_F[i], As_bolt_F[i], A_bolt_F[i]=get_bolt_props(Bolt_F[i])
                Lb_F[i]=2*t_F[i]+d_Bolt_F[i]*0.8
                leff[i]=2*np.pi*D_int_F[i]/2
                alpha_r_F[i]=4*(e1_F[i]/t_F[i])**3*As_bolt_F[i]/(Lb_F[i]*leff[i])

                epsilon_F[i]=e2_F[i]*np.min([1,np.sqrt(alpha_r0_F[i]/alpha_r_F[i])])
                fub_F[i],fyb_F[i]=get_fu_fy(Bolt_Steel_F[i])
                fub_F[i]=fub_F[i]*10**6
                fyb_F[i]=fyb_F[i]*10**6   
                if Steel_F[i]!="0":
                    fu_F[i],fy_F[i]=get_fu_fy(Steel_F[i])
                    fu_F[i]=fu_F[i]*10**6
                    fy_F[i]=fy_F[i]*10**6
                else:
                    fu_F[i]=Sigma_u[int(ELE_F[0])]
                    fy_F[i]=Sigma[int(ELE_F[0])]
                n_F[i]=np.min([(2*e2_F[i]+epsilon_F[i])/3,epsilon_F[i]+0.74*t_F[i]])
                aws_F[i]=0.7*t_Tube_F[i]
                m_F[i]=e1_F[i]-t_Tube_F[i]/2-0.8*np.sqrt(2)*aws_F[i]
                Rb_F[i]=D_Bolt_Flange_F[i]/2
                R_F[i]=D_int_F[i]/2
                if (D_Ext_F[i]/2)**2-m_F[i]**2<0:
                    NTRd[i]=0
                    Flag_Data_F[i]=1
                else:
                    elinha_F[i]=np.sqrt((D_Ext_F[i]/2)**2-m_F[i]**2)-Rb_F[i]
                    FtRd[i]=0.9*As_bolt_F[i]*fub_F[i]/GammaM2


                    k1_F[i]=np.log((Rb_F[i])/(D_int_F[i]/2))
                    mplRd_F[i]=(t_F[i]**2*fy_F[i])/(4*GammaM0)
                    

                    Re_F[i]= Rb_F[i]+n_F[i]
                    k3_F[i]=np.log(Re_F[i]/R_F[i])
                    Bride[i]="Creuse"
                    if Bride[i]=="Creuse":
                        NT1Rd[i]=np.min([2*np.pi*mplRd_F[i]*(1+(1/k1_F[i])),2*np.pi*mplRd_F[i]*N_Bolt_F[i]*np.min([2,1+2/np.pi*elinha_F[i]/m_F[i]])])
                        NT2Rd[i]=2*np.pi*mplRd_F[i]+N_Bolt_F[i]*FtRd[i]*(1-(k1_F[i]/k3_F[i]))
                    elif Bride[i]=="Pleine":
                        NT1Rd[i]=np.min([2*np.pi()*mplRd_F[i]*(1+(2/k1_F[i])),2*np.pi*mplRd_F[i]*N_Bolt_F[i]*np.min([2,2/np.pi*elinha_F[i]/m_F[i]])])
                        NT2Rd[i]=2*np.pi()*mplRd_F[i]*(1+(1/k3_F[i]))+N_Bolt_F[i]*FtRd[i]*(1-(k1_F[i]/k3_F[i]))
                    NT3Rd[i]=N_Bolt_F[i]*FtRd[i]
                    NT4Rd[i]=Area[int(ELE_F)]*Sigma[int(ELE_F)]/GammaM0
                    NTRd[i]=np.min([NT1Rd[i],NT2Rd[i],NT3Rd[i],NT4Rd[i]])
            else:
                fub_F[i],fyb_F[i]=get_fu_fy(Bolt_Steel_F[i])
                fub_F[i]=fub_F[i]*10**6
                fyb_F[i]=fyb_F[i]*10**6   
                if Steel_F[i]!="0":
                    fu_F[i],fy_F[i]=get_fu_fy(Steel_F[i])
                    fu_F[i]=fu_F[i]*10**6
                    fy_F[i]=fy_F[i]*10**6
                else:
                    fu_F[i]=Sigma_u[int(ELE_F[0])]
                    fy_F[i]=Sigma[int(ELE_F[0])]
                d_Bolt_F[i], d_0_F[i], As_bolt_F[i], A_bolt_F[i]=get_bolt_props(Bolt_F[i])
                FtRd[i]=0.9*As_bolt_F[i]*fub_F[i]/GammaM2
                NT3Rd[i]=N_Bolt_F[i]*FtRd[i]
                NT4Rd[i]=Area[int(ELE_F)]*Sigma[int(ELE_F)]/GammaM0
                dm_F[i]=nut_across_flats(Bolt_F[i])
                BpRd[i]=0.6*np.pi*dm_F[i]*t_F[i]*fu_F[i]/GammaM2
                NTRd[i]=np.min([NT3Rd[i],NT4Rd[i],BpRd[i]])
            if NTRd[i]!=0:
                Ratio_F[i]=F_Ele_F[i]/NTRd[i]
            else:
                Ratio_F[i]=0
            if Ratio_F[i]>1:
                Check_F[i]="K0"
            elif Ratio_F[i]<1 and Flag_Data_F[i]==0:
                Check_F[i]="OK"
            else:
                Check_F[i]="NA"
        else:
            Ratio_F[i]=0
            Check_F[i]="NA"
    Check_F[0]="NA"
    F_Out = np.column_stack((Ext_Per_F,D_between_Bolt_F,t_F,Steel_F,Bolt_F,N_Bolt_F,Bolt_Steel_F,Stiffeners_F,N_Stiff_F,t_Stiff_F,L_Stiff_F,h_Stiff_F,D_Ext_F,D_Bolt_Flange_F,D_int_F,t_Tube_F,aws_F, e1_F,e2_F,alpha_r0_F,d_Bolt_F,Lb_F,leff,alpha_r_F,epsilon_F,fub_F,fyb_F,n_F,m_F,elinha_F,FtRd,Rb_F,R_F,Re_F,k3_F,NT1Rd,NT2Rd,NT3Rd,NT4Rd,BpRd,NTRd,Ratio_F,Check_F))
    F_Out_csv=Flanged_Shaft_Output(Ratio_F,Check_F,nos_vento,divisor,D_Ext_F,t_F,Steel_F,D_Bolt_Flange_F,N_Bolt_F,Bolt_F,Bolt_Steel_F,NTRd,F_Ele_F)
    F_Out_Exp=Flanged_Shaft_Exp_Output(Ratio_F,Check_F,nos_vento, divisor,
                         N_Bolt_F, Bolt_F, Bolt_Steel_F, NTRd, F_Ele_F,
                         D_Ext_F, t_F, Steel_F,N_Bolt_F ,D_Bolt_Flange_F,
                         e2_F, e1_F, alpha_r0_F, d_Bolt_F, d_0_F,
                         As_bolt_F, A_bolt_F, Lb_F, leff, alpha_r_F,
                         epsilon_F, n_F, aws_F, m_F, Rb_F, R_F,
                         elinha_F, FtRd, k1_F, mplRd_F, Re_F, k3_F,
                         Bride, NT1Rd, NT2Rd, NT3Rd, NT4Rd)
    return Ratio_F,F_Out_csv,F_Out_Exp


def Flanged_Foundation_Angle(Flange_Foundation_Angle,Trusstype,F_tracao):

    if Flange_Foundation_Angle[0,0]=="No":
        i_Montante=Trusstype=="Leg"
        F_montante=np.max(F_tracao[i_Montante])
        _,_,Abolt,_=get_bolt_props(Flange_Foundation_Angle[0,4])
        Nbolt=int(Flange_Foundation_Angle[0,5])
        fub, fyb=get_fu_fy(Flange_Foundation_Angle[0,6])
        fyb=fyb*10**6
        Nrd=fyb*Nbolt*Abolt
        ratio=F_montante/Nrd
        if ratio<1:
            Check="OK",
        else:
            Check="KO"
        row = np.array([[
                    "Embedded",
                    "Flange Widht (mm)",
                    "Flange Thickness (mm)",
                    "Flange Steel",
                    "Bolt Diameter",
                    "Number of bolts",
                    "Bolt Steel",
                    "Ned (kN)",
                    "Nrd (kN)",
                    "Ratio",
                    "Check"
        ],  
                    [
                    Flange_Foundation_Angle[0,0],
                    Flange_Foundation_Angle[0,1],
                    Flange_Foundation_Angle[0,2],
                    Flange_Foundation_Angle[0,3],
                    Flange_Foundation_Angle[0,4],
                    Flange_Foundation_Angle[0,5],
                    Flange_Foundation_Angle[0,6],
                    np.round(F_montante*10**-3,2),
                    np.round(Nrd*10**-3,2),
                    np.round(ratio,2),
                    Check]], dtype=object)

    elif Flange_Foundation_Angle[0,0]=="Yes":
        ratio=0
        row = np.array([[
                    "Embedded",
                    "Flange Widht (mm)",
                    "Flange Thickness (mm)",
                    "Flange Steel",
                    "Bolt Diameter",
                    "Number of bolts",
                    "Bolt Steel",
                    "Ned (kN)",
                    "Nrd (kN)",
                    "Ratio",
                    "Check"
        ],  
                    [
                    Flange_Foundation_Angle[0,0],
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--"]], dtype=object)        
    else:
        ratio=0
        row = np.array([[
                    "Embedded",
                    "Flange Widht (mm)",
                    "Flange Thickness (mm)",
                    "Flange Steel",
                    "Bolt Diameter",
                    "Number of bolts",
                    "Bolt Steel",
                    "Ned (kN)",
                    "Nrd (kN)",
                    "Ratio",
                    "Check"
        ],  
                    [
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--",
                    "--"]], dtype=object)        

    return ratio, row