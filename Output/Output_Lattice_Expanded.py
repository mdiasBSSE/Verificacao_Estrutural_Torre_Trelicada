import numpy as np
import pandas as pd
def Expanded_Output(Name_Case,nos_vento,):
    Expanded_out="Expanded Output"
    for i in range(1,len(Name_Case)):
        Expanded_Output.append(f"Case - {Name_Case[i]}")
        Expanded_Output.append(" ")
        for j in range(nos_vento):
            Name_Enc_Out = f"Shaft { nos_vento.shape[0] -i}"
            row = np.array([
                Name_Enc_Out,

            ], dtype=object)

            Enc_Out_csv = np.vstack((Enc_Out_csv, row))


def Buckling_Exp_Output(Class_Enc,nos_vento,TrussType,Truss_Out,Size_Profile,Elements_Matrix,Lambda_Enc_eff,Lambda_1,NTracRd,NbRd,Troco,Ratio_Enc,F_Compressao,Comprimento_Enc,Ratio_Class,Lim_Class,Lambda_Lim,k_Enc,alpha_Enc,Phi_Enc,Chi_Enc,Aeff,Red,eta_Enc,Lambda_Check):
    Enc_Out_csv=np.array(["Member", "Profile","Buckling Lenght", "Steel", "Ratio Class", "Lim Class", "Section Class","Lambda 1","Lambda Lim", "Lambda Check", "Slenderness",
                          "k","alpha","Phi","Chi","Aeff","Red","Eta","NcRd","NbRd","Ned","Ratio"])
    ID_Ele_vec= np.arange(0, len(Class_Enc) + 1)

    for j in range(3):
        for i in range(1, nos_vento.shape[0]):
            # Agrupa quando j==2 ("Horizontal Bar")
            if j == 2:
                mask_type = (TrussType == "Horizontal Bar") | (TrussType == "External Manual Bar") | (TrussType == "Internal Manual Bar")
            else:
                mask_type = (TrussType == Truss_Out[j])

            i_Out_Zone = (Troco == i) & mask_type
            indices = np.where(i_Out_Zone)[0]

            if len(indices) != 0:
                iCrit_local = np.argmax(Ratio_Enc[indices])
                iCrit_abs = indices[iCrit_local]


                # nome do grupo de saída
                Name_Enc_Out = f"{Truss_Out[j]} - Shaft { nos_vento.shape[0] -i}"
                row = np.array([
                    Name_Enc_Out,
                    Size_Profile[iCrit_abs],
                    round(Comprimento_Enc[iCrit_abs]*1000,3),
                    Elements_Matrix[iCrit_abs - 1, 7],
                    round(Ratio_Class[iCrit_abs],3),
                    round(Lim_Class[iCrit_abs],3),
                    int(Class_Enc[iCrit_abs]),
                    round(Lambda_1[iCrit_abs],3),
                    int(Lambda_Lim[iCrit_abs]),
                    Lambda_Check[iCrit_abs],
                    round(Lambda_Enc_eff[iCrit_abs] * Lambda_1[iCrit_abs], 1),
                    round(k_Enc[iCrit_abs],3),
                    round(alpha_Enc[iCrit_abs],3),
                    round(Phi_Enc[iCrit_abs],3),
                    round(Chi_Enc[iCrit_abs],3),
                    round(Aeff[iCrit_abs]*10**6,3),
                    round(Red[iCrit_abs],3),
                    round(eta_Enc[iCrit_abs],3),
                    round(NTracRd[iCrit_abs] * 1e-3, 1),
                    round(NbRd[iCrit_abs] * 1e-3, 1),
                    round(np.abs(F_Compressao[iCrit_abs]) * 1e-3, 1),
                    round(Ratio_Enc[iCrit_abs], 3)
                ], dtype=object)

                Enc_Out_csv = np.vstack((Enc_Out_csv, row))
    return Enc_Out_csv

def Connection_Bracing_Exp_Output(
    TrussType, Troco, nos_vento, Truss_Out, Ratio_Lig, t_gusset, e1_gusset, e2_gusset,
    p1_gusset, p2_gusset, fy_Lig, N_Bolt_Lig, Bolt_Steel, Bolt, Fele, FLig_Rd, Ele_Lig,
    d_Bolt, d_0, As_Bolt, A_Bolt, FvRd, Ratio_bolt_lig, alpha_d_end, alpha_d_inner,
    k1_edge, k1_inner, alpha_b_end, alpha_b_inner, FbRd_Edge_Inner, FbRd_Edge_End,
    FbRd, Anv_block, Ant_block, Veff1Rd, Veff2Rd, beta_Cant, NuRd):

    # Cabeçalhos completos
    Lig_Out_csv = np.array([[
        "Member",
        "Plate thickness [mm]",
        "e1 [mm]",
        "e2 [mm]",
        "p1 [mm]",
        "p2 [mm]",
        "Steel Plate",
        "Number of bolts",
        "Bolts Size",
        "Bolts steel",
        "d_bolt [mm]",
        "d_0 [mm]",
        "As_bolt [mm²]",
        "A_bolt [mm²]",
        "FvRd [kN]",
        "Ratio bolt-lig",
        "αd_end",
        "αd_inner",
        "k1_edge",
        "k1_inner",
        "αb_end",
        "αb_inner",
        "FbRd_Edge_Inner [kN]",
        "FbRd_Edge_End [kN]",
        "FbRd [kN]",
        "Anv_block [mm²]",
        "Ant_block [mm²]",
        "Veff1Rd [kN]",
        "Veff2Rd [kN]",
        "βCant",
        "NuRd [kN]",
        "Fele [kN]",
        "FLig_Rd [kN]",
        "Ratio"
    ]], dtype=object)

    TrussType_Lig_Out = TrussType[TrussType != "Leg"]
    Troco_Lig_Out = Troco[TrussType != "Leg"]

    for j in [1, 2]:
        for i in range(1, nos_vento.shape[0]):
            # Agrupa quando j == 2 ("Horizontal Bar")
            if j == 2:
                mask_type = (
                    (TrussType_Lig_Out == "Horizontal Bar") |
                    (TrussType_Lig_Out == "External Manual Bar") |
                    (TrussType_Lig_Out == "Internal Manual Bar")
                )
            else:
                mask_type = (TrussType_Lig_Out == Truss_Out[j])

            i_Out_Zone = (Troco_Lig_Out == i) & mask_type
            indices = np.where(i_Out_Zone)[0]

            if len(indices) != 0:
                iCrit_local = np.argmax(Ratio_Lig[indices])
                iCrit_abs = indices[iCrit_local]
                Name_Enc_Out = f"{Truss_Out[j]} - Shaft {nos_vento.shape[0] - i}"

                row = np.array([
                    Name_Enc_Out,
                    np.round(t_gusset[iCrit_abs]*1e3, 1),
                    np.round(e1_gusset[iCrit_abs]*1e3, 1),
                    np.round(e2_gusset[iCrit_abs]*1e3, 1),
                    np.round(p1_gusset[iCrit_abs]*1e3, 1),
                    np.round(p2_gusset[iCrit_abs]*1e3, 1),
                    f"S{int(fy_Lig[iCrit_abs]*1e-6)}",
                    int(N_Bolt_Lig[iCrit_abs]),
                    Bolt_Steel[iCrit_abs],
                    Bolt[iCrit_abs],
                    np.round(d_Bolt[iCrit_abs]*1e3, 1),
                    np.round(d_0[iCrit_abs]*1e3, 1),
                    np.round(As_Bolt[iCrit_abs]*1e6, 1),
                    np.round(A_Bolt[iCrit_abs]*1e6, 1),
                    np.round(FvRd[iCrit_abs]*1e-3, 1),
                    np.round(Ratio_bolt_lig[iCrit_abs], 3),
                    np.round(alpha_d_end[iCrit_abs], 3),
                    np.round(alpha_d_inner[iCrit_abs], 3),
                    np.round(k1_edge[iCrit_abs], 3),
                    np.round(k1_inner[iCrit_abs], 3),
                    np.round(alpha_b_end[iCrit_abs], 3),
                    np.round(alpha_b_inner[iCrit_abs], 3),
                    np.round(FbRd_Edge_Inner[iCrit_abs]*1e-3, 1),
                    np.round(FbRd_Edge_End[iCrit_abs]*1e-3, 1),
                    np.round(FbRd[iCrit_abs]*1e-3, 1),
                    np.round(Anv_block[iCrit_abs]*1e6, 1),
                    np.round(Ant_block[iCrit_abs]*1e6, 1),
                    np.round(Veff1Rd[iCrit_abs]*1e-3, 1),
                    np.round(Veff2Rd[iCrit_abs]*1e-3, 1),
                    np.round(beta_Cant[iCrit_abs], 2),
                    np.round(NuRd[iCrit_abs]*1e-3, 1),
                    np.round(Fele[int(Ele_Lig[iCrit_abs])]*1e-3, 1),
                    np.round(FLig_Rd[iCrit_abs]*1e-3, 1),
                    np.round(Ratio_Lig[iCrit_abs], 3)
                ], dtype=object)

                Lig_Out_csv = np.vstack((Lig_Out_csv, row))

    return Lig_Out_csv

def Connection_Shaft_Exp_Output(TrussType, Troco, Ratio_Lig_B, nos_vento, divisor,
                                t_B, e1_B, e2_B, p1_B, p2_B, Corte_B,
                                fy_B, N_Bolt_B, Bolt_Steel_B, Bolt_B,
                                F_Ele_B, FLig_Rd_B,
                                d_Bolt_B, d_0_B, As_Bolt_B, A_Bolt_B, FvRd_B,
                                Ratio_bolt_B, alpha_d_end_B, alpha_d_inner_B,
                                k1_edge_B, k1_inner_B, alpha_b_end_B, alpha_b_inner_B,
                                FbRd_Edge_Inner_B, FbRd_Edge_End_B, FbRd_B,
                                Anv_block_B, Ant_block_B, Veff1Rd_B, Veff2Rd_B,
                                beta_Cant_B, NuRd_B
                            ):
    import numpy as np

    # Filtra apenas os elementos do tipo "Leg"
    TrussType_Lig_Out = TrussType[TrussType == "Leg"]
    Troco_Lig_Out = Troco[TrussType == "Leg"]

    # Cabeçalho completo (igual ao da função Bracing)
    B_Out_csv = np.array([[
        "Member",
        "Plate thickness [mm]",
        "e1 [mm]",
        "e2 [mm]",
        "p1 [mm]",
        "p2 [mm]",
        "Number of plates",
        "Steel Plate",
        "Number of bolts",
        "Bolts Size",
        "Bolts steel",
        "d_bolt [mm]",
        "d_0 [mm]",
        "As_bolt [mm²]",
        "A_bolt [mm²]",
        "FvRd [kN]",
        "Ratio bolt-lig",
        "αd_end",
        "αd_inner",
        "k1_edge",
        "k1_inner",
        "αb_end",
        "αb_inner",
        "FbRd_Edge_Inner [kN]",
        "FbRd_Edge_End [kN]",
        "FbRd [kN]",
        "Anv_block [mm²]",
        "Ant_block [mm²]",
        "Veff1Rd [kN]",
        "Veff2Rd [kN]",
        "βCant",
        "NuRd [kN]",
        "Fele [kN]",
        "FLig_Rd [kN]",
        "Ratio"
    ]], dtype=object)

    if len(Ratio_Lig_B) != 1:
        for i in range(1, nos_vento.shape[0]):
            inicio = (i - 1) * divisor
            fim = inicio + divisor
            bloco = Ratio_Lig_B[inicio:fim]

            # índice crítico
            i_local = np.argmax(bloco)
            i_global = inicio + i_local

            Name_Enc_Out = f"Leg - Shaft {nos_vento.shape[0] - i}"

            row = np.array([
                Name_Enc_Out,
                np.round(t_B[i_global]*1e3, 1),
                np.round(e1_B[i_global]*1e3, 1),
                np.round(e2_B[i_global]*1e3, 1),
                np.round(p1_B[i_global]*1e3, 1),
                np.round(p2_B[i_global]*1e3, 1),
                int(Corte_B[i_global]),
                f"S{int(fy_B[i_global]*1e-6)}",
                int(N_Bolt_B[i_global]),
                Bolt_B[i_global],
                Bolt_Steel_B[i_global],
                np.round(d_Bolt_B[i_global]*1e3, 1),
                np.round(d_0_B[i_global]*1e3, 1),
                np.round(As_Bolt_B[i_global]*1e6, 1),
                np.round(A_Bolt_B[i_global]*1e6, 1),
                np.round(FvRd_B[i_global]*1e-3, 1),
                np.round(Ratio_bolt_B[i_global], 3),
                np.round(alpha_d_end_B[i_global], 3),
                np.round(alpha_d_inner_B[i_global], 3),
                np.round(k1_edge_B[i_global], 3),
                np.round(k1_inner_B[i_global], 3),
                np.round(alpha_b_end_B[i_global], 3),
                np.round(alpha_b_inner_B[i_global], 3),
                np.round(FbRd_Edge_Inner_B[i_global]*1e-3, 1),
                np.round(FbRd_Edge_End_B[i_global]*1e-3, 1),
                np.round(FbRd_B[i_global]*1e-3, 1),
                np.round(Anv_block_B[i_global]*1e6, 1),
                np.round(Ant_block_B[i_global]*1e6, 1),
                np.round(Veff1Rd_B[i_global]*1e-3, 1),
                np.round(Veff2Rd_B[i_global]*1e-3, 1),
                np.round(beta_Cant_B[i_global], 2),
                np.round(NuRd_B[i_global]*1e-3, 1),
                np.round(F_Ele_B[i_global]*1e-3, 1),
                np.round(FLig_Rd_B[i_global]*1e-3, 1),
                np.round(Ratio_Lig_B[i_global], 3)
            ], dtype=object)

            B_Out_csv = np.vstack((B_Out_csv, row))

    return B_Out_csv

def Flanged_Shaft_Exp_Output(Ratio_F, Check_F,nos_vento, divisor,
                         N_Bolt_F, Bolt_F, Bolt_Steel_F, NTRd, F_Ele_F,
                         D_Ext_F, t_F, Steel_F,N_bolt, D_Bolt_Flange_F,
                         e2_F, e1_F, alpha_r0_F, d_Bolt_F, d_0_F,
                         As_bolt_F, A_bolt_F, Lb_F, leff, alpha_r_F,
                         epsilon_F, n_F, aws_F, m_F, Rb_F, R_F,
                         elinha_F, FtRd, k1_F, mplRd_F, Re_F, k3_F,
                         Bride, NT1Rd, NT2Rd, NT3Rd, NT4Rd):
    F_Out_csv=np.array([["Member","Flange Diameter","Flange thickness","Flange Steel",
              "Bolts Inscription","Number of bolts","Bolts Diameter","Bolts steel",
              "Ned","Nrd","Ratio","As_bolt","A_bolt","Lb","leff","alpha_r",
              "epsilon","n","aws","m","Rb","R","elinha","FtRd","k1","mplRd",
              "Re","k3","Bride","NT1Rd","NT2Rd","NT3Rd","NT4Rd",
              "N_Bolt_F","Bolt_F","Bolt_Steel_F","NTRd","F_Ele_F","Ratio_F","Check"]])
    if len(Ratio_F)!=1:
        for i in range(1, nos_vento.shape[0]):
            inicio = (i-1) * (divisor)+1
            
            fim = inicio + divisor
            bloco = Ratio_F[inicio:fim]
            i_local = np.argmax(bloco)           # posição dentro do bloco
            i_global = inicio + i_local   

            Name_Enc_Out = f"Leg - Shaft {i}"
            row = np.array([
                Name_Enc_Out,
                np.round(D_Ext_F[i_global]*10**3,0),
                np.round(t_F[i_global]*10**3,0),
                Steel_F[i_global],
                int(N_bolt[i_global]),
                np.round(D_Bolt_Flange_F[i_global]*10**3,0),
                np.round(e2_F[i_global]*10**3,1),
                np.round(e1_F[i_global]*10**3,1),
                np.round(alpha_r0_F[i_global],3),
                np.round(d_Bolt_F[i_global]*10**3,1),
                np.round(d_0_F[i_global]*10**3,1),
                np.round(As_bolt_F[i_global]*10**6,1),
                np.round(A_bolt_F[i_global]*10**6,1),
                np.round(Lb_F[i_global]*10**3,1),
                np.round(leff[i_global]*10**3,1),
                np.round(alpha_r_F[i_global],3),
                np.round(epsilon_F[i_global],3),
                np.round(n_F[i_global],3),
                np.round(aws_F[i_global]*10**3,1),
                np.round(m_F[i_global],4),
                np.round(Rb_F[i_global]*10**3,1),
                np.round(R_F[i_global]*10**3,1),
                np.round(elinha_F[i_global]*10**3,1),
                np.round(FtRd[i_global]*10**-3,1),
                np.round(k1_F[i_global],3),
                np.round(mplRd_F[i_global],3),
                np.round(Re_F[i_global]*10**3,1),
                np.round(k3_F[i_global],3),
                Bride[i_global],
                np.round(NT1Rd[i_global]*10**-3,1),
                np.round(NT2Rd[i_global]*10**-3,1),
                np.round(NT3Rd[i_global]*10**-3,1),
                np.round(NT4Rd[i_global]*10**-3,1),
                N_Bolt_F[i_global],
                Bolt_F[i_global],
                Bolt_Steel_F[i_global],
                np.round(NTRd[i_global]*10**-3,1),
                np.round(F_Ele_F[i_global]*10**-3,1),
                np.round(Ratio_F[i_global], 3),
                Check_F[i_global]
                ], dtype=object)
            F_Out_csv = np.vstack((row,F_Out_csv))
        F_Out_csv=np.flip(F_Out_csv,axis=0)
        F_Out_csv= np.vstack((F_Out_csv[0],F_Out_csv[1:][::-1]))
    
    return F_Out_csv


def Shaft_Wind(Case_Out_csv,NameCase,At,Af,Ac,Acsup,As,Indice_Cheios,cf0f,cf0c,cf0csup,cfs0,cfs,c0_Vento_For,CsCd,Pressao_troco_medio,FmW,FTW,L_lin_equiv,Ka,FmW_lin,FTW_lin,nos_vento):
    Case_Out_csv = np.vstack((Case_Out_csv, [NameCase,"","","","","","","","","","","","","","","","","","","",""]))
    Case_Out_csv = np.vstack((Case_Out_csv, ["Name_Enc_Out","At","Af","Ac","Acsup","As","Indice_Cheios","cf0f","cf0c","cf0csup","cfs0","cfs","c0_Vento_For","CsCd","Pressao_troco_medio","FmW","FTW","L_lin_equiv","Ka","FmW_lin","FTW_lin"]))
    for i in range(1,nos_vento.shape[0]):


        Name_Enc_Out = f"Shaft {nos_vento.shape[0] - i}"

        row = np.array([
            Name_Enc_Out,
            np.round(At[i],3),
            np.round(Af[i],3),
            np.round(Ac[i],3),
            np.round(Acsup[i],3),
            np.round(As[i],3),
            np.round(Indice_Cheios[i],3),
            np.round(cf0f[i],3),
            np.round(cf0c[i],3),
            np.round(cf0csup[i],3),
            np.round(cfs0[i]),
            np.round(cfs[i],3),
            np.round(c0_Vento_For[i],3),
            np.round(CsCd,3),
            np.round(Pressao_troco_medio[i],1),
            np.round(FmW[i],1),
            np.round(FTW[i],1),
            np.round(L_lin_equiv[i],3),
            np.round(Ka[i],2),
            np.round(FmW_lin[i],1),
            np.round(FTW_lin[i],1)
        ], dtype=object)
        Case_Out_csv = np.vstack((Case_Out_csv, row))
    
    return Case_Out_csv

def Output_Areas(Matrix, NameCaso, n_decimais=3):
    """
    Retorna uma matriz NumPy com cabeçalho na primeira linha e valores arredondados.

    Parâmetros:
      Matrix: lista de listas ou array NumPy com dados numéricos
      NameCaso: lista de strings com os nomes das colunas
      n_decimais: casas decimais para arredondar os valores (default=3)
    """
    # Converte para NumPy e arredonda os valores
    Matrix = np.array(Matrix, dtype=float)
    Matrix = np.round(Matrix, n_decimais)
    
    # Converte para object para permitir o cabeçalho string
    Matrix_obj = Matrix.astype(object)
    
    # Cria a primeira linha com o cabeçalho
    header = np.array([NameCaso], dtype=object)
    
    # Junta cabeçalho + matriz
    resultado = np.vstack((header, Matrix_obj))
    
    return resultado

def salvar_matrizes_csv(nome_arquivo, *matrizes):
    """
    Salva várias matrizes ou DataFrames em um único arquivo CSV.
    Mantém cabeçalhos quando o input é um DataFrame.
    """
    
    with open(nome_arquivo, 'w', encoding='utf-8') as f:
        for i, matriz in enumerate(matrizes):
            # --- Caso seja DataFrame ---
            if isinstance(matriz, pd.DataFrame):
                matriz.to_csv(f, index=False)
            
            # --- Caso seja array numpy ---
            else:
                matriz = np.array(matriz, dtype=object)
                if np.issubdtype(matriz.dtype, np.number):
                    np.savetxt(f, matriz, delimiter=',', fmt='%.6f')
                else:
                    np.savetxt(f, matriz, delimiter=',', fmt='%s')
            
            # linha em branco entre matrizes
            if i < len(matrizes) - 1:
                f.write('\n')