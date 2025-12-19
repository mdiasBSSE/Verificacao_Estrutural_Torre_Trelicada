import numpy as np
from Sections_Materials.Section_Properties import get_bolt_props,get_fu_fy,nut_across_flats
from Output.Output_Lattice_Reduced import Connection_Bracing_Output,Connection_Shaft_Output,Flanged_Shaft_Output
from Output.Output_Lattice_Expanded import Connection_Bracing_Exp_Output,Connection_Shaft_Exp_Output,Flanged_Shaft_Exp_Output
from Utilities.Utilities_ops import elementos_do_no_filt
from Connections.Bolted_Connection_Portugal import Bolted_Connection_Bracing_Portugal,Bolted_Connection_Shaft_Portugal,Flanged_Bolted_Connection_Portugal,Flanged_Foundation_Angle_Portugal
from Connections.Bolted_Connection_France import Bolted_Connection_Bracing_France,Bolted_Connection_Shaft_France,Flanged_Bolted_Connection_France,Flanged_Foundation_Angle_France
from Connections.Bolted_Connection_Spain import Bolted_Connection_Bracing_Spain,Bolted_Connection_Shaft_Spain,Flanged_Bolted_Connection_Spain,Flanged_Foundation_Angle_Spain
from Connections.Bolted_Connection_Italy import Bolted_Connection_Bracing_Italy,Bolted_Connection_Shaft_Italy,Flanged_Bolted_Connection_Italy,Flanged_Foundation_Angle_Italy

def Bolted_Connection_Bracing(Pais,Dig_Hor_Matrix,F_Tracao,F_Compressao,TrussType,Dim,GammaM2,ExpVento,Esp,Area,Leg_Type,Sigma_u,Sigma,GammaM0,Troco,nos_vento):    
    if Pais in ["Portugal","Portugal-RSA"]:
        Ratio_Lig,Lig_Out_csv,Lig_Out_Exp=Bolted_Connection_Bracing_Portugal(Dig_Hor_Matrix,F_Tracao,F_Compressao,TrussType,Dim,GammaM2,ExpVento,Esp,Area,Leg_Type,Sigma_u,Sigma,GammaM0,Troco,nos_vento)
        return Ratio_Lig,Lig_Out_csv,Lig_Out_Exp
    elif Pais in ["France"]:
        Ratio_Lig,Lig_Out_csv,Lig_Out_Exp=Bolted_Connection_Bracing_France(Dig_Hor_Matrix,F_Tracao,F_Compressao,TrussType,Dim,GammaM2,ExpVento,Esp,Area,Leg_Type,Sigma_u,Sigma,GammaM0,Troco,nos_vento)
        return Ratio_Lig,Lig_Out_csv,Lig_Out_Exp
    elif Pais in ["Spain"]:
        Ratio_Lig,Lig_Out_csv,Lig_Out_Exp=Bolted_Connection_Bracing_Spain(Dig_Hor_Matrix,F_Tracao,F_Compressao,TrussType,Dim,GammaM2,ExpVento,Esp,Area,Leg_Type,Sigma_u,Sigma,GammaM0,Troco,nos_vento)
        return Ratio_Lig,Lig_Out_csv,Lig_Out_Exp
    elif Pais in ["Italy"]:
        Ratio_Lig,Lig_Out_csv,Lig_Out_Exp=Bolted_Connection_Bracing_Italy(Dig_Hor_Matrix,F_Tracao,F_Compressao,TrussType,Dim,GammaM2,ExpVento,Esp,Area,Leg_Type,Sigma_u,Sigma,GammaM0,Troco,nos_vento)
        return Ratio_Lig,Lig_Out_csv,Lig_Out_Exp
def Bolted_Connection_Shaft(Pais,Shaft_Bolted_Matrix,TrussType,F_Tracao,F_Compressao,GammaM2,Sigma_u,Sigma,GammaM0,Area,Esp,Troco,nos_vento,divisor):
    if Pais in ["Portugal","Portugal-RSA"]:
        Ratio_Lig_B,B_Out_csv,B_Out_Exp=Bolted_Connection_Shaft_Portugal(Shaft_Bolted_Matrix,TrussType,F_Tracao,F_Compressao,GammaM2,Sigma_u,Sigma,GammaM0,Area,Esp,Troco,nos_vento,divisor)
        return Ratio_Lig_B,B_Out_csv,B_Out_Exp
    elif Pais in ["France"]:
        Ratio_Lig_B,B_Out_csv,B_Out_Exp=Bolted_Connection_Shaft_France(Shaft_Bolted_Matrix,TrussType,F_Tracao,F_Compressao,GammaM2,Sigma_u,Sigma,GammaM0,Area,Esp,Troco,nos_vento,divisor)
        return Ratio_Lig_B,B_Out_csv,B_Out_Exp
    elif Pais in ["Spain"]:
        Ratio_Lig_B,B_Out_csv,B_Out_Exp=Bolted_Connection_Shaft_Spain(Shaft_Bolted_Matrix,TrussType,F_Tracao,F_Compressao,GammaM2,Sigma_u,Sigma,GammaM0,Area,Esp,Troco,nos_vento,divisor)
        return Ratio_Lig_B,B_Out_csv,B_Out_Exp
    elif Pais in ["Italy"]:
        Ratio_Lig_B,B_Out_csv,B_Out_Exp=Bolted_Connection_Shaft_Italy(Shaft_Bolted_Matrix,TrussType,F_Tracao,F_Compressao,GammaM2,Sigma_u,Sigma,GammaM0,Area,Esp,Troco,nos_vento,divisor)
        return Ratio_Lig_B,B_Out_csv,B_Out_Exp
    return Ratio_Lig_B,B_Out_csv,B_Out_Exp

def Flanged_Bolted_Connection(Pais,Shaft_Flange_Matrix,TrussType,F_Tracao,Dim,Esp,Sigma_u,Sigma,GammaM2,GammaM0,Area,nos_vento,divisor):
    if Pais in ["Portugal","Portugal-RSA"]:
        Ratio_F,F_Out_csv,F_Out_Exp=Flanged_Bolted_Connection_Portugal(Shaft_Flange_Matrix,TrussType,F_Tracao,Dim,Esp,Sigma_u,Sigma,GammaM2,GammaM0,Area,nos_vento,divisor)
        return Ratio_F,F_Out_csv,F_Out_Exp
    if Pais in ["France"]:
        Ratio_F,F_Out_csv,F_Out_Exp=Flanged_Bolted_Connection_France(Shaft_Flange_Matrix,TrussType,F_Tracao,Dim,Esp,Sigma_u,Sigma,GammaM2,GammaM0,Area,nos_vento,divisor)
        return Ratio_F,F_Out_csv,F_Out_Exp
    if Pais in ["Spain"]:
        Ratio_F,F_Out_csv,F_Out_Exp=Flanged_Bolted_Connection_Spain(Shaft_Flange_Matrix,TrussType,F_Tracao,Dim,Esp,Sigma_u,Sigma,GammaM2,GammaM0,Area,nos_vento,divisor)
        return Ratio_F,F_Out_csv,F_Out_Exp
    if Pais in ["Italy"]:
        Ratio_F,F_Out_csv,F_Out_Exp=Flanged_Bolted_Connection_Italy(Shaft_Flange_Matrix,TrussType,F_Tracao,Dim,Esp,Sigma_u,Sigma,GammaM2,GammaM0,Area,nos_vento,divisor)
        return Ratio_F,F_Out_csv,F_Out_Exp

def Flanged_Foundation_Angle(Pais,Flange_Foundation_Angle,Trusstype,F_tracao):
    if Pais in ["Portugal","Portugal-RSA"]:
        ratio, row=Flanged_Foundation_Angle_Portugal(Flange_Foundation_Angle,Trusstype,F_tracao)
    if Pais in ["France"]:
        ratio, row=Flanged_Foundation_Angle_France(Flange_Foundation_Angle,Trusstype,F_tracao)
    if Pais in ["Spain"]:
        ratio, row=Flanged_Foundation_Angle_Spain(Flange_Foundation_Angle,Trusstype,F_tracao)
    if Pais in ["Italy"]:
        ratio, row=Flanged_Foundation_Angle_Italy(Flange_Foundation_Angle,Trusstype,F_tracao)
    return ratio, row