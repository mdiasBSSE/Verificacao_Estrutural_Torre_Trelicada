import numpy as np
def Foundation_Calc_Basic(Foundation_Matrix,FundMethod,MfundMax,RatioFundCalc,SoilSelfWeight,AllowedSoilTension_SLS,AllowedSoilTension_ULS,Mfund):
    if Foundation_Matrix is not None and Foundation_Matrix.size > 0:
        Foundation_Type=Foundation_Matrix[0,1]
        Found_Slab_Height=Foundation_Matrix[1,1]
        Found_Slab_Widht=Foundation_Matrix[2,1]
        Found_Pier_Height=Foundation_Matrix[3,1]
        Found_Pier_Widht=Foundation_Matrix[4,1]
        if FundMethod==True:
            if Foundation_Type=="Slab":
                if Found_Slab_Height==0 or Found_Slab_Widht==0:
                    FundMethod==False
            elif Foundation_Type=="Slab and Pier":
                if Found_Slab_Height==0 or Found_Slab_Widht==0 or Found_Pier_Height==0 or Found_Pier_Widht==0:
                    FundMethod==False    
            else:
                FundMethod==False 
        if FundMethod==True:
            MfundMax=MfundMax*10**3
            if Mfund>np.max(MfundMax):
                RatioFund=Mfund/MfundMax[-1]
            else:
                for i in range(len(MfundMax)-1):
                    if MfundMax[i]<Mfund and Mfund<=MfundMax[i+1]:
                        decliveFun=(RatioFundCalc[i+1]-RatioFundCalc[i])/(MfundMax[i+1]-MfundMax[i])
                        bFun=RatioFundCalc[i]-decliveFun*MfundMax[i]
                        RatioFund=decliveFun*Mfund+bFun
            #if Flag_Fund==1:
            #    CheckFund="KO"
            #    FundResults=[int(Mfund),"No data",CheckFund]
            RatioFund=np.round(RatioFund,2)
            if RatioFund>1.00:
                CheckFund="KO"
            else:
                CheckFund="OK"
            Fund_Out_Geo_Csv = np.array([["Foundation"," "],
                ["Slab Width (m)", Found_Slab_Widht],
                ["Slab Height (m)", Found_Slab_Height],
                ["Pier Width (m)", Found_Pier_Widht],
                ["Pier Height (m)", Found_Pier_Height],
                ["Soil Self-Weight (kN/m3)", SoilSelfWeight],
                ["Allowed Soil Tension (SLS) (Pa)", AllowedSoilTension_SLS],
                ["Allowed Soil Tension (ULS) (Pa)", AllowedSoilTension_ULS] 
            ], dtype=object)
            Fund_Out_Ratio_Csv=np.array([["Ratio Foundation",RatioFund]])

        else:
            Fund_Out_Geo_Csv = np.array([["Foundation"," "],
                ["Slab Width (m)", "--"],
                ["Slab Height (m)", "--"],
                ["Pier Width (m)", "--"],
                ["Pier Height (m)", "--"],
                ["Soil Self-Weight (kN/m3)", "--"],
                ["Allowed Soil Tension (SLS) (Pa)","--"],
                ["Allowed Soil Tension (ULS) (Pa)", "--"]  
            ], dtype=object)
            RatioFund=0
            Fund_Out_Ratio_Csv=np.array([["Ratio Foundation","NA"]])

    return RatioFund,Fund_Out_Geo_Csv,Fund_Out_Ratio_Csv



    