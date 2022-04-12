# Metku Winrami adapter

## TODO
Write C++ class that implements Winrami interface*
  * CWinramiControl:
    - Model building:
        - AddGP_CRD
        - AddGP_GP
        - AddGP_Vali
        - AddVali
        - UnSelectALL
        - SelectOsa
        - SelectPoint
        - SelectPoints
        - AddKTapaus
        - AddYhdistely
        - AddTasainenKuorma
        - GenValiSelection
        - SetReunaEhto
        - SetStabilityOne
        - SetFrameStatus
        - AddEmbProfile
        - ASYNC_AddEmbProfile
        - SetYhdistelyParams
        - SetNivelöi
        - SetLuotettavuusKerroin2015
        - AddFastLiitos
        - AddNPitDefConst
        - GenValiSelection_AddPisteKuorma

    - Basics:
        - ASYNC_Calculate
        - SetNow
        - GetOsaPaino
        - GetOminaisReaktio
        - GetMomForValiNew
        - GetMaxU
        - GetMaxW
        - GetNodePosMom2015
        - GetPosN
        - GetPosV

    - Design of members and joints:
        - SetVapaaVali2014
        - SetVAhvistukset2014
        - SetAutomationScanPituus (or Basic? Check!)
        - ASYNC_GetMItoitusStrings
        - SetAutomationForceBeg
        - GetLiitosStringsNew

    - Fire design:
        - SetPaloluokka
        - GetOsaKriittinen

    - Cost calculation:
        - GetNonProdCost

*small changes in API are accepted
