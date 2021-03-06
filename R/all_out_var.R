#' Return all possible STICS outputs for var.mod
#'
#' @description Helper function to print the list of all possible variables to set as output
#' from the STICS model.
#'
#' @seealso [set_out_var()]
#'
#' @examples
#' library(sticRs)
#' All_vars= all_out_var()
#'
#' @export
#'
all_out_var= function(){
  c("abso(n)",
    "age_prairie",
    "airg(n)",
    "albedolai",
    "allocfruit",
    "ammomes",
    "amptcultmat",
    "anit(n)",
    "anit_engrais(n)",
    "anit_uree(n)",
    "anoxmoy",
    "AZamm(1)",
    "AZamm(2)",
    "AZamm(3)",
    "AZamm(4)",
    "AZamm(5)",
    "azlesd",
    "AZnit(1)",
    "AZnit(2)",
    "AZnit(3)",
    "AZnit(4)",
    "AZnit(5)",
    "azomes",
    "bouchon",
    "Cb",
    "Cbmulch",
    "cdemande",
    "cEdirect",
    "cEdirecttout",
    "cep",
    "ces",
    "cestout",
    "cet",
    "cet_from_lev",
    "cetm",
    "Cetmtout",
    "cetp",
    "chargefruit",
    "Chuma",
    "Chumi",
    "Chumt",
    "cintermulch",
    "cinterpluie",
    "Cmulch",
    "Cmulchdec",
    "Cmulchnd",
    "CNgrain",
    "Cnondec(1)",
    "Cnondec(10)",
    "Cnondec(2)",
    "Cnondec(3)",
    "Cnondec(4)",
    "Cnondec(5)",
    "Cnondec(6)",
    "Cnondec(7)",
    "Cnondec(8)",
    "Cnondec(9)",
    "CNplante",
    "co2(n)",
    "CO2hum",
    "CO2res",
    "CO2sol",
    "codebbch_output",
    "concNO3les",
    "concNO3sol(1)",
    "concNO3sol(2)",
    "concNO3sol(3)",
    "concNO3sol(4)",
    "concNO3sol(5)",
    "condenit",
    "couvermulch",
    "cpluie",
    "cprecip",
    "cpreciptout",
    "Cr",
    "Crac",
    "Cresiduprofil(1)",
    "Cresiduprofil(10)",
    "Cresiduprofil(2)",
    "Cresiduprofil(3)",
    "Cresiduprofil(4)",
    "Cresiduprofil(5)",
    "Cresiduprofil(6)",
    "Cresiduprofil(7)",
    "Cresiduprofil(8)",
    "Cresiduprofil(9)",
    "crg",
    "crgtout",
    "CsurNres_pature",
    "ctairtout",
    "ctcult",
    "ctculttout",
    "ctetptout",
    "ctmoy",
    "Ctousresidusprofil",
    "cum_et0",
    "cum_et0_from_lev",
    "cum_immob",
    "cumlracz",
    "cumraint",
    "cumrg",
    "cumvminh",
    "cumvminr",
    "da(1)",
    "da(2)",
    "day_after_sowing",
    "day_cut",
    "deltai(n)",
    "deltaz",
    "demande",
    "densite",
    "densiteequiv",
    "dfol",
    "diftemp1intercoupe",
    "diftemp2intercoupe",
    "dltags",
    "dltaisen",
    "dltams(n)",
    "dltamsen",
    "dltaremobil",
    "dltmsrac_plante",
    "drain",
    "drain_from_plt",
    "drain_from_lev",
    "drat",
    "drlsenmortalle",
    "dtj(n)",
    "dureehumec",
    "dureeRH",
    "durvie(n)",
    "ebmax",
    "ebmax_gr",
    "eai",
    "Edirect",
    "efda",
    "efdensite",
    "efdensite_rac",
    "efNrac_mean",
    "em_N2O",
    "em_N2Oden",
    "em_N2Onit",
    "emd",
    "emulch",
    "eo",
    "eop",
    "eos",
    "ep",
    "epc_recal(1)",
    "epc_recal(2)",
    "epc_recal(3)",
    "epc_recal(4)",
    "epc_recal(5)",
    "epsib",
    "esol",
    "et",
    "et0",
    "etm",
    "etpp(n)",
    "exces(1)",
    "exces(2)",
    "exces(3)",
    "exces(4)",
    "exces(5)",
    "exobiom",
    "exofac",
    "exofac1moy",
    "exofac2moy",
    "exolai",
    "fapar",
    "fco2",
    "fco2s",
    "fgelflo",
    "fixmaxvar",
    "fixpot",
    "fixreel",
    "flurac",
    "flusol",
    "fpari",
    "fpari_gr",
    "fpft",
    "fpv(n)",
    "FsNH3",
    "fstressgel",
    "ftemp",
    "fxa",
    "fxn",
    "fxt",
    "fxw",
    "gel1",
    "gel1_percent",
    "gel2",
    "gel2_percent",
    "gel3",
    "gel3_percent",
    "H2Orec",
    "H2Orec_percent",
    "hauteur",
    "dominant",
    "varrapforme",
    "Hmax",
    "Hnappe",
    "Hpb",
    "Hph",
    "HR(1)",
    "HR(2)",
    "HR(3)",
    "HR(4)",
    "HR(5)",
    "HR_vol_1_10",
    "HR_vol_1_30",
    "HR_vol_121_150",
    "HR_vol_151_180",
    "HR_vol_31_60",
    "HR_vol_61_90",
    "HR_vol_91_120",
    "huile",
    "huile_percent",
    "humair",
    "humair_percent",
    "humidite",
    "humidite_percent",
    "humirac_mean",
    "iamfs",
    "idebdess",
    "idebdorms",
    "idrps",
    "ifindorms",
    "iflos",
    "igers",
    "ilans",
    "ilaxs",
    "ilevs",
    "imats",
    "imontaisons",
    "infil_recal(1)",
    "infil_recal(2)",
    "infil_recal(3)",
    "infil_recal(4)",
    "infil_recal(5)",
    "inn",
    "inn1intercoupe",
    "inn1moy",
    "inn2intercoupe",
    "inn2moy",
    "innlai",
    "inns",
    "innsenes",
    "inous",
    "intermulch",
    "interpluie",
    "iplts",
    "irazo(n)",
    "ircarb(n)",
    "irecs",
    "irrigjN",
    "irrigN",
    "isens",
    "izrac",
    "lai(n)",
    "lai(ao)",
    "lai(as)",
    "lai_mx_av_cut",
    "laimax",
    "laisen(n)",
    "largeur",
    "leaching_from_plt",
    "leaching_from_lev",
    "leai",
    "lessiv",
    "LRACH(1)",
    "LRACH(2)",
    "LRACH(3)",
    "LRACH(4)",
    "LRACH(5)",
    "lracsentot",
    "mabois",
    "maenfruit",
    "mafeuil",
    "mafeuil_kg_ha",
    "mafeuiljaune",
    "mafeuiltombe",
    "mafeuilverte",
    "mafrais",
    "mafruit",
    "mafruit_kg_ha",
    "masec(n)",
    "masec_kg_ha",
    "masec_mx_av_cut",
    "masecneo",
    "masectot",
    "masecveg",
    "matigestruc",
    "matigestruc_kg_ha",
    "matuber",
    "mortalle",
    "mortmasec",
    "mortreserve",
    "MSexporte",
    "msjaune",
    "msneojaune",
    "msrac(n)",
    "msrec_fou",
    "MSrecycle",
    "msresjaune",
    "N_mineralisation",
    "N_volatilisation",
    "Nb",
    "nbfeuille",
    "nbinflo_recal",
    "nbj0remp",
    "nbjechaudage",
    "nbjgel",
    "nbjpourdecirecolte",
    "nbjpourdecisemis",
    "Nbmulch",
    "NCbio",
    "Ndenit",
    "Nexporte",
    "nfruit(1)",
    "nfruit(2)",
    "nfruit(3)",
    "nfruit(4)",
    "nfruit(5)",
    "nfruit(nboite)",
    "nfruit(nboite-1)",
    "nfruitnou",
    "Nhuma",
    "Nhumi",
    "Nhumt",
    "nitetcult(n)",
    "nitrifj",
    "Nmineral_from_plt",
    "Nmineral_from_lev",
    "Nmulchdec",
    "Nmulchnd",
    "Nnondec(1)",
    "Nnondec(10)",
    "Nnondec(2)",
    "Nnondec(3)",
    "Nnondec(4)",
    "Nnondec(5)",
    "Nnondec(6)",
    "Nnondec(7)",
    "Nnondec(8)",
    "Nnondec(9)",
    "nodn",
    "Norgeng",
    "Nr",
    "Nrac",
    "Nrecycle",
    "Nresiduprofil(1)",
    "Nresiduprofil(10)",
    "Nresiduprofil(2)",
    "Nresiduprofil(3)",
    "Nresiduprofil(4)",
    "Nresiduprofil(5)",
    "Nresiduprofil(6)",
    "Nresiduprofil(7)",
    "Nresiduprofil(8)",
    "Nresiduprofil(9)",
    "Ntousresidusprofil",
    "numcoupe",
    "numcult",
    "Nvolat_from_plt",
    "Nvolat_from_lev",
    "Nvoleng",
    "offrenod",
    "p1000grain",
    "pdsfruit(1)",
    "pdsfruit(2)",
    "pdsfruit(3)",
    "pdsfruit(4)",
    "pdsfruit(5)",
    "pdsfruit(nboite)",
    "pdsfruit(nboite-1)",
    "pdsfruitfrais",
    "penfruit",
    "pfeuil(n)",
    "pfeuiljaune",
    "pfeuilverte(n)",
    "phoi",
    "pHvol",
    "pousfruit",
    "poussracmoy",
    "precip",
    "precipjN",
    "precipN",
    "preserve",
    "profexteau",
    "profextN",
    "profnappe",
    "psibase",
    "ptigestruc",
    "QCapp",
    "QCO2hum",
    "QCO2mul",
    "QCO2res",
    "QCO2sol",
    "QCplantetombe",
    "QCprimed",
    "QCrac",
    "QCresorg",
    "QCressuite",
    "QCrogne",
    "Qdrain",
    "Qdraincum",
    "Qem_N2O",
    "Qem_N2Oden",
    "Qem_N2Onit",
    "Qfix",
    "Qles",
    "Qlesd",
    "Qminh",
    "Qminr",
    "qmulch",
    "QNapp",
    "QNdenit",
    "QNdenit_from_plt",
    "QNdenit_from_lev",
    "QNexport",
    "QNgrain",
    "Qnitrif",
    "QNorgeng",
    "QNplante",
    "QNplante_mx_av_cut",
    "QNplantetombe",
    "QNprimed",
    "QNrac",
    "QNresorg",
    "QNressuite",
    "QNrogne",
    "QNvoleng",
    "QNvolorg",
    "qres_pature",
    "Qressuite",
    "Qressuite_tot",
    "ra_recal",
    "raint",
    "ras",
    "Ratm",
    "rc",
    "rdif",
    "remobilj",
    "remontee",
    "rendementsec",
    "resmes",
    "resperenne",
    "resrac",
    "rfpi",
    "rfvi",
    "rlj",
    "rltot",
    "rmaxi",
    "rnet",
    "rnetS",
    "rnet_plant",
    "rombre",
    "rsoleil",
    "RsurRU",
    "RsurRUrac",
    "RU",
    "ruissel",
    "ruisselsurf",
    "ruisselt",
    "runoff_from_plt",
    "runoff_from_lev",
    "RUrac",
    "saturation",
    "senfac",
    "sla",
    "SoilAvW",
    "SoilN",
    "SoilNM",
    "SoilWatM",
    "som_HUR",
    "som_sat",
    "somcour",
    "somcourdrp",
    "somcourfauche",
    "somcourmont",
    "somdifftculttair",
    "somtemp",
    "somudevair",
    "somudevcult",
    "somupvtsem",
    "sourcepuits",
    "spfruit",
    "splai",
    "stemflow",
    "str1intercoupe",
    "str2intercoupe",
    "stu1intercoupe",
    "stu2intercoupe",
    "sucre",
    "sucre_percent",
    "surf(ao)",
    "surf(as)",
    "swfac",
    "swfac1moy",
    "swfac2moy",
    "tairveille",
    "tauxcouv(n)",
    "tcult",
    "tcult_tairveille",
    "tcultmax",
    "tcultmin",
    "tempeff",
    "tetp(n)",
    "tetstomate",
    "teturg",
    "tmax(n)",
    "tmaxext(n)",
    "tmin(n)",
    "tminext(n)",
    "tmoy(n)",
    "tmoyext(n)",
    "tmoyIpltJuin",
    "tmoyIpltSept",
    "tncultmat",
    "tnhc",
    "tnrc",
    "totapN",
    "totapNres",
    "totir",
    "tpm(n)",
    "trg(n)",
    "trgext(n)",
    "trr(n)",
    "TS(1)",
    "TS(2)",
    "TS(3)",
    "TS(4)",
    "TS(5)",
    "turfac",
    "turfac1moy",
    "turfac2moy",
    "tustress",
    "tvent(n)",
    "udevair",
    "udevcult",
    "ulai(n)",
    "upvt(n)",
    "vitmoy",
    "xmlch1",
    "zrac",
    "tsol(10)",
    "hur_10_vol")
}


#' Find output variable names for STICS
#'
#' @description Helper function that return the name used as input
#' for the STICS model with a partial match.
#'
#' @param Var Character vector with a (partial) STICS output variable name
#'
#' @details The function understand [base::regex()] as input.
#'
#' @seealso [all_out_var()]
#'
#' @examples
#' library(sticRs)
#' find_STICS_var("lai")
#'
#' @export
#'
find_STICS_var= function(Var){
  All_vars= all_out_var()
  All_vars[grep(Var, All_vars,ignore.case = TRUE)]
}
