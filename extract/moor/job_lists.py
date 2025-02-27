"""
Module to create dicts for multiple (or single) mooring extractions.
"""

def get_sta_dict(job_name):
    
    # specific job definitions
    
    if job_name == 'willapa_bc': # Willapa Bay Center PCSGA Mooring
        sta_dict = {
            'wbc': (-123.9516, 46.6290)
            }
            
    elif job_name == 'mickett_1':
        sta_dict = {
        'ORCA_Hansville': (-122.6270, 47.9073),
        'ORCA_Hoodsport': (-123.1126, 47.4218),
        'ORCA_Point_Wells': (-122.3972, 47.7612),
        'Central_Main_Stem_Hood_Canal': (-122.989507, 47.574352),
        'North_Central_Main_Basin': (-122.440755, 47.825099)
        }
            
    elif job_name == 'mickett_2':
        sta_dict = {
        'Carr_Inlet_ORCA': (-122 - 43.8/60, 47 + 16.8/60),
        'East_of_Fox_Island': (-122 - 35.158/60, 47 + 13.185/60)
        }
        
    elif job_name == 'stoll_corals':
        sta_dict = {
        'Carson_D01_Lopez': (-122.8728, 48.36816),
        'Carson_D02_Admiralty': (-122.7883, 48.19252),
        'Carson_D04_Admiralty': (-122.8166, 48.19764),
        'Carson_D05_Keystone': (-122.6576, 48.12828),
        'Carson_D07_NorthAdmiralty': (-122.8898, 48.22245),
        'Carson_D08_Canada': (-123.149, 48.36136),
        'USNM_19270_Canada': (-123.233, 48.35),
        'USNM_92626_Admiralty': (-122.80, 48.1917),
        'USNM_19228_Dungeness': (-123.189, 48.225),
        'USNM_19272_Admiralty': (-122.817, 48.20),
        'USNM_92620_Lopez': (-122.85, 48.3667),
        }
            
    elif job_name == 'stoll_obs':
        sta_dict = {
        'DOE_SJF002': (-123.025, 48.25),
        'DOE_ADM002': ( -122.8417151, 48.1875056),
        'DOE_ADM001': ( -122.616715, 48.0300056),
        'WOAC_STN21': (-122.8504, 48.1883),
        'WOAC_STN20': (-122.6848, 48.142),
        'WOAC_STN19': (-122.6318, 48.0915),
        }
            
    elif job_name == 'Kelly':
        # note I pushed two of the locations a bit West to get off the landmask
        sta_dict = {
        'Seal_Rock': (-122.87004, 47.70557),
        'Little_Dewatto': (-123.08612-.005, 47.44489),
        'Red_Bluff': (-123.10438-.007, 47.41625)
        }
            
    elif job_name == 'jazzy':
        sta_dict = {
        'Middle_Bank': (-123.09651, 48.40935),
        'East_Bank': (-122.97376, 48.30042),
        'Upright_Channel': (-122.923005, 48.55410),
        'Blakely_Orcas': (-122.82880, 48.58790),
        'Rosario_Strait': (-122.74001, 48.64631),
        'North_Station': (-123.04166, 48.58330),
        'South_Station': (-122.94330, 48.42000),
        'Hein_Bank': (-123.03940, 48.35825)
        }
        
    elif job_name == 'ooi':
        sta_dict = {
            'CE01':(-124.095, 44.6598), # Oregon Inshore (25 m)
            'CE02':(-124.304, 44.6393), # Oregon Shelf (80 m)
            'CE04':(-124.956, 44.3811), # Oregon Offshore (588 m)
            'PN01A':(-125.3983, 44.5096), # Slope Base (2905 m)
        }
        
    elif job_name == 'erika_esci491w2022':
        sta_dict = {
        'Olympia': (-122.9165, 47.0823),
        'Tacoma': (-122.4758, 47.3015),
        'Seattle_West_Point': (-122.4435, 47.6813),
        'Bellingham': (-122.5519, 48.7348),
        'Central_Hood_Canal': (-122.9895, 47.5744),
        'Skokomish': (-123.1272, 47.3639),
        'Hein_Bank': (-123.0394, 48.35825),
        'Admiralty': (-122.6949, 48.1370),
        'Everett': (-122.2806, 47.9855)
        }
        
    elif job_name == 'scoot':
        sta_dict= {
        'F_090_WEH': (-126.92523957836097, 50.86003686529567),
        'F_087_SIM': (-126.88867577473951, 50.8779837900623),
        'F_083_CEC': (-126.71014703214891, 50.86003686529567),
        'M_083_087': (-126.79852585292794, 50.89607323113631),
        'F_084_CYP': (-126.65795536471262, 50.84223133403236),
        'M_083_084': (-126.69268072600828, 50.86003686529567),
        'F_088_SIR': (-126.60638135698926, 50.84223133403236),
        'M_084_089': (-126.6235045170437, 50.86003686529567),
        'F_082_BRE': (-125.34911668789157, 50.28660059365715),
        'F_089_VEN': (-125.33704071737607, 50.300007509695305),
        'F_086_RAZ': (-125.01765650844104, 50.327117635547935),
        'E_079_TOF': (-126.04793822349058, 49.21393673803513),
        'E_129_STG': (-125.12767093555838, 50.0962109697609),
        'E_004_NOO': (-126.60638135698926, 49.57825722433716),
        'E_002_ESP': (-126.9806322902372, 49.829746612851736),
        'E_016_RIV': (-126.13824776832129, 49.66598688018571),
        'F_004_CON': (-126.47181502357597, 49.65693905574285),
        'F_022_GOR': (-126.43883566700802, 49.64795343794384),
        'F_016_MUC': (-126.3414536354723, 49.63902959909838),
        'M_023_014': (-126.3414536354723, 50.607192072687155),
        'F_008_BAR': (-125.2773744875855, 50.31351066854337),
        'F_005_AHL': (-124.15671501359827, 49.77957267482829),
        'F_013_VAN': (-123.85335983884296, 49.675097341923546),
        'F_011_SAL': (-123.83316138272168, 49.62136556218934),
        'F_006_NEW': (-123.65810809633727, 49.64795343794384),
        'E_006_RIV': (-123.5436501783167, 49.693507914816124)}
        
    elif job_name == 'kastner':
        sta_dict = {
        'Canal_Mouth': (-122.637493, 47.928439),
        'Bridge': (-122.621784, 47.858911),
        'Joint': (-122.819507, 47.669810),
        'Dabob_Bay_Entrance': (-122.860989, 47.693196),
        'Dabob_Bay_Head': (-122.805422, 47.808231),
        'Duckabush_River': (-122.908291, 47.633095),
        'Hama_Hama': (-123.026320, 47.534954),
        'Lilliwaup': (-123.088399, 47.456583),
        'ORCA_Hoodsport': (-123.1126, 47.4218),
        'Skokomish': (-123.127835, 47.363217),
        'Sisters_Point': (-123.022404, 47.358448),
        'ORCA_Twanoh': (-123.0083, 47.375),
        'Head': (-122.893559, 47.411036)
        }
        
    elif job_name == 'ROMS_update':
        sta_dict = {
        'ORCA_Hoodsport': (-123.1126, 47.4218),
        'CE02':(-124.304, 44.6393) # Oregon Shelf (80 m)
        }

    elif job_name == 'alpe':
        sta_dict = {
        'superplot': (0,45.4),
        'Estuary_Mouth': (0,45.1)
        }

    elif job_name == 'alpe2':
        sta_dict = {
        'superplot': (0,45.3),
        'Estuary_Mouth': (0,45.1)
        }

    # 07/07/2022 added job to extract moorings for wind experiment analyses
    elif job_name == 'wind_ana':
        sta_dict = {
        'eastern_moor': (1.0,44.85),
        'central_moor': (0.0,44.85),
        'western_moor': (-1.0,44.85)
        }

    # 2022.11.09 extract mooring from Birch Bay treatment plant in LiveOcean
    elif job_name == 'birchbay_wwtp':
        sta_dict = {
        'wwtp': (-122.8030401205365,48.8975820886066),
        'shifted-wwtp': (-122.80977293924359,48.8975820886066)
        }

    # 2023.01.13 extract mooring from Oak Harbor Lagoon treatment plant in LiveOcean
    elif job_name == 'oakharborlagoon_wwtp':
        sta_dict = {
        'wwtp': (-122.60105555932371,48.28559664274097),
        }

    # 2023.03.14 extract mooring wwtp4 in idealized estuary
    elif job_name == 'wwtp4':
        sta_dict = {
        'wwtp': (0.6001938503633347,44.34291446936327),
        }

    # 2023.05.23 preliminary bottom DO in LiveOcean
    elif job_name == 'pennlynch':
        sta_dict = {
        'PennCove':  (-122.714423,48.226958),
        'LynchCove': (-123.0083,47.3750) # orca buoy location
        }
    
    # 2023.05.25 preliminary bottom DO in LiveOcean
    elif job_name == 'dabobportorchard':
        sta_dict = {
        'DabobBay':  (-122.808721,47.785261),
        'PortOrchard': (-122.577043,47.690751)
        }

    # 2023.07.17 orca buoy locations (from Erin's job lists)
    elif job_name == 'orca':
        sta_dict = {
        'CI': (-122.7300, 47.2800),
        'PW': (-122.3972, 47.7612),
        'NB': (-122.6270, 47.9073),
        'DB': (-122.8029, 47.8034),
        'HP': (-123.1126, 47.4218),
        'TW': (-123.0083, 47.3750)
        }

    # 2024.02.05 exploring influence of WWTP nutrients in Puget Sound
    elif job_name == 'noWWTPNtest':
        sta_dict = {
        'main_basin': (-122.445051, 47.615672),
        'holmes_harbor': (-122.528470, 48.055797),
        'hood_canal': (-122.813348, 47.667696), 
        'admiralty_sill': (-122.700874, 48.138784)
        }

    # random point near river in hcal model
    elif job_name == 'hcal':
        sta_dict = {
        'near_river': (-122.6,47.2),
        'middle': (-122.6,47.5),
        }

    # 2024.07.03 Extraction locations for 21 inlets in Puget Sound, for the year 2014
    elif job_name == 'twentyoneinlets':
        sta_dict = {
        'similk': (-122.566493, 48.433323),       # 01. Similk Bay
        'oak': (-122.643288, 48.282798),          # 02. Oak Harbor
        'crescent': (-122.594159, 48.284254),     # 03. Crescent Harbor
        'penn': (-122.688449, 48.230183),         # 04. Penn Cove
        'portsusan': (-122.407245, 48.142236),    # 05. Port Susan
        'killsut': (-122.710107, 48.058841),      # 06. Killsut Harbor
        'holmes': (-122.530311, 48.052192),       # 07. Holmes Harbor
        'dabob': (-122.821167, 47.763285),        # 08. Dabob Bay
        'dyes': (-122.686548, 47.616928),         # 09. Dyes Inlet
        'sinclair': (-122.641701, 47.549999),     # 10. Sinclair Inlet (Ecology monitoring station SIN001)
        'elliot': (-122.368301, 47.596668),       # 11. Elliot Bay (Ecology monitoring station ELB015) 
        'lynchcove': (-122.928299, 47.398331),    # 12. Lynch Cove (Ecology monitoring station HCB007)  
        'lynchcove2': (-123.023300, 47.356670),   #                (Ecology monitoring station HCB004) - updated 2024.11.16
        'case': (-122.807292, 47.332090),         # 13. Case Inlet
        'carr': (-122.708297, 47.276668),         # 14. Carr Inlet (Ecology monitoring station CRR001) - updated 2024.11.16
        'quartermaster': (-122.470469, 47.377077),# 15. Quartermaster Harbor
        'commencement': (-122.448303, 47.290001), # 16. Commencement Bay (Ecology monitoring station CMB003)
        'hammersley': (-123.076698, 47.213329),   # 17. Hammersley Inlet (Ecology monitoring station OAK004) 
        'totten': (-123.019997, 47.121670),       # 18. Totten Inlet (Ecology monitoring station TOT002)
        'eld': (-122.949787, 47.115069),          # 19. Eld Inlet
        'budd': (-122.916702, 47.091671),         # 20. Budd Inlet (Ecology monitoring station BUD005)
        'henderson': (-122.834437, 47.143149),    # 21. Henderson Inlet
        }

    # 2024.11.16 Extraction locations Ecology CTD stations in Main Basin
    elif job_name == 'mainbasin_EcolCTD':
        sta_dict = {
        'ADM002': (-122.8416976928711,48.1875),
        'PTH005': (-122.76329803466797,48.08332824707031),
        'ADM001': (-122.61669921875,48.029998779296875),
        'ADM003': (-122.48179626464844,47.87916946411133),
        'PSB003': (-122.44170379638672,47.65999984741211),
        'EAP001': (-122.37999725341797,47.41667175292969)
        }
            
            
    else:
        print('Unsupported job name!')
        a = dict()
        return a
        
    return sta_dict