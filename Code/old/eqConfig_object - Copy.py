import numpy as np
import collections

    
def setup(self,config,inside=False):
    conformal = True
    dSOL=0.005  # [m]
    Nsol=51
    dRfw = 0.25  # [m]
    div_ex = 0.3
    TFopp = 'V'  # TF shape optimisation 'L'==length, 'V'==volume
    sheild_base = -1
    tfw = 0.1  # first wall thickness
    tBBsupport = 0.1  # blanket support
    #BS = np.array([1.075,1.56])  # blanket+sheild (in/out)
    BS = np.array([0.78-0.2,1.304-0.2])  # blanket+sheild (in/out) #PMI xls
    BS -= (tfw+tBBsupport)
    BBfrac = np.array([0.8,0.65])
    BBfrac = np.array([1,1])
    BB = list(BBfrac*BS)  # breeding blanket thickness (in/out)
    sheild = list((1-BBfrac)*BS)  # neutron sheilding
    #sheild = [0.2,0.6]
    #VV = 0.32  # vacumn vessel thickness (in/out)
    VV = [0.587,1.2]  # vacumn vessel thickness (in/out)
    graze = 1.5*np.pi/180  # toroidal grazing angle
    dPlate = 0.2  # target plate length
    nTF = 18  # number of TF coils
    Jmax = 30e6  # 15e6  # max current density
    trim = [0.93,0.81]  # trim fraction (in/out)
    config = config
    targets = collections.OrderedDict()
   
    if config == 'vde':
        dRfw = 0.25  # [m]
        TFopp = 'L'
        inflate = 0
        div_ex = 1
        dPlate = 0.5  # target plate length
        sheild_connect=[0.15,0.9]
        sheild_base = -1
        Dsheild =[0,1]  # inner outer 0-1
        Dsheild =[]  # inner outer 0-1
        trim = [0.71,0.74]  # trim fraction (in/out)
        filename = '../eqdsk/vde.eqdsk'
        dataname = 'SX7'
        coils = {'internal':{'id':[],'dR':[],'dZ':[]},
                 'external':{'id':list(range(10,16)),'dR':[],'dZ':[]}}
        targets['inner'] = {'L2D':[1.1],'open':False,'graze':graze,'dR':0.1}
        targets['outer'] = {'L2D':[3.55],'open':False,'graze':graze,'dR':0.1}
    
    if config == 'SX8_IDM':
        dRfw = 0.25  # [m]
        TFopp = 'L'
        inflate = 0
        div_ex = 1
        dPlate = 0.5  # target plate length
        sheild_connect=[0.15,0.9]
        sheild_base = -1
        Dsheild =[0,1]  # inner outer 0-1
        Dsheild =[]  # inner outer 0-1
        trim = [0.71,0.74]  # trim fraction (in/out)

        filename = '../eqdsk/2015_SX_ext_coils_eqdk_2MK6XX_v1_0.eqdsk'  # IDM

        dataname = 'SX8_IDM'
        coils = {'internal':{'id':[],'dR':[],'dZ':[]},
                 'external':{'id':list(range(10,16)),'dR':[],'dZ':[]}}
        targets['inner'] = {'L2D':[1.1],'open':False,
                                 'graze':graze,'dR':0}
        targets['outer'] = {'L2D':[3.55],'open':False,
                                 'graze':graze,'dR':0} 
                                 
    if config == 'SX8':
        dRfw = 0.25  # [m]
        TFopp = 'L'
        inflate = 0
        div_ex = 1
        dPlate = 0.5  # target plate length
        sheild_connect=[0.15,0.9]
        sheild_base = -1
        Dsheild =[0,1]  # inner outer 0-1
        Dsheild =[]  # inner outer 0-1
        trim = [0.71,0.74]  # trim fraction (in/out)
        filename = '../eqdsk/Equil_AR3d1_16coils_bt_1d03li_0d8_I'
        filename += 'pl_20d25_SX_on_SFgeom_NOEDDY_EOF_fine_iter5_v3.eqdsk'
        
        #filename = '../eqdsk/Equil_AR3d1_16coils_bt_1d03li_0d8_I'
        #filename += 'pl_20d25_SX_on_SFgeom_NOEDDY_EOF_fine_iter5_v3_new.eqdsk'
        filename = '../eqdsk/2015_SX_ext_coils_eqdk_2MK6XX_v1_0.eqdsk'  # IDM
        #filename = '../eqdsk/Equil_AR3d1_16coils_bt_1d03li_0d8_Ipl_'
        #filename += '20d25_SX_on_SFgeom_NOEDDY_EOF_fine_iter3_v5_old.eqdsk'

        dataname = 'SX8'
        coils = {'internal':{'id':[],'dR':[],'dZ':[]},
                 'external':{'id':list(range(10,16)),'dR':[],'dZ':[]}}
        targets['inner'] = {'L2D':[1.1],'open':False,
                                 'graze':graze,'dR':0}
        targets['outer'] = {'L2D':[3.55],'open':False,
                                 'graze':graze,'dR':0}   

    if config == 'SX7':
        dRfw = 0.25  # [m]
        TFopp = 'L'
        inflate = 0
        div_ex = 1
        dPlate = 0.5  # target plate length
        sheild_connect=[0.15,0.9]
        sheild_base = -1
        Dsheild =[0,1]  # inner outer 0-1
        Dsheild =[]  # inner outer 0-1
        trim = [0.71,0.74]  # trim fraction (in/out)
        filename = '../eqdsk/Equil_AR3d1_16coils_bt_1d03li_0d8_I'
        filename += 'pl_20d25_SX_on_SFgeom_NOEDDY_EOF_fine_iter5_v2.eqdsk'
        dataname = 'SX7'
        coils = {'internal':{'id':[],'dR':[],'dZ':[]},
                 'external':{'id':list(range(10,16)),'dR':[],'dZ':[]}}
        targets['inner'] = {'L2D':[1.1],'open':False,'graze':graze,'dR':0.1}
        targets['outer'] = {'L2D':[3.55],'open':False,'graze':graze,'dR':0.1}

    if config == 'SX6':
        dRfw = 0.25  # [m]
        TFopp = 'L'
        inflate = 0
        div_ex = 1
        dPlate = 0.5  # target plate length
        sheild_connect=[0.15,0.9]
        sheild_base = -1
        Dsheild =[0,1]  # inner outer 0-1
        Dsheild =[]  # inner outer 0-1
        trim = [0.71,0.74]  # trim fraction (in/out)
        filename = '../eqdsk/Equil_AR3d1_16coils_bt_1d03li_0d8_I'
        filename += 'pl_20d25_SX_on_SFgeom_NOEDDY_EOF_fine_iter4_v3.eqdsk'
        dataname = 'SX6'
        coils = {'internal':{'id':[],'dR':[],'dZ':[]},
                 'external':{'id':list(range(10,16)),'dR':[],'dZ':[]}}
        targets['inner'] = {'L2D':[1.1],'open':False,'graze':graze,'dR':0.1}
        targets['outer'] = {'L2D':[3.55],'open':False,'graze':graze,'dR':0.1}


    if config == 'SX5':
        dRfw = 0.25  # [m]
        TFopp = 'L'
        inflate = 0
        div_ex = 0.3
        dPlate = 0.5  # target plate length
        sheild_connect=[0.2,1]
        sheild_base = 0.1
        Dsheild =[0,1]  # inner outer 0-1
        Dsheild =[]  # inner outer 0-1
        trim = [0.71,0.74]  # trim fraction (in/out)
        filename = '../eqdsk/Equil_AR3d1_16coils_bt_1d03li_0d8_I'
        filename += 'pl_20d25_SX_on_SFgeom_NOEDDY_EOF_fine_iter3_v5.eqdsk'
        dataname = 'SX5'
        coils = {'internal':{'id':[],'dR':[],'dZ':[]},
                 'external':{'id':[],'dR':[],'dZ':[]}}
        targets['inner'] = {'L2D':[1.1],'open':False,'graze':graze,'dR':0.105}
        targets['outer'] = {'L2D':[4.15],'open':False,'graze':graze,'dR':0.105}

    elif config == 'SX4':
        dRfw = 0.4  # [m]
        TFopp = 'L'
        inflate = 0
        div_ex = 0.1
        dPlate = 0.25  # target plate length
        sheild_connect=[0.25,1]
        sheild_base = 0.1
        Dsheild =[0,1]  # inner outer 0-1
        Dsheild =[]  # inner outer 0-1
        trim = [0.715,0.72]  # trim fraction (in/out)
        filename = '../eqdsk/Equil_AR3d1_16coils_bt_1d03li_0d8_I'
        filename += 'pl_20d25_SX_on_SFgeom_NOEDDY_EOF_fine_iter2.eqdsk'
        dataname = 'SX4'
        coils = {'internal':{'id':[],'dR':[],'dZ':[]},
                 'external':{'id':[],'dR':[],'dZ':[]}}
        targets['inner'] = {'L2D':[1.1],'open':False,'graze':graze,'dR':0.1}
        targets['outer'] = {'L2D':[4.15],'open':False,'graze':graze,'dR':0.8}

    elif config == 'SX3':
        TFopp = 'L'
        inflate = 0
        div_ex = 0.35
        dPlate = 0.2  # target plate length
        sheild_connect=[0.3,1]
        sheild_base = 0.1
        Dsheild =[0,1]  # inner outer 0-1
        Dsheild =[]  # inner outer 0-1
        trim = [0.76,0.72]  # trim fraction (in/out)
        filename = '../eqdsk/Equil_AR3d1_16coils_bt_1d03li_0d8_I'
        filename += 'pl_20d25_SX_on_SFgeom_NOEDDY_EOF_fine_v2.eqdsk'

        dataname = 'SX3'
        coils = {'internal':{'id':[],'dR':[],'dZ':[]},
                 'external':{'id':[],'dR':[],'dZ':[]}}

        targets['inner'] = {'L2D':[1.5],'open':False,'graze':graze,'dR':0.1}
        targets['outer'] = {'L2D':[3],'open':False,'graze':graze,'dR':0.3}

    elif config == 'SX2':
        TFopp = 'L'
        inflate = 0.125
        div_ex = 0.2
        dPlate = 0.4  # target plate length
        sheild_connect=[0.25,0.9]
        sheild_base = 0.1
        Dsheild =[0.045,0.68]  # inner outer 0-1
        Dsheild =[]  # inner outer 0-1
        trim = [0.58,0.48]  # trim fraction (in/out)
        filename = '../eqdsk/redtosx4_1.eqdsk'
        dataname = 'SX2'
        coils = {'internal':{'id':[],'dR':[],'dZ':[]},
                 'external':{'id':[],'dR':[],'dZ':[]}}
        if inside:
            dataname += 'inside'
            targets['inner'] = {'L2D':[3],'open':False,'graze':graze}
            targets['outer'] = {'L2D':[5.5],'open':True,'graze':graze}
        else:
            targets['inner'] = {'L2D':[0.7],'open':False,'graze':graze,'dR':0.1}
            targets['outer'] = {'L2D':[6.45],'open':False,'graze':graze,'dR':0.3}

    elif config == 'SX':
        inflate = 0.2
        TFopp = 'V'
        div_ex = 2
        dPlate = 0.5  # target plate length
        sheild_connect=[0.25,0.9]
        Dsheild =[0.04,0.58]  # inner outer 0-1
        Dsheild =[]  # inner outer 0-1
        trim = [0.62,0.59]  # trim fraction (in/out)
        filename = '../eqdsk/Super-X_equilibrium_2MC43F_v1_0.eqdsk'
        dataname = 'SX'
        #coils = {'internal':{'id':[11,12,12],'dR':[0.6],'dZ':[]},
        #         'external':{'id':list(range(0,11))+[13],'dR':[],'dZ':[]}}
        coils = {'internal':{'id':[],'dR':[],'dZ':[]},
                      'external':{'id':[],'dR':[],'dZ':[]}}
        targets['inner'] = {'L2D':[1.4],'open':False,
                                 'graze':graze,'dR':0}
        targets['outer'] = {'L2D':[5.2],'open':False,
                                 'graze':graze,'dR':0}
            
    elif config == 'SXex':
        inflate = 0
        TFopp = 'L'
        div_ex = 1.3
        #dPlate = 5.5  # target plate length
        sheild_connect=[0.2,0.9]
        Dsheild =[]  # inner outer 0-1
        sheild_base = 0  # 0.25
        trim = [0.63,0.68]  # trim fraction (in/out)
        #filename = '../eqdsk/Super-X_equilibrium_2MC43F_v1_0.eqdsk'
        filename = './plot_data/SXex.eqdsk'
        dataname = 'SXex'
        #coils = {'internal':{'id':[11,12,12],'dR':[0.6],'dZ':[]},
        #         'external':{'id':list(range(0,11))+[13],'dR':[],'dZ':[]}}
        coils = {'internal':{'id':[],'dR':[],'dZ':[]},
                      'external':{'id':[],'dR':[],'dZ':[]}}
        targets['inner'] = {'L2D':[1.0],'open':False,
                                 'graze':graze,'dR':0.15,
                                 'dPlate':0.25}
        targets['outer'] = {'L2D':[5.615],'open':False,  # 5.607
                                 'graze':graze,'dR':0,
                                 'dPlate':1.0}  # 7.5

    elif config == 'X':
        dPlate = 0.2  # target plate length
        div_ex = 0.2
        Dsheild =[0.12,0.6]  # inner outer 0-1
        Dsheild =[]  # inner outer 0-1
        #trim = [0.62,0.65]  # trim fraction (in/out)
        trim = [0.81,0.81]  # trim fraction (in/out)
        sheild_connect=[0,1]
        sheild_base = -1
        #filename = '../eqdsk/Equil_AR3d1_bt_1d03li_0d8_Ipl_'+\
        #'20d25_XD_N_2LAT7R_v1_0.eqdsk'
        filename = '../eqdsk/Equil_AR3d1_XD_2015_Invess_5d0MA_EOF.eqdsk'
        filename = '../eqdsk/Equil_AR3d1_XD_2015_Invess_5d0MA_EOF_v2.eqdsk'
        filename = '../eqdsk/2015_XD_eqdk_2MJ4JW_v2_1.eqdsk'
        filename = '../eqdsk/Equil_AR3d1_16coils_v3_bt_'
        filename += '1d03li_0d8_Ipl_20d25_SF_NOEDDY_SOF_FINAL3_v9.eqdsk'
        filename = '../eqdsk/Equil_AR3d1_16coils_v3_bt_1d03li_'
        filename += '0d8_Ipl_20d25_SF_NOEDDY_SOF_FINAL3_v12.eqdsk'
        filename = '../eqdsk/Equil_AR3d1_16coils_v3_bt_1d03li_0d8_'
        filename += 'Ipl_20d25_XDext_NOEDDY_SOF_FINAL3_v11.eqdsk'
        filename = '../eqdsk/Equil_AR3d1_16coils_v3_bt_1d03li_'
        filename += '0d8_Ipl_20d25_XDext_v4_NOEDDY_SOF_FINAL3_v1.eqdsk'
        dataname = 'X'
        coils = {'internal':{'id':[],'dR':[],'dZ':[]},  # 12  0.85*2
                 'external':{'id':range(10,16),'dR':[],'dZ':[]}}  # list(range(0,13))
        targets['inner'] = {'L2D':[0.95],'open':False,  #0.85
                                 'graze':graze,'dR':-1}
        targets['outer'] = {'L2D':[1.92],'open':False,  # 1.9,1.95
                                 'graze':1.0*np.pi/180,'dR':0} #2.433 386
         
    elif config == 'Xex':
        dPlate = 0.2  # target plate length
        div_ex = 0.1
        Dsheild =[0.12,0.6]  # inner outer 0-1
        Dsheild =[]  # inner outer 0-1
        #trim = [0.62,0.65]  # trim fraction (in/out)
        trim = [0.78,0.7]  # trim fraction (in/out)
        sheild_connect=[0,1]
        sheild_base = 0.05
        #filename = '../eqdsk/Equil_AR3d1_bt_1d03li_0d8_Ipl_'+\
        #'20d25_XD_N_2LAT7R_v1_0.eqdsk'
        filename = '../eqdsk/Equil_AR3d1_XD_2015_Invess_5d0MA_EOF.eqdsk'
        filename = '../eqdsk/Equil_AR3d1_XD_2015_Invess_5d0MA_EOF_v2.eqdsk'
        filename = '../eqdsk/2015_XD_eqdk_2MJ4JW_v2_1.eqdsk'
        dataname = 'X'
        coils = {'internal':{'id':[],'dR':[],'dZ':[]},  # 12  0.85*2
                 'external':{'id':[],'dR':[],'dZ':[]}}  # list(range(0,13))
        targets['inner'] = {'L2D':[0.8],'open':False,'graze':graze,'dR':0.1}
        targets['outer'] = {'L2D':[3.12],'open':False,'graze':graze,'dR':0.11} #2.433
                    
    elif config == 'SN':
        TFopp = 'L'
        dPlate = 0.5  # target plate length
        div_ex = 0.25
        trim = [0.88,0.95]  # trim fraction (in/out)
        Dsheild =[0.2,0.55]  # inner outer 0-1
        Dsheild =[]
        sheild_connect=[0,1]
        sheild_base = -1
        #filename = '../eqdsk/Equil_AR3d1_bt_1d03li_0d8_Ipl'+\
        #'_20d25_EOF_2M56U8_v1_0.eqdsk'
        filename = '../eqdsk/2015_SN_eqdk_2MG9A4_v1_0.eqdsk'
        filename = '../eqdsk/2015_SN_eqdk_2MG9A4_v1_0_IDM.eqdsk'
        dataname = 'SND'
        coils = {'internal':{'id':[],'dR':[],'dZ':[]},
                 'external':{'id':[0,4],'dR':[],'dZ':[]}}  #list(range(0,11))
        targets['inner'] = {'L2D':[0.6],'open':False,
                                 'graze':graze,'dR':0.0}
        targets['outer'] = {'L2D':[0.65],'open':False,
                                 'graze':graze,'dR':0.0}
    
    elif config == 'SF':
        filename = '../eqdsk/Equil_AR3d1_bt_1d03li_0d8_Ipl_'+\
        '20d25_26co_2L2R9C_v1_0.eqdsk'
        dataname = 'SFD'
        coils = {'internal':{'id':[],'dR':[],'dZ':[]},
                 'external':{'id':list(range(0,10)),'dR':[],'dZ':[]}} 
        targets['inner1'] = {'L2D':[2],'open':False,'graze':graze}
        targets['inner2'] = {'L2D':[2],'open':False,'graze':graze}
        targets['outer1'] = {'L2D':[2],'open':False,'graze':graze}
        targets['outer2'] = {'L2D':[2],'open':False,'graze':graze} 
    
    elif config == 'SFm':
        TFopp = 'V'
        div_ex = 0.18
        dPlate = 0.35 # target plate length
        sheild_connect=[0,1]  # 0.3,0.95
        #Dsheild =[0,0.84]  # inner outer 0-1
        Dsheild =[]  # inner outer 0-1
        trim = [0.75,0.7]  # trim fraction (in/out)
        
        filename = '../eqdsk/Equil_AR3d1_16coils_SFminus_v4_2015'+\
        '_09_bt_1d03li_0d8_Ipl_20d25_SOF.eqdsk'
        #filename = '../eqdsk/2015_SFminus_eqdsk_2MCQRZ_v1_0_IDM.eqdsk'

        dataname = 'SFm'
        coils = {'internal':{'id':[],'dR':[],'dZ':[]},
                 'external':{'id':list(range(0,16)),'dR':[],'dZ':[]}} 
        targets['inner1'] = {'L2D':[1.1+0.52],'open':True,
                                  'graze':graze,'dR':0.0}
        targets['inner2'] = {'L2D':[1.2-0.08],'open':False,
                                  'graze':graze,'dR':-1}
        targets['outer1'] = {'L2D':[1.65-0.7],'open':False,
                                  'graze':graze,'dR':-1}
        targets['outer2'] = {'L2D':[1.1+0.15],'open':True,
                                  'graze':graze,'dR':0.0}  

    elif config == 'SFp':
        TFopp = 'V'
        div_ex = 0.18
        dPlate = 0.35  # target plate length
        sheild_connect=[0.22,1]
        #Dsheild =[0,0.84]  # inner outer 0-1
        Dsheild =[]  # inner outer 0-1
        trim = [0.75,0.75]  # trim fraction (in/out)
        
        filename = '../eqdsk/Equil_AR3d1_16coils_SFplus_v4_2015'+\
        '_09_bt_1d03li_0d8_Ipl_20d25_SOF.eqdsk'

        dataname = 'SFp'
        coils = {'internal':{'id':[],'dR':[],'dZ':[]},
                 'external':{'id':list(range(0,10)),'dR':[],'dZ':[]}} 
        targets['inner1'] = {'L2D':[1.1],'open':True,'graze':graze,'dR':0.1}
        targets['inner2'] = {'L2D':[1.2],'open':False,'graze':graze,'dR':-1}
        targets['outer1'] = {'L2D':[1.65],'open':False,'graze':graze,'dR':-1}
        targets['outer2'] = {'L2D':[1.1],'open':True,'graze':graze,'dR':0.1}  
        
    keys = ['conformal','dSOL','Nsol','dRfw','div_ex','TFopp',
             'sheild_base','tfw','tBBsupport','BS','BBfrac','BB',
             'sheild','VV','graze','dPlate','nTF','Jmax','trim',
             'config','targets']
    config = {}
    for key in keys():
        config[key] = eval(key)
    print(config)
        
    
                    
                 
        