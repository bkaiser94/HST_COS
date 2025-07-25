"""
File to contain values and data that I'm going to be calling in multiple files, and is just easier to keep here.
local_lightcurve.cos.py does actually look here now for the masking values, so no need to update anywhere other than here (unless you already made the lightcurves before updating this file.)
"""

#refpath= 'HST_COS_ref_files/' #path to be used by local_lightcurve/utils.py for expand_refname()
refpath='/Users/BenKaiser/Desktop/HST_COS/HST_COS_ref_files/' #Contains the reference files for both STIS and COS
#degree of masking level: 0 means no mask, 1 means small mask (usually G130M setting), 2 means large mask (G140L setting)
mask_deg=1





cos_refresh_rate = 0.032 #seconds
lpos_list= [56132, 57063, 58028] #mjd dates of lifetime position changes (added LP-POS 4)
#gain_change_list=Time( ['2009-05-11', '2009-08-12', '2011-03-08', '2012-03-26', '2012-07-23', '2013-06-24', '2014-07-21', '2014-11-03', '2015-02-09'], scale='utc')
gain_change_list= ['2009-05-11', '2009-08-12', '2011-03-08', '2012-03-26', '2012-07-23', '2013-06-24', '2014-07-21', '2014-11-03', '2015-02-09']
#gain_change_list_mjd = gain_change_list.mjd

lyman_mask_list = [
    [0,0],
    [1214,1217],
    [1206,1226]]

oxygen_mask_list= [
    [0,0],
    [1301,1308],
    [1295,1313]]

nitrogen_mask_list = [
    [0,0],
    [1198,1202],
    [1195,1207]]

lyman_mask= lyman_mask_list[mask_deg]
oxygen_mask= oxygen_mask_list[mask_deg]
nitrogen_mask= nitrogen_mask_list[mask_deg]


G160M_1600A_wave_list= [[1452, 1570],[1615,1760],[1452,1760]]
G160M_1577A_wave_list= [[1387,1557],[1578,1742],[1387,1742]] #Based only on the GW-Librae observation which only included 4 exposures all in FP-POS 3
G130M_1222A_wave_list = [[1077,1201],[1233,1354],[1077,1354],[1250,1271]]
#G130M_1291A_wave_list=[[1142, 1268],[1298,1422],[1142,1422], [1142,1200],[1250,1268]]
G130M_1291A_wave_list=[[1142, 1268],[1298,1422],[1142,1422], [1142,1200],[1250,1268],[1354,1360]] #added one for G227-5
#G130M_1291A_wave_list=[ [1214,1217], [1301,1308], [1198,1202]] #done to track geocoronal emission separately. Need mask deg 0

G130M_1300A_wave_list= [[1250,1271],[1150,1200],[1300,1425],[1150,1425],[1150,1271]]
G130M_1309A_wave_list = [[1162,1286],[1318,1442],[1162,1442],[1250,1271], [1162,1200]]
G130M_1327A_wave_list= [[1182,1307],[1335,1460],[1182,1460]] 
G140L_1105A_wave_list= [[1180,1300],[1300,1500],[1500,1700],[1180,1800]]
G130M_1096A_wave_list= [[955, 1075],[1106,1232], [955,1232]]
G230L_2950A_wave_list=[[1715,3118],[1715,1922],[2805,3118]]


wave_limit_list= G160M_1600A_wave_list+G130M_1291A_wave_list+ G130M_1300A_wave_list + G140L_1105A_wave_list+G130M_1096A_wave_list + G160M_1577A_wave_list+G130M_1327A_wave_list+G130M_1309A_wave_list+G130M_1222A_wave_list+G230L_2950A_wave_list



stis_dict={
    '':''
    
    
    
    
    }




silicon_lines= [1190, 1193, 1195, 1197, 1207, 1260, 1265, 1304, 1309] #silicon lines in the UV

emission_lines=[
    ['C III', 1175],
    ['N V', 1242],
    ['C II', 1335],
    ['Si IV', 1400],
    ['C IV', 1550]] # from page 5 of HST/COS II

seg_gap_dict={'G130M':{
    '1096': [1072,1104],
    '1222':[1201,1233],
    '1291':[1268,1298],
    '1300':[1278,1306],
    '1309':[1286,1318],
    '1096':[1075,1106],
    '1327':[1307,1335]},
             'G140L':{'1105': [1179,1180]}, 
             'G160M': {'1600': [1569,1619],
                       '1577': [1557,1578]},
             'G230L':{
                 '2950':[1922,2805]}
             }
