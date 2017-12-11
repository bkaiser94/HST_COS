"""
File to contain values and data that I'm going to be calling in multiple files, and is just easier to keep here.
These values are not called by lightcurve.cos.py, so you have to manually make sure the masks match up for the airglow values (the wavelength limits are called by the script that calls cos.py, so you're fine there.)
These are the masking values in the corresponding spectra that get plotted. The actual light curve masking values are located in cos.py, and have to be separately updated manually at the moment.
"""
#from astropy.time import Time
#booleans for whether or not to actually do the masking of the different lines
do_lyman = True #this one just does a much smaller mask when disabled. (It would be the only thing plotted otherwise)
do_oxygen= True #smaller mask when False (good for G130M setup)
do_nitrogen= True
cos_refresh_rate = 0.032 #seconds
lpos_list= [56132, 57063] #mjd dates of position changes
#gain_change_list=Time( ['2009-05-11', '2009-08-12', '2011-03-08', '2012-03-26', '2012-07-23', '2013-06-24', '2014-07-21', '2014-11-03', '2015-02-09'], scale='utc')
gain_change_list= ['2009-05-11', '2009-08-12', '2011-03-08', '2012-03-26', '2012-07-23', '2013-06-24', '2014-07-21', '2014-11-03', '2015-02-09']
#gain_change_list_mjd = gain_change_list.mjd
#wave_limit_list =  [[1130,1430],[1130,1800],[1130,1850],[1130,1900]] #old one as of 2017-11-06
#wave_limit_list =  [[930,1430], [1105,1232], [1172,1430],[1172,1800]] #new ones. Masks the weird bump in the ~1140 range with the G140L, and allows full spectra for G130M (new as of 2017-11-06)
#wave_limit_list =  [[1172,1800],[1150,1280],[1300,1435],[1300,1500],[1500,1700],[1105,1200]] #new ones as of 2017-11-27
#wave_limit_list =  [[1180,1300],[1300,1500],[1500,1700],[1180,1800],[1250,1271],[1150,1200],[1300,1425],[1150,1425],[1150,1271],[1452, 1570],[1615,1760]] #new ones as of 2017-11-28
G160M_1600A_wave_list= [[1452, 1570],[1615,1760]]
G130M_1291A_wave_list= [[1180,1800],[1250,1271],[1150,1200],[1300,1425],[1150,1425],[1150,1271]]
G130M_1300A_wave_list= G130M_1291A_wave_list #these are the same for now
G140L_1105A_wave_list= [[1180,1300],[1300,1500],[1500,1700],[1180,1800]]
G130M_1096A_wave_list= [[955, 1060],[1106,1232]]
wave_limit_list= G160M_1600A_wave_list+G130M_1291A_wave_list+ G130M_1300A_wave_list + G140L_1105A_wave_list+G130M_1096A_wave_list

#wave_limit_list =  [[1250,1279]] #broad lines in target
#wave_limit_list= [[1150, 1425]]
if do_lyman:
    lyman_mask= [1206, 1226]
else:
    lyman_mask= [1214, 1217]
if do_oxygen:
    oxygen_mask= [1295, 1313] #airglow wavelengths to be filtered according to lightcurve
else:
    oxygen_mask= [1298,1306]
if do_nitrogen:
    nitrogen_mask= [1195, 1207] #hopefully this helps the dimmer targets... hopefully
else:
    nitrogen_mask= [1200,1200]


silicon_lines= [1190, 1193, 1195, 1197, 1207, 1260, 1265, 1304, 1309] #silicon lines in the UV
