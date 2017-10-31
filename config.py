"""
File to contain values and data that I'm going to be calling in multiple files, and is just easier to keep here.
"""
#from astropy.time import Time

lpos_list= [56132, 57063] #mjd dates of position changes
#gain_change_list=Time( ['2009-05-11', '2009-08-12', '2011-03-08', '2012-03-26', '2012-07-23', '2013-06-24', '2014-07-21', '2014-11-03', '2015-02-09'], scale='utc')
gain_change_list= ['2009-05-11', '2009-08-12', '2011-03-08', '2012-03-26', '2012-07-23', '2013-06-24', '2014-07-21', '2014-11-03', '2015-02-09']
#gain_change_list_mjd = gain_change_list.mjd
wave_limit_list =  [[1130,1430],[1130,1800],[1130,1850],[1130,1900]]
lyman_mask= [1208, 1225]
oxygen_mask= [1295, 1312] #airglow wavelengths to be filtered according to lightcurve

silicon_lines= [1190, 1193, 1195, 1197, 1207, 1260, 1265, 1304, 1309] #silicon lines in the UV
