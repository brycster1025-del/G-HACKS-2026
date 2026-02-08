import novatel_edie as ne

import numpy as np





#TIME_LOG="#TIMEA,USB1,0,50.5,FINESTEERING,2209,515163.000,02000020,9924,16809;VALID,-2.501488425e-09,6.133312031e-10,-17.99999999630,2022,5,13,23,5,45000,VALID*1100ad64\r\n"
#BESTPOS_LOG='#BESTPOSA,USB1,0,58.5,FINESTEERING,2209,502061.000,02000020,cdba,16809;SOL_COMPUTED,PPP,51.15043706870,-114.03067882331,1097.3462,-17.0001,WGS84,0.0154,0.0139,0.0288,"TSTR",11.000,0.000,43,39,39,38,00,00,7f,37*52483ac5\r\n'
filename="RTK_Data2"# Put file path containing log here
parser=ne.FileParser(filename)


best_lat=0 #initialize all variables
best_long=0
best_hgt=0
best_lat_dev=100
best_long_dev=100
best_hgt_dev=100
for log in parser:
    if isinstance(log,ne.messages.BESTPOS):
        lat=log.latitude
        long=log.longitude
        height=log.orthometric_height
        hgt_std_dv=log.height_std_dev
        long_std_dev=log.longitude_std_dev
        lat_std_dev=log.latitude_std_dev
        if lat_std_dev<= best_lat_dev:
            best_lat=lat
            best_lat_dev=lat_std_dev
        if long_std_dev<=best_long_dev:
            best_long=long
            best_long_dev=long_std_dev
        if hgt_std_dv<= best_hgt_dev:
            best_hgt=height
            best_hgt_dev=hgt_std_dv
print("Best individual coordinates in the geodetic frame are (lat,long,height)" ,best_lat,best_long,best_hgt)
print("Corresponding standard deviations are",best_lat_dev,best_long_dev,best_hgt_dev)

