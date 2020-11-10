# MAJOR UPDATE (AUGUST 14): We were doing median fits twice, this has been corrected in this latest version.

# run this as " python TransientHunter.py runlist" where 'runlist' contains a list of runs.


from time import time
from datetime import datetime

import os

##get_ipython().magic(u'cd /raid/biggams/tsenlin/')


# In[1]:

import numpy as np
import matplotlib.pyplot as plt
#import iminuit
#import probfit
import PedVarAnalysisModule as pva
#from scipy.signal import gaussian , general_gaussian

# Load modules related to SQL.
import pymysql
import pandas as pd
import tables



import sys
#sys.path.append('/homes/jivaro/gabrielc/PyVBF/')


#Load Astropy modules
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS, Galactic, FK4, FK5
from astropy.vo.client import conesearch

from astroquery.simbad import Simbad

# import stuff to make patches
import matplotlib
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches

class connect_to_VDB:
    def __init__(self):
        pass
    def __enter__(self):
        self.con = pymysql.connect(host= 'romulus.ucsc.edu',
                                user='readonly',
                                db  ='VERITAS',
                                password='')
        return self.con
    def __exit__(self,*args):
        self.con.close()



# Load in list of runs.
try:
    File = open(sys.argv[1],'r')

except:
    print ("File " + sys.argv[1] + " not found.")
    sys.exit(37)

PicTempNames = []

for ii in File.readlines() :
    PicTempName = ii.split()
    PicTempNames.append(PicTempName)
    
    
from IPython.core.pylabtools import figsize
figsize(15,10)
plt.rcParams['axes.labelsize']   = 25
plt.rcParams['axes.titlesize']   = 30
plt.rcParams['xtick.labelsize']  = 15
plt.rcParams['ytick.labelsize']  = 15
plt.rcParams['legend.fontsize']  = 30
plt.rcParams['lines.markersize'] = 15 



#Load list of run ids and l3 rates.
id_vs_l3 = np.loadtxt('runid_to_l3rate')
# generate dictionary that maps runID to l3 rate
id_to_l3 = {str(int(id_vs_l3[ii][0])) : int(id_vs_l3[ii][1]) for ii in range(0,len(id_vs_l3[:,1]))}

#Set SIMBAD parameters

customSimbad = Simbad()

#Modify simbad so that it tells  the RA&DEC and magnitudes of stars nearby a certain pixe.

customSimbad.add_votable_fields('ra(deg)','dec(deg)','flux(B)','flux(V)')
#customSimbad.add_votable_fields('flux(B)','flux(V)')
customSimbad.remove_votable_fields('coordinates')


# ## Set pixel coordinates on camera plane and set location of Mt.Hopkins


channel_num = 500
coords =  np.zeros((channel_num,3))



# get pixel coordinates from file
file = open('Pixel_Coordinates','r')
lines = file.readlines()[1:] #
for line in lines:
	line = line.strip()
	columns = line.split(  )
#	print line
#	print columns[1]
#	print coords[int(columns[0])][1]
	coords[int(columns[0])-1][0] = (columns[0])
	coords[int(columns[0])-1][1] = (columns[1])
	coords[int(columns[0])-1][2] = (columns[2])



# Now make another list of pixel coordinates in the polar coordinate system

polar_coords = np.zeros((500,3))
polar_coords[:,0] = coords[:,0]-1
polar_coords[:,1] = np.sqrt(coords[:,1]**2+coords[:,2]**2)
polar_coords[:,2] = np.arctan2(coords[:,2],coords[:,1])*180/np.pi

mt_hopkins = EarthLocation(lat=31.675011*u.deg, lon=-110.952293*u.deg, height=1270*u.m)


# In[3]:

def GetRunInfo(RunID, verbose = False):
    sql_cmd_t = '''SELECT info.source_id,info.data_start_time, info.data_end_time,
       obs.ra, obs.decl, obs.epoch, info.offsetRa, info.offsetDEC
       FROM tblRun_Info AS info, tblObserving_Sources as obs
       WHERE info.run_id = {} AND obs.source_id = info.source_id '''.format(RunID)

    with connect_to_VDB() as con:
        table = pd.read_sql_query(sql_cmd_t,con)

    source = table['source_id']
    time_start = table['data_start_time']
    time_end = table['data_end_time']
    ra = table['ra']
    dec = table['decl']
    offset_ra = table['offsetRa']
    offset_dec = table['offsetDEC']
    #print table


    # Get telescopes pointing in radians.
    tel_time_start = time_start.values
    tel_time_end = time_end.values
    tel_RA = ra.values + offset_ra.values
    tel_DEC = dec.values + offset_dec.values

    if verbose == True:
        print (source.values[0], RunID ,tel_time_start[0], tel_RA[0]*180/np.pi, tel_DEC[0]*180/np.pi)
    return  tel_time_start[0], tel_time_end[0] , tel_RA[0]    , tel_DEC[0]

def StarAltAz_StartEnd(star_RA = 0, star_DEC = 0, runID = 59823):
    '''Star AltAZ StartEnd takes in RA and DEC coordinates and gives you the corresponding coordinates in the camera plane at the beginning, middle, and at the end of a run.

    returns star_camera_map_start ,star_altaz_middle, star_altaz_end
    '''
    time_start, time_end, telRA, telDEC = GetRunInfo(runID)
    time_middle = time_start + np.timedelta64(10,'m')
    time_middle = Time(str(time_middle))  
    time_start = Time(str(time_start))


    time_end = Time(str(time_end))



    sky_coordinates = SkyCoord('{}  {} '.format(telRA , telDEC), unit = (u.rad , u.rad ) )
    centre_altaz = sky_coordinates.transform_to(AltAz(obstime=time_start, location = mt_hopkins) )

    star = SkyCoord(star_RA ,star_DEC , unit="deg")
    star_altaz = star.transform_to(AltAz(obstime=time_start, location = mt_hopkins) ) 
    star_camera_map_start = [-np.cos(star_altaz.alt.rad)*((star_altaz.az.deg)  -(centre_altaz.az.deg)),-( (star_altaz.alt.deg) - (centre_altaz.alt.deg))]


    centre_altaz = sky_coordinates.transform_to(AltAz(obstime=time_middle, location = mt_hopkins) )
    star_altaz = star.transform_to(AltAz(obstime= time_middle , location = mt_hopkins) ) 
    star_camera_map_middle = [-np.cos(star_altaz.alt.rad)*((star_altaz.az.deg)  -(centre_altaz.az.deg)),-( (star_altaz.alt.deg) - (centre_altaz.alt.deg))]



    centre_altaz = sky_coordinates.transform_to(AltAz(obstime=time_end, location = mt_hopkins) )
    star_altaz = star.transform_to(AltAz(obstime= time_end , location = mt_hopkins) ) 
    star_camera_map_end = [-np.cos(star_altaz.alt.rad)*((star_altaz.az.deg)  -(centre_altaz.az.deg)),-( (star_altaz.alt.deg) - (centre_altaz.alt.deg))]


    return star_camera_map_start,star_camera_map_middle ,star_camera_map_end


def StarFinder(run_num):
    time_start, time_end, telRA, telDEC = GetRunInfo(run_num)
    channels_w_stars = []

    yearmonthtime = time_start.astype('M8[D]')
    yearmonthtime = str(yearmonthtime)


    #get current and voltage values...
    print('Loading current samples...')
    scale = 5
    print('done!')

    # grab all stars falling inside the camera plane

    print('Finding stars in the camera FOV...')
    sky_coordinates = SkyCoord('{}  {} '.format(telRA , telDEC), unit = (u.rad , u.rad ) )
    try:		
        result = customSimbad.query_region(sky_coordinates, radius='1.75 degrees')
    except KeyboardInterrupt:
        quit()
    except:
        print ("error!")
        #continue
    result = np.array(result)
    print('done!')
    #print(result)
    # Go through every star in the FOV above a certain brightness.
    try:
		for ii in range(0,len(result)): # iterates over stars.

		    if result[ii][3] < 10: # filters on  blue band magnitude less than 10.

		#this is a workaround exclusive to reedbuck if using python 3, the first entry is  ''b'Star_name''' rather than just ''Star_name''
		#				original_result = str(np.copy(result[ii][0]))
		#				result[ii][0] = original_result
		#				result[ii][0] = str(result[ii][0][2:-1])

		        starcoords = customSimbad.query_object('{}'.format(result[ii][0]))
		        starcoords_start, starcoords_middle, starcoords_end =StarAltAz_StartEnd(star_RA =starcoords[0][1],   star_DEC = starcoords[0][2],  runID= run_num)

		        xstar_start = starcoords_start[0]
		        ystar_start =starcoords_start[1]
		        rstar_start = np.sqrt(xstar_start**2+ystar_start**2)
		        thetastar_start = np.arctan2(ystar_start,xstar_start)*180/np.pi

		        xstar_end = starcoords_end[0]
		        ystar_end =starcoords_end[1]
		        rstar_end = np.sqrt(xstar_end**2+ystar_end**2)
		        thetastar_end = np.arctan2(ystar_end,xstar_end)*180/np.pi

		        xstar_middle = starcoords_middle[0]
		        ystar_middle =starcoords_middle[1]
		        rstar_middle = np.sqrt(xstar_start**2+ystar_start**2)
		        thetastar_middle = np.arctan2(ystar_start,xstar_start)*180/np.pi

		        vector1 = np.array(starcoords_end) - np.array(starcoords_start) 
		        vector2 = np.array(starcoords_middle) - np.array(starcoords_start)

		        # Check only the channels that fall within the angle of the star at the beginning and the end of the run
		        # Make sure to adjust which angle (initial or final) is the upper or lower range, this depends on wether the arc is clockwise or counterclockwise.
		        crossprod = np.cross(vector1,vector2)

		        # Check the channels that fall near the radious of the star ( radius from center)
		            # The radius range is scaled about 0.045 degrees, which ensures that a star with a 68%psf radius of 0.03 will fall inside of a pixel if    r_star - pixel < pixel_r - pixelPSF
		        try:
		            if 3 < result[ii][3] < 7:
		                radius_range = polar_coords[:,0][ 0.045*5 >= np.abs(rstar_start - polar_coords[:,1])]
		            elif  result[ii][3] < 3:
		                radius_range = polar_coords[:,0][ 0.15*3 >= np.abs(rstar_start - polar_coords[:,1])]
		            elif 7 <= result[ii][3] < 9:
		                radius_range = polar_coords[:,0][ 0.045*3 >= np.abs(rstar_start - polar_coords[:,1])]
		            elif  9 <= result[ii][3]:
		                radius_range = polar_coords[:,0][ 0.045*0.8 >= np.abs(rstar_start - polar_coords[:,1])]
		            # Check only the channels that fall within the angle of the star at the beginning and the end of the run


		            if crossprod >= 0:
		                #print('clockwise path!') # If clockwise use end angle as upper limit, start angle as lower limit
		                if thetastar_end > 90 and thetastar_start < -90:  #this accounts for the case where the star crosses 3rd quadrant to 2nd.           
		                    theta_range = np.where((polar_coords[:,2] >= thetastar_end) | (polar_coords[:,2] <= thetastar_start) )
		                else:
		                    theta_range = np.where((polar_coords[:,2] >= thetastar_end) & (polar_coords[:,2] <= thetastar_start) )
		            else:
		                #print('counter clockwise path!') # If counterclockwise, flip limits.
		                if thetastar_end < -90 and thetastar_start > 90:  #this accounts for the case where the star crosses 2nd quadrant to 3rd.  
		                    theta_range = np.where((polar_coords[:,2] <= thetastar_end) | (polar_coords[:,2] >= thetastar_start) )
		                else:
		                    theta_range = np.where((polar_coords[:,2] >= thetastar_start) & (polar_coords[:,2] <= thetastar_end) )



		############
		            #theta_range = np.where((polar_coords[:,2] > thetastar_end) & (polar_coords[:,2] < thetastar_start) )
		            # Create a list of pixels that are near the star's radius and within the angles.
		            # radius_range has all the pixels that fall inside the search radius of every pixel in a circle. theta_range has all the pixels that fall inside the angular arc of te star.
		            #calibration_channels
		            calibration_channels =np.intersect1d(radius_range,theta_range)
		        except ValueError:
		            pass

		        if len(calibration_channels) > 0: #Check that the star fell entirely within some channels
		#					print('Star {} with blue magnitude {} falls within channels {}'.format( result[ii][0], result[ii][3],calibration_channels[ calibration_channels <= 459 ]))


		            for jj in calibration_channels: 
		                jj = int(jj)
		#							jj = int(jj)

		#					plt.errorbar(x,y[jj] , y_err[jj], ms = 8, fmt = 'o')

		                channels_w_stars.append(jj)


    except TypeError:
        pass

    return np.array(list(set(channels_w_stars))), yearmonthtime


# In[78]:




def RunData(filename, channel_num = 499, ):

    f = PyVBF.PyVBFreader(filename) # reads the file.
    RunNum = filename[-10:-5]
    print 'Opening Run number',  RunNum
    numPacket  = f.get_numPackets() #Get the number of packets in the cvbf file.
    print "Total packets " , numPacket

    charge1 = charge2 = charge3 = charge4 = np.zeros((channel_num,numPacket))
    #	charge1 = charge2 = charge3 = charge4 = np.zeros((channel_num,(numPacket-1)/2))
    #for ii in range(1,numPacket):  # Go through one bin at a time, (bins are called 'slices' in this code).
    time1 = time()

    for ii in range(1,numPacket):
        f.go_to_packet(ii)

        try:
            f.loadEvent(0) # Loads a single telescope
    #		event_type = f.getRawEventType()
            samples = f.getSamples() # This gets all the samples stored in every pixel of a telescope, for a single event.
            charge1[:,ii] = np.var(samples[:,:],axis=1) ## add the variance.
            del samples
        except Exception:
            charge1[:,ii] = int(0)


    list = [ x for x in range(0,500)]
    charge1[:,0] = list

    del f
    time2 = time()

    print 'time passed = {} seconds'.format(time2-time1)
    return charge1


# Load data

def readxy(Filename):
    file = open(Filename, 'r')
    Column0 = []
    Column1 = []
##	Column2 = []
    for line in file:
        line = line.strip()
        columns = line.split( )
        Column0.append(float(columns[0]))
        Column1.append(float(columns[1]))
##		Column2.append(float(columns[2]))
    return np.array(Column0), np.array(Column1),



def GetCurrentVoltage(RunID):
    # Get currents from here.
    RunID = RunID
    #51692 Note: airplane/satellite at 16:20ish in. Affected L3.
    # 52166 Note: Spike at 4 minuotes from UFO.
    sql_cmd_t = '''SELECT hv.db_start_time,hv.channel, hv.current_meas, hv.voltage_meas
       FROM tblHV_Telescope0_Status AS hv,tblRun_Info AS info
       WHERE info.run_id = {} AND hv.db_start_time >= info.data_start_time
       AND hv.db_start_time BETWEEN info.data_start_time AND info.data_end_time
       AND hv.db_end_time   BETWEEN info.data_start_time AND info.data_end_time'''.format(RunID)

    current = np.zeros((499,80))
    voltage = np.zeros((499,80))

    with connect_to_VDB() as con:
        table = pd.read_sql_query(sql_cmd_t,con)
    for ii in range(0,499):
        chan = table[ table['channel'] == ii+1]
        tmp= chan['current_meas']
        current[ii,:len(chan['current_meas'].values)] = chan['current_meas'].values
        voltage[ii,:len(chan['voltage_meas'].values)] = chan['voltage_meas'].values

    return current, voltage , len(chan['current_meas'].values) , len(chan['voltage_meas'].values)





# Go through all the points in the graph. At every point, grab the next "window_size" points and the previous "window_size" points
# Take the median of that set and that becomes your new point. This should wash out long transients like stars and keep shorter ones.

def median_fit(file_data,window_size = 10):
    '''MULTI-ARRAYS ONLY.
        Go through all the points in the graph. At every point, grab the next "number_of_points  = window_size" points and the previous " number_of_points =window_size" points
     Take the median of that set and that becomes your new point. This should wash out long transients like stars and keep shorter ones.
    '''
    WS = 10
    window_size = WS

    step_size = 1

    data = np.copy(file_data)
#    new_vals = np.copy(scharge[:,0])

    max_step = len(data[0])-window_size # Position of final point.
    step_num = np.floor_divide(len(data[0]),step_size)  # total number of steps.
    median_data = np.zeros((500,len(data[0])))


    for step_number in np.arange(0,len(data[0])-window_size,step_size):
    # Make sure you go up the last point you can do the median calculation
        if step_number - window_size >=0:
            median_data[:,step_number] = (np.median(data[:,step_number-window_size:step_number+window_size], axis = 1))
        else:
            median_data[:,step_number] = 0

    median_removed_data = data - median_data
    return median_data, median_removed_data





def SIMBADTransientHunter(filter_var,filter_err_var, charge, channel_list ,  window_size = 200, WS = 200, telID = 1):

    '''SIMBADTransientHuner analyses channels without pathologies. It will save a channel's signal  if
    that channel as  at least 3 points above 3-sigma that are no more than 3 bins away from the middle of that region. 
     It does not look at channels with stars in them, these channels are vetoed out.

     '''



    window_size = WS
#  Apply a de-trend to every channel. This removes the slow-varying shape of the signal and keeps a baseline. Window size is set to 200 so that it guarantees that the de-trend will not remove a transient from the baseline signal.
    median_data,baseline_data = median_fit(filter_var,window_size=WS)
    # Create a list of the lower bound	, that is we want to make sure that the lower bound of a bin value is above 3-sigma
    base_var_minus_errvar = baseline_data[:, WS : -WS ] - filter_err_var[:, WS : -WS ]


    for jj in channel_list:
            trigger_on = False
            if np.all(I[jj][:LI] > 0) == True and np.min(filter_var[jj] < 100): # Checks that the current is not dead at any time in that pixel.


                # Try to find points that are above 3-sigma (in their lower bound) and that are also close together. All points above 3-sgima are saved to "med_pick_of_the_litter")
                med_pick_of_the_litter = np.argwhere(np.abs(base_var_minus_errvar[jj])>=3*np.std(baseline_data[jj ,WS: -WS ]))



                neighbour_check_array = []

                for ii in med_pick_of_the_litter[:,0]: # go through the dataset containing the bin number of all points above 3-sigma.
                    distance_check_array = np.array([x - ii for x in med_pick_of_the_litter[:,0]]) #iterate over the array, each iteration subtracts the ii'th entry. Entries close to that bin will turn to numbers like 3,2 or 1.
                    neighbour_check_array = distance_check_array[ (np.abs(distance_check_array) > 0) & (np.abs(distance_check_array) <= 3)] #Check that there are entries close together.

                    if len(neighbour_check_array) >= 2 and trigger_on == False: # Check that there's at least three entries close together. 'trigger_on' checks that a transient's already been found on this channel.
                        trigger_on = True # flip the value of 'trigger_on'

                        print('Found Something in Channel {}, saving data to ''TransientRuns'' folder'.format(jj))
                        if not os.path.exists('TransientRuns/{}/{}'.format(date_run_on, runID)):
                            os.makedirs('TransientRuns/{}/{}'.format(date_run_on, runID) )
                        np.savetxt('TransientRuns/{}/{}/Run{}_Tel{}_Ch{}_TRANSIENT_CHARGE.txt'.format(date_run_on,runID,  runID,  telID ,jj), np.around(charge[jj], decimals =2) ,fmt='%3.2f')


            else:
                print('HV was turned off in channel {} , we will not search there.'.format(jj))



##################################


def SIMBADTransientHunter2(filter_var, filter_err_var, charge, channel_list, rundate ,  window_size = 200, WS = 200, telID = 1, fake_l3_rate = 200, bin_size = 10, runID = 0):

    '''SIMBADTransientHunter2 is more advanced than the first. It keeps track of channels with transients and keeps a log of every run/telescope/channel with a transient. analyses channels without pathologies. It will save a channel's signal  if either conditon is true:
    1. that channel as  at least 3 points above 3-sigma that are no more than 3 bins away from the middle of that region. 
     It does not look at channels with stars in them, these channels are vetoed out.

    returns T1_on, T1_off, T1_peaks
     '''
    try:
        fake_l3_rate = id_to_l3[runID]
        print('l3 rate is ',fake_l3_rate)
    except KeyError:
        pass

    T1_on = []
    T1_off = []

    T1_peaks =  np.zeros(500)
    T1_all = [ii for ii in range(0,499)]


    window_size = WS
#  Apply a de-trend to every channel. This removes the slow-varying shape of the signal and keeps a baseline. Window size is set to 200 so that it guarantees that the de-trend will not remove a transient from the baseline signal.
    median_data,baseline_data = median_fit(filter_var,window_size=WS)
    # Create a list of the lower bound	, that is we want to make sure that the lower bound of a bin value is above 3-sigma
    base_var_minus_errvar = baseline_data[:, WS : -WS ] - filter_err_var[:, WS : -WS ]


    for jj in channel_list:
            trigger_on = False
            if np.all(I[jj][:LI] > 0) == True and np.min(filter_var[jj] < 100): # Checks that the current is not dead at any time in that pixel.


                # Try to find points that are above 3-sigma (in their lower bound) and that are also close together. All points above 3-sgima are saved to "med_pick_of_the_litter")
                med_pick_of_the_litter = np.argwhere(np.abs(base_var_minus_errvar[jj])>=3*np.std(baseline_data[jj ,WS: -WS ]))



                neighbour_check_array = []

                for ii in med_pick_of_the_litter[:,0]: # go through the dataset containing the bin number of all points above 3-sigma.
                    distance_check_array = np.array([x - ii for x in med_pick_of_the_litter[:,0]]) #iterate over the array, each iteration subtracts the ii'th entry. Entries close to that bin will turn to numbers like 3,2 or 1.
                    neighbour_check_array = distance_check_array[ (np.abs(distance_check_array) > 0) & (np.abs(distance_check_array) <= 3)] #Check that there are entries close together.

                if len(neighbour_check_array) >= 2 and trigger_on == False: # Check that there's at least three numbers close together.
                    # If this check passes, congratulations... you found a transient.
                    trigger_on = True
                    print('Found Something in Tel {} Channel {}, saving data to ''TransientRuns'' folder'.format(telID, jj))

                    index1 = np.argwhere( (np.abs(distance_check_array) > 0) & (np.abs(distance_check_array) <= 3)) #

                        #Find the bin numbers of the bins that triggered the Transient Hunter.

                    trigger_bins = [int(med_pick_of_the_litter[kk]) for kk in index1]
                    peak_vars = filter_var1[jj,trigger_bins[0]+WS-10:trigger_bins[-1]+WS+10]
                    try:
                        peak_bin = int(np.argwhere(filter_var1[jj] == np.max(peak_vars)))
                        peak_value = np.max(peak_vars)
                        #convert bin number to approximate event, and then to approximate timing. Note that using the L3 rate to find the true timing of transients such as shooting stars will undersell how fast the transient really is...
                        peak_approx_event = (peak_bin*bin_size*3)
                        peak_approx_seconds = peak_approx_event/fake_l3_rate
        





                        file = open('TerrestrialTransientLog','a')
                                    # File labels go as follows:  date, Run ID,  telID, channel number, largest variance value, bin number of largest variance, approximate event number, approximate seconds.
                        file.write('{} {} {} {} {} {} {} {} \n'.format(rundate,runID, telID, int(jj), peak_value ,peak_bin,peak_approx_event,peak_approx_seconds ))
                        file.close()

                        T1_on.append(jj)
                        T1_peaks[jj] = peak_bin
                    except TypeError:
                        pass
    #Obtain a list of all channels that did not trigger the transient hunter.
            else:
                print('HV was turned off in channel {} , we will not search there.'.format(jj))
    T1_off = list(set(T1_all)  -set(T1_on))
    return T1_on,T1_off,T1_peaks



def Map_UFO(data_on, data_off, data_peaks,date_run_on = 'd19921124', tel_num = 1, l3 = 200, run_num = 0, set_filter = False ):
    if len(data_on) > 0:
        # see if you can get the l3 rate from the list. If that fails it remains set to 200 Hz which is approximately right for V5 runs.
        try:
            l3 = id_to_l3[run_num]
        except KeyError:
            pass
        print('l3 rate is {}'.format(l3))
        pixels = np.loadtxt('pixels.csv')
        tel_1 = pixels[pixels[:,0] ==1]

        plt.close()
        plt.clf()


        data_on = np.array(data_on)
        data_off = np.array(data_off)
        data_peaks = np.array(data_peaks)


        fig, ax = plt.subplots()

        timing_array = (data_peaks[data_peaks > 0]- np.min(data_peaks[data_peaks > 0 ]))*30/l3




        if set_filter == True and len(data_on) >= 3:
            # Throw out transient events that acurred 500 seconds after the first event
            kick_out_array = (np.argwhere(timing_array > 50))
            kick_out_array = kick_out_array.flatten()
            Npeaks = np.delete(timing_array,kick_out_array)
            Non = np.delete(data_on,kick_out_array)
            Noff = np.append(data_off,data_on[kick_out_array])

            Noff =np.sort(Noff)

            data_on = Non
            data_off = Noff
            timing_array = Npeaks
        elif set_filter == True and len(data_on) < 3:
            return


        patches =[]

        for pixel in tel_1[data_on]:
            circle = mpatches.Circle((pixel[4],pixel[5]), pixel[6], color='r', fill='black') # create a circle for every pixel on the camera.
            patches.append(circle) # circles are stored in 'patches' array.
            ax.annotate('{0}'.format(int(pixel[1]-1)), xy=(pixel[4], pixel[5] ), 
                        xytext=(pixel[4]-0.00, pixel[5] ),
                        size = 8 , style = 'oblique' , horizontalalignment = 'center')

        collection = PatchCollection(patches, alpha=0.7 , facecolor = 'blue',cmap=matplotlib.cm.inferno )

        ax.add_collection(collection)
        ax.autoscale_view()




        patches =[]



        for pixel in tel_1[data_off]:
            circle = mpatches.Circle((pixel[4],pixel[5]), pixel[6], color='r', fill='black') # create a circle for every pixel on the camera.
            patches.append(circle) # circles are stored in 'patches' array.
            ax.annotate('{0}'.format(int(pixel[1]-1)), xy=(pixel[4], pixel[5] ), 
                        xytext=(pixel[4]-0.00, pixel[5] ),
                        size = 8 , style = 'oblique' , horizontalalignment = 'center')

        collection2 = PatchCollection(patches, alpha=0.7, edgecolor = 'black', facecolor = 'white',  )

        ax.add_collection(collection2)
        ax.autoscale_view()


        ax.set_xlabel('deg')
        ax.set_ylabel('deg')

        #creata an array with the timing of the transients. It is zeroed at the earliest transient.

        collection.set_array(timing_array)
        cbar = fig.colorbar(collection, ax=ax , )
        #cmap=matplotlib.cm.jet
        cbar.set_label('t(s)')
        ax.autoscale_view()

        if not os.path.exists('transientmaps/{}/{}'.format(date_run_on, run_num)):
            os.makedirs('transientmaps/{}/{}'.format(date_run_on, run_num) )
        if set_filter == False:
                ax.set_title('T{} Timing of Transient Peaks '.format(tel_num))
                plt.savefig('transientmaps/{}/{}/{}_T{}_Transient_Time_Map.png'.format(date_run_on,run_num,run_num,tel_num), dpi = 100)
        if set_filter == True:
                ax.set_title('T{} Timing of Transient Peaks (w. time cut) '.format(tel_num))
                plt.savefig('transientmaps/{}/{}/{}_T{}_TIMECUT_Transient_Time_Map.png'.format(date_run_on,run_num,run_num,tel_num), dpi = 100)

        plt.clf()
	plt.close()
	print('''Camera map w. transients saved to 'transientmaps' folder ''')
    else:
        print('No transients to show.')


#############################################################
#   #                                                     # #
#  #            Code Starts Running Here                 #  #
# #                                                     #   #
#############################################################

date_run_on = 'd'+datetime.today().strftime('%d%m%y')

print('did we make it this far?')
# Make a directory to indicate the date that the user ran this code.
if not os.path.exists('TransientRuns/%s' % date_run_on):
    os.makedirs('TransientRuns/%s' % date_run_on)

if not os.path.exists('StarRuns/%s' % date_run_on):
    os.makedirs('StarRuns/%s' % date_run_on)

print(date_run_on)

for i in range(len(PicTempNames)):
    FileName = PicTempNames[i][0]
    runID = FileName[-21:-16]
    tel_num = FileName[-14]
    print 'now reading {}'.format(FileName)

    time_start = time()
    time1 = time()
    print('runID =  {} \n'.format(runID))
    print('telID =  {} \n'.format(tel_num))
    I,V,LI,VI = GetCurrentVoltage(runID)


    median_charge1 = np.array( pd.read_hdf(   FileName   ,'fixed'))
    #median_charge1, num_events1 = pva.TakeMedian(charge1)
    filter_var1, filter_err_var1 , number_of_bins1 , bin_size1 = pva.ped_var_filter(median_charge1,bin_size = 10)
    time2 = time()
    print('Time to load data = {} seconds '.format(time2-time1))
    print('Telescope data has been loaded, median filtered, binned, and average filtered. Starting the star finder. \n')
    time1 = time()
    print('Stars have been searched for, all pixels with stars will NOT be searched for transients!')
    channels_with_stars, rundate1 = StarFinder(runID)
    all_channels = [ii for ii in range(0,499)]
    channels_list1 = list(set(all_channels)  -set(channels_with_stars))
    #SIMBADTransientHunter(filter_var1,filter_err_var1, charge1, channels_list1,  window_size = 200 ,  telID = tel_num)
    Tel_on,Tel_off,Tel_peaks = SIMBADTransientHunter2(filter_var1,filter_err_var1, median_charge1 , channels_list1, rundate1, window_size = 200 ,  telID = tel_num,bin_size = bin_size1, runID = runID)
    del   median_charge1, filter_var1, filter_err_var1, number_of_bins1 , bin_size1 ,
    time_end = time()
    print( '''Run Finished. Time Passed =  {} seconds \n \n Starting New Run \n'''.format(time_end-time_start))

print('''Code has reached its end. Goodbye.''')
