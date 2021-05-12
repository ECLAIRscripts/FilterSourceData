import numpy as np
import xarray as xr
import multiprocessing as mp
from multiprocessing import Pool
from netCDF4 import Dataset
from numba import autojit, prange
import os
import matplotlib
from matplotlib import pyplot as plt
from metpy.calc import pressure_to_height_std
from metpy.units import units
import LES2emu

def calc_LTS(tpot,levels):
	index_750 = min(range(len(levels)), key=lambda i: abs(levels[i]-750))
	index_1000 = min(range(len(levels)), key=lambda i: abs(levels[i]-1000))
	return tpot[int(index_750)-1]-tpot[int(index_1000)-1]

def calcInversion(a):
	return np.max(a)-np.min(a)

def checkProfiles(tpot_profile_lq,levels,LWP_cloudtop,LWP,iwp,iwp_total):
	lts=calc_LTS(tpot_profile_lq,levels)
	#q_inv = np.max(total_water[cloud_top_ind:-1])-np.min(total_water[cloud_top_ind:-1])
	#if((lts > 18.55)): # and (lts < 15)):
	if(True):
		if((LWP_cloudtop > (LWP+iwp_total)*0.5) and (LWP_cloudtop > 0.1*iwp)):
		#	if((1.0  < q_inv)): # and (q_inv < total_water[cloud_top_ind])):
			return True
	return False

def calcLWP(q,levels):
	#clacl pressure differences = dp
	dp = levels[1::] - levels[0:-1] 
	g_rcp = 1.0/9.80665
	pdpg = dp*g_rcp
	lwp = 0
	#lwp2 = np.trapz(q,levels )*g_rcp
	for i in range(0,len(q)):
		lwp = lwp + q[i]*g_rcp*(levels[i+1]-levels[i])

	#print(lwp,len(levels),len(q))
	#input('press.')
	return lwp
	#input()



def checkContit(profile,cloud_top_ind):
	tmp = profile[cloud_top_ind::] #from cloud top to second last layer (second closes to ground)
	idx_clouds = np.where(tmp >= 0.01)[0] #index where there is clouds
	idx_noclouds = np.where(tmp < 0.01)[0] #no clouds from bottom from first cloud to top.
	
	if((idx_noclouds.size > 0) and (idx_clouds.size > 0.0)): #check that there is dry layer above
		try:
			cloud_bottom = np.max(idx_clouds)
			cloud_top = idx_noclouds[np.max(np.where(idx_noclouds <= cloud_bottom))]+1
			return cloud_top+cloud_top_ind,cloud_bottom+cloud_top_ind
		except:
			return 0,0
		#print('cloud top ',cloud_top,'       cloud bottom   ',cloud_bottom)

	return 0,0

def setfont(ax):
	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(28)
	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(28)

def plotProfiles(tpot_profile_lq,total_water,pbl,lat,lons,times,cloud_top_ind,month,levels,ind700,clw,bottom,inv_ind,cdnc,q_pbl,q_inv,lwp,ind):
	lat=np.round(lat,2)
	lons=np.round(lons,2)
	title = 'LAT: '+str(lat)+' LON: '+str(lons)+' TIME: '+str(times)
	#make directory if does not exist
	dir1=str(lat)+'_'+str(lons)
	if not os.path.exists('/lustre/tmp/nordling/eclair/EXPR/LOW_LEVEL_CONTINUES_NOFOG_NOICE_NOABOVECLOUDS/figures/'+ind):
    		os.makedirs('/lustre/tmp/nordling/eclair/EXPR/LOW_LEVEL_CONTINUES_NOFOG_NOICE_NOABOVECLOUDS/figures/'+ind)
	if not os.path.exists('/lustre/tmp/nordling/eclair/EXPR/LOW_LEVEL_CONTINUES_NOFOG_NOICE_NOABOVECLOUDS/figures/'+ind+'/'+dir1):
    		os.makedirs('/lustre/tmp/nordling/eclair/EXPR/LOW_LEVEL_CONTINUES_NOFOG_NOICE_NOABOVECLOUDS/figures/'+ind+'/'+dir1)
	if not os.path.exists('/lustre/tmp/nordling/eclair/EXPR/LOW_LEVEL_CONTINUES_NOFOG_NOICE_NOABOVECLOUDS/figures/'+ind+'/'+dir1+'/'+str(month)):
    		os.makedirs('/lustre/tmp/nordling/eclair/EXPR/LOW_LEVEL_CONTINUES_NOFOG_NOICE_NOABOVECLOUDS/figures/'+ind+'/'+dir1+'/'+str(month))

	save_dir='/lustre/tmp/nordling/eclair/EXPR/LOW_LEVEL_CONTINUES_NOFOG_NOICE_NOABOVECLOUDS/figures/'+ind+'/'+dir1+'/'+str(month)
	#print('LEVLES',levels[cloud_top_ind],cloud_top_ind)

	font = {'family': 'serif',
		'color':  'darkred',
		'weight': 'normal',
		'size': 16,
		}
	f, (ax1, ax2) = plt.subplots(1,2,figsize=(25,25))
	ax3 = ax2.twiny()
	#ax4 = ax2.twiny()
	f.suptitle(title) 
	ax1.plot(tpot_profile_lq,levels,linewidth=3.0,label='tpot')
	#ax1.text(2, 0.65,'q_pbl: '+str(round(q_pbl,4))+' q_inv: '+str(round(q_inv,4)), fontdict=font)
	ax1.plot([np.min(tpot_profile_lq[cloud_top_ind::])-5,np.max(tpot_profile_lq[cloud_top_ind::])+5],[levels[cloud_top_ind],levels[cloud_top_ind]],linewidth=3.0,label='cloud top '+str(round(q_pbl,4)))
	ax1.plot([np.min(tpot_profile_lq[cloud_top_ind::])-5,np.max(tpot_profile_lq[cloud_top_ind::])+5],[levels[inv_ind],levels[inv_ind]],linewidth=3.0,label='inversion top '+str(round(q_inv,4)))
	ax1.plot([np.min(tpot_profile_lq[cloud_top_ind::])-5,np.max(tpot_profile_lq[cloud_top_ind::])+5],[levels[ind700],levels[ind700]],linewidth=3.0,label='700hpa level '+str(round(lwp,4)))
	ax1.plot([np.min(tpot_profile_lq[cloud_top_ind::])-5,np.max(tpot_profile_lq[cloud_top_ind::])+5],[levels[bottom],levels[bottom]],linewidth=3.0,label='cloud bottom')
	ax1.legend()
	ax2.plot(total_water,levels,linewidth=3.0,label='total water')
	ax3.plot(clw,levels,'r',linewidth=3.0,label='clw')
	#ax4.plot(cdnc,levels,'r',linewidth=3.0,label='cdnc')
	ax2.plot([np.min(total_water[cloud_top_ind::]),np.max(total_water[cloud_top_ind::])],[levels[cloud_top_ind],levels[cloud_top_ind]],linewidth=3.0,label='cloud top')
	ax2.legend()
	ax3.xaxis.set_ticks_position("top")
	ax3.xaxis.set_label_position("top")
	ax3.spines["bottom"].set_position(("axes", -0.15))
	matplotlib.rcParams.update({'font.size': 28})
	plt.rc('axes', labelsize=28)
	ax1.invert_yaxis()
	ax2.invert_yaxis()
	#ax2.set_xlim([-0.2,1.2])
	ax1.set_xlim([np.min(tpot_profile_lq[ind700::])-2,np.max(tpot_profile_lq[ind700::])+2])
	ax1.set_ylim([levels[-1],650])
	ax2.set_ylim([levels[-1],650])
	ax1.set_title('Potential temperature')
	ax2.set_title('Total water '+str(round(q_pbl,4))+' '+str(round(q_inv,4)))
	#ax3.set_title('Cloud water')
	ax1.set_ylabel('Pressure (hPa)')
	ax1.set_xlabel('K')
	ax2.set_xlabel('g/kg)')
	setfont(ax1)
	setfont(ax2)
	f.savefig(save_dir+'/'+str(times)+'.png')
	plt.close(f)
	np.savetxt(save_dir+'/'+str(times)+'.txt',clw)
	

def analyise_month(month):
	print('READING DATA FROM MONTH........',str(month))
	m=str(month)
	if(month < 10):
		m = '0'+m

	ifile = xr.open_dataset('../../../emulator_on/emulator_on_2001'+m+'.01_echam.nc')
	tpot = np.copy(ifile['tpot'].to_masked_array()[:,:,:,:])
	cloud_water = np.copy(ifile['xl'].to_masked_array()[:,:,:,:])
	q = np.copy(ifile['q'].to_masked_array()[:,:,:,:])
	aps = np.copy(ifile['aps'].to_masked_array()[:,:,:])
	iwp_= np.copy(ifile['xivi'].to_masked_array()[:,:,:])
	ice= np.copy(ifile['xi'].to_masked_array()[:,:,:])
	lats = np.copy(ifile['lat'].to_masked_array()[:])
	lons = np.copy(ifile['lon'].to_masked_array()[:])
	times_ =np.copy(ifile['time'].to_masked_array()[:])
	a = np.copy(ifile['hyam'][:].to_masked_array())
	b = np.copy(ifile['hybm'][:].to_masked_array())
	a_ = np.copy(ifile['hyai'][:].to_masked_array())
	b_ = np.copy(ifile['hybi'][:].to_masked_array())
	accretion1_int = np.copy(ifile['accretion1_int'][:,:,:].to_masked_array())
	accretion2_int= np.copy(ifile['accretion2_int'][:,:,:].to_masked_array())
	autoconv_int= np.copy(ifile['autoconv_int'][:,:,:].to_masked_array())
	rainevap_int= np.copy(ifile['rainevap_int'][:,:,:].to_masked_array())
	cos_mu = np.copy(ifile['cos_mu'][:,:,:].to_masked_array()) #CHANGE BACK"!!!!!!!!!!! to COSMU
	sea_ice = np.copy(ifile['seaice'].to_masked_array()[:,:,:])
	land = np.copy(ifile['slm'].to_masked_array()[:,:,:])
	ifile.close()
	
	ifile = xr.open_dataset('../../../emulator_on/emulator_on_2001'+m+'.01_tracer.nc')
	cdnc = np.copy(ifile['CDNC'][:,:,:,:].to_masked_array())
	particles = np.copy(ifile['NUM_AS'][:,:,:,:].to_masked_array())+np.copy(ifile['NUM_CS'][:,:,:,:].to_masked_array())+np.copy(ifile['NUM_KS'][:,:,:,:].to_masked_array())
	num_as = np.copy(ifile['NUM_AS'][:,:,:,:].to_masked_array())
	num_cs = np.copy(ifile['NUM_CS'][:,:,:,:].to_masked_array())
	num_ks = np.copy(ifile['NUM_KS'][:,:,:,:].to_masked_array())
	ifile.close()

	ifile = xr.open_dataset('../../../emulator_on/emulator_on_2001'+m+'.01_ham.nc')
	rdry_AS = np.copy(ifile['rdry_AS'][:,:,:].to_masked_array())
	B_AS = np.copy(ifile['B_AS'][:,:,:].to_masked_array())
	ifile.close()

	ifile = xr.open_dataset('../../../emulator_on/emulator_on_2001'+m+'.01_activ.nc')
	W= np.copy(ifile['W'][:,:,:].to_masked_array())
	ifile.close()	

	strato = np.zeros((tpot.shape[2],tpot.shape[3]))
	tpot_inv = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	q_inv = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	tpot_pbl = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	pbl = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	pbl_m = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	lwp = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	lwp_1 = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	lwp_2 = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	lts_ = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	cdnc_ = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	q_mean = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	aero = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	num_as_ = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	num_cs_ = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	num_ks_ = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	mask = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))

	q_pbl = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	q_pbl1 = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))	
	q_pbl2 = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))	
	cos_mu_= np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	clw_max= np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	cloud_top_= np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	cloud_base_= np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	B_AS_ = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	rdry_AS_ = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	rdry_AS_eff = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))

	accretion1_int_ = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	accretion2_int_ = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	autoconv_int_ = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	rainevap_int_ = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))
	W_ = np.zeros((tpot.shape[0],tpot.shape[2],tpot.shape[3]))

	#print(strato.shape)
	#levels=ifile['plev'].to_masked_array()[1:9]/100
	prof_count=0
	prof_lim=2000	
	print('START ANALYSIS FOR..........',month)
	profile_count=0
	for i in range(0,tpot.shape[0]):
		for x in range(0,96):
			for y in range(0,192):
				if((sea_ice[i,x,y] > 0) and (land[i,x,y] > 0)):
					continue
				#rh_profile = rh[i,:,x,y]
				#cloud_profie= clouds[i,:,x,y]
				cloud_water_profile = cloud_water[i,:,x,y]*1000
				ice_profile = ice[i,:,x,y]*1000
				tpot_profile=tpot[i,:,x,y]#-(2.501e6/1005.7)*(cloud_water_profile/1000)
				#tpot_profile=tpot[i,:,x,y]
				IWP_total = iwp_[i,x,y]*1000
				total_water = cloud_water_profile + q[i,:,x,y]*1000
				levels = (a+b*aps[i,x,y])/100
				levels_ = (a_+b_*aps[i,x,y])/100
				ind700 = min(range(len(levels)), key=lambda i: abs(levels[i]-700))
				LWP=calcLWP(cloud_water_profile,levels_*100)
	
				#check if we have low level clouds below 700hPa and dry layer above that
				#input()
				#if(True):
				#print(np.where(cloud_water_profile[ind700::] > 0.01))
				if(((np.where(cloud_water_profile[ind700::] > 0.01)[0]).size >= 1) and (cloud_water_profile[len(cloud_water_profile)-1] < 0.01)):# and (ice_profile < 0.0001)):
					#print("find clouds")
					#if((checkContit(cloud_water_profile,ind700))):
					if(True):
						cloud_top_ind,cloud_bottom_ind = checkContit(cloud_water_profile,ind700)
						#print('Cloud top:',cloud_top_ind,'Cloudbase:',cloud_bottom_ind)
						if((cloud_top_ind==0) or (cloud_bottom_ind == 0)):
							continue
						#print('find continues cloud')
						#cloud top index is max_cloud water under 700hpa level
						#cloud_top_ind = np.min(np.where(cloud_water_profile[ind700::] > 0.0001))+ind700
						#cloud_bottom_ind = np.max(np.where(cloud_water_profile[ind700::] > 0.0001))+ind700
						LWP_cloudtop = calcLWP(cloud_water_profile[cloud_top_ind::],levels_[cloud_top_ind::]*100)	#calc lwp under low-level cloud
						iwp = calcLWP(ice_profile[cloud_top_ind::],levels_[cloud_top_ind::]*100)
						#check 
						#if(checkProfiles(tpot_profile,levels,LWP_cloudtop,LWP,iwp,IWP_total) and (cloud_top_ind > ind700)):
						#print('Cloud top:',cloud_top_ind,'Cloudbase:',cloud_bottom_ind)
						#print('LWP_cloudtop:',LWP_cloudtop)
						#print('(LWP+IWP_total)*0.5',(LWP+IWP_total)*0.5)
						#print('LWP_cloudtop',LWP_cloudtop)
						#print('0.1*iwp',0.1*iwp)
						#plotProfiles(tpot_profile,total_water,0,lats[x],lons[y],times_[i],cloud_top_ind,month,levels,ind700,cloud_water_profile,cloud_bottom_ind,0,0,0,0,0,'neg')	
						if((LWP_cloudtop > (LWP+IWP_total)*0.5) and (LWP_cloudtop > 0.1*iwp)):							
							strato[x,y] =strato[x,y]+ 1
							mask[i,x,y] = 1
							#print('FIND STATOCUMULS CLOUD',month,strato[x,y])
							inv_ind = cloud_top_ind-2
							inv_low = np.min((cloud_bottom_ind+2,46))
							#print('indwex:',inv_ind,inv_low,ind700)
							tpot_inv[i,x,y] = calcInversion(tpot_profile[inv_ind:inv_low])
							pbl[i,x,y] = (aps[i,x,y] - (levels_[cloud_top_ind]*100))/100.0
							pbl_m[i,x,y] = -1*(pressure_to_height_std(aps[i,x,y]*units.Pa) - pressure_to_height_std((levels_[cloud_top_ind]*100)*units.Pa)).magnitude*1000
							tpot_pbl[i,x,y] = np.min(tpot_profile[inv_ind:inv_low])
							#lwp_1[i,x,y] = LWP_cloudtop*np.mean(cloud_profie[cloud_top_ind::])
							#lwp_2[i,x,y] = LWP_cloudtop*(clw_max[i,x,y]/np.mean(total_water[inv_ind+2::]))
							lwp[i,x,y] = LWP_cloudtop
							q_mean[i,x,y] = np.mean(total_water[inv_ind+2::])

							q_inv[i,x,y] = calcInversion(total_water[inv_ind:inv_low])
							cloud_top_[i,x,y] = pressure_to_height_std((levels_[cloud_top_ind]*100)*units.Pa).magnitude*1000
							cloud_base_[i,x,y] = pressure_to_height_std((levels_[cloud_bottom_ind ]*100)*units.Pa).magnitude*1000
							q_pbl[i,x,y] = LES2emu.solve_rw_lwp(aps[i,x,y], tpot_pbl[i,x,y], lwp[i,x,y]*0.001, pbl_m[i,x,y]*0.001 )*1000.
							#q_pbl1[i,x,y] = LES2emu.solve_rw_lwp(aps[i,x,y], tpot_pbl[i,x,y], lwp_1[i,x,y]*0.001, pbl_m[i,x,y]*0.001 )*1000.
							#q_pbl2[i,x,y] = LES2emu.solve_rw_lwp(aps[i,x,y], tpot_pbl[i,x,y], lwp_2[i,x,y]*0.001, pbl_m[i,x,y]*0.001 )*1000.
							aero[i,x,y] = particles[i,46,x,y]
							clw_max[i,x,y] = np.max(cloud_water_profile[cloud_top_ind::])



							if(cloud_top_ind != cloud_bottom_ind):
								cdnc_[i,x,y] = np.nanmax(cdnc[i,(cloud_top_ind)-1:(cloud_bottom_ind)+1	,x,y])*1e-6
							else:
								cdnc_[i,x,y] = np.nanmax(cdnc[i,cloud_top_ind,x,y])*1e-6
							num_as_[i,x,y] = num_as[i,46,x,y]
							num_cs_[i,x,y] = num_cs[i,46,x,y]
							num_ks_[i,x,y] = num_ks[i,46,x,y]
							cos_mu_[i,x,y] = cos_mu[i,x,y]
							B_AS_[i,x,y] = B_AS[i,46,x,y]
							rdry_AS_[i,x,y] = rdry_AS[i,46,x,y]
							rdry_AS_eff[i,x,y] = (B_AS_[i,x,y]/1.008)**(1.0/3.0)*rdry_AS_[i,x,y]
							accretion1_int_[i,x,y] = accretion1_int[i,x,y]
							accretion2_int_[i,x,y] = accretion2_int[i,x,y]
							autoconv_int_[i,x,y] = autoconv_int[i,x,y]
							rainevap_int_[i,x,y] = rainevap_int[i,x,y]
							W_[i,x,y] = W[i,46,x,y]
							'''
							#if(prof_lim < prof_count):
							if(((q_pbl[i,x,y] - q_inv[i,x,y]) > 0) and (q_pbl[i,x,y] > 0)):
								#if(cdnc_[i,x,y] < 10):
								#plot profile
									plotProfiles(tpot_profile,total_water,pbl,lats[x],lons[y],times_[i],cloud_top_ind,month,levels,ind700,cloud_water_profile,cloud_bottom_ind,inv_ind,cdnc[i,:,x,y],q_pbl[i,x,y],q_inv[i,x,y],lwp[i,x,y],'pos')

							if(((q_pbl[i,x,y] - q_inv[i,x,y]) < 0) and (q_pbl[i,x,y] > 0)):
								#if(cdnc_[i,x,y] < 10):
									#plot profile
								plotProfiles(tpot_profile,total_water,pbl,lats[x],lons[y],times_[i],cloud_top_ind,month,levels,ind700,cloud_water_profile,cloud_bottom_ind,inv_ind,cdnc[i,:,x,y],q_pbl[i,x,y],q_inv[i,x,y],lwp[i,x,y],'neg')		
								prof_count =0
							prof_count+=1
							'''
	
	print('SAVE FILE',month) 
	rootgrp = Dataset("/lustre/tmp/nordling/eclair/EXPR//EMULATOR_INPUT/LOW_LEVEL_cold7"+str(month)+".nc", "w", format="NETCDF4")
	time = rootgrp.createDimension("time", tpot.shape[0])
	lat = rootgrp.createDimension("lat", tpot.shape[2])
	lon = rootgrp.createDimension("lon", tpot.shape[3])
	times = rootgrp.createVariable("time","f8",("time",))
	latitudes = rootgrp.createVariable("lat","f8",("lat",))
	longitudes = rootgrp.createVariable("lon","f8",("lon",))
	strato_ = rootgrp.createVariable("strato","f8",("lat","lon",))

	tpot_inv_ = rootgrp.createVariable("tpot_inv","f8",("time","lat","lon",))
	q_inv_ = rootgrp.createVariable("q_inv","f8",("time","lat","lon",))
	pbl_ = rootgrp.createVariable("pbl","f8",("time","lat","lon",))
	pbl_m_ = rootgrp.createVariable("pbl_m","f8",("time","lat","lon",))
	lwp_ = rootgrp.createVariable("lwp","f8",("time","lat","lon",))
	lwp__1 = rootgrp.createVariable("lwp1","f8",("time","lat","lon",))
	lwp__2 = rootgrp.createVariable("lwp2","f8",("time","lat","lon",))
	tpot_pbl_ = rootgrp.createVariable("tpot_pbl","f8",("time","lat","lon",))
	aerosols = rootgrp.createVariable("particles","f8",("time","lat","lon",))
	cdnc__ = rootgrp.createVariable("cdnc","f8",("time","lat","lon",))
	cos_mu__= rootgrp.createVariable("cos_mu","f8",("time","lat","lon",))
	mask_= rootgrp.createVariable("mask","f8",("time","lat","lon",))
	clw_max_= rootgrp.createVariable("clw_max","f8",("time","lat","lon",))
	q_pbl_ = rootgrp.createVariable("q_pbl","f8",("time","lat","lon",))
	q_pbl_1 = rootgrp.createVariable("q_pbl1","f8",("time","lat","lon",))
	q_pbl_2 = rootgrp.createVariable("q_pbl2","f8",("time","lat","lon",))

	B_AS__ = rootgrp.createVariable("B_AS","f8",("time","lat","lon",))
	rdry_AS__ = rootgrp.createVariable("rdry_AS","f8",("time","lat","lon",))
	rdry_AS_eff_ = rootgrp.createVariable("rdry_AS_eff","f8",("time","lat","lon",))

	cloud_top__= rootgrp.createVariable("cloud_top","f8",("time","lat","lon",))
	cloud_base__= rootgrp.createVariable("cloud_base","f8",("time","lat","lon",))

	as_ = rootgrp.createVariable("as","f8",("time","lat","lon",))
	cs_ = rootgrp.createVariable("cs","f8",("time","lat","lon",))
	ks_ = rootgrp.createVariable("ks","f8",("time","lat","lon",))

	W__ = rootgrp.createVariable("W","f8",("time","lat","lon",))
	accretion1_int__ = rootgrp.createVariable("accretion1_int","f8",("time","lat","lon",))
	accretion2_int__ = rootgrp.createVariable("accretion2","f8",("time","lat","lon",))
	autoconv_int__ = rootgrp.createVariable("autoconv","f8",("time","lat","lon",))
	rainevap_int__ = rootgrp.createVariable("rainevap","f8",("time","lat","lon",))

	
	times[:] = times_
	latitudes[:] = lats
	longitudes[:] = lons
	strato_[:,:] = strato
	cos_mu__[:,:,:] = cos_mu_ 
	mask_[:,:,:] = mask
	clw_max_[:,:,:] = clw_max
	cloud_top__[:,:,:] = cloud_top_ 
	cloud_base__[:,:,:] = cloud_base_ 

	q_pbl_[:,:,:] = q_pbl
	q_pbl_1[:,:,:] = q_pbl1
	q_pbl_2[:,:,:] = q_pbl2

	B_AS__[:,:,:] = B_AS_
	rdry_AS__[:,:,:] = rdry_AS_
	rdry_AS_eff_[:,:,:] =rdry_AS_eff

	pbl_m_[:,:,:] = pbl_m
	tpot_inv_[:,:,:] = tpot_inv
	q_inv_[:,:,:] = q_inv
	pbl_[:,:,:] = pbl
	lwp_[:,:,:] = lwp
	lwp__1[:,:,:] = lwp_1
	lwp__2[:,:,:] = lwp_2
	tpot_pbl_[:,:,:] = tpot_pbl
	aerosols[:,:,:] = aero
	cdnc__[:,:,:] = cdnc_
	as_[:,:,:] = num_as_
	cs_[:,:,:] = num_cs_
	ks_[:,:,:] = num_ks_
	accretion1_int__[:,:,:] = accretion1_int_
	accretion2_int__[:,:,:] = accretion2_int_
	autoconv_int__[:,:,:] = autoconv_int_
	rainevap_int__[:,:,:] = rainevap_int_
	W__[:,:,:] = W_

	rootgrp.close()
	
def main():
	print('Main function')
	#analyise_month(1)
	m=[1,2,3,4,5,6,7,8,9,10,11,12]
	#m=[12]
	pool = Pool(processes= 12)
	for i in pool.imap_unordered(analyise_month, m):
		print(i)
	#jobs = []
	#for i in m:
	#	p = mp.Process(target=analyise_month, args=(i,))
	#	jobs.append(p)
	#	p.start()

main()
