from astropy.io import fits
from astropy.wcs import WCS
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import os
import sys


from datetime import datetime
startTime = datetime.now()

#do something

print datetime.now() - startTime



print "analyzing "
print sys.argv[2]
#run     python muse_v1.py ADP.2016-06-15T12:54:04.046.fits 1


data_path="/data/highzgal/sourabh/VLT_MUSE/"
primary_path="/data/highzgal/sourabh/VLT_MUSE/analysis/output_files/"


file_shift="%s/region_files_corr/file_shift.txt"%(primary_path)
shift_x=0.0
shift_y=0.0
#shift_x=0.0
#shift_y=0.0

x3=[]
with open(file_shift, "r") as f3:              #
        for line in f3:
		line = line.strip()
		line = line.split()
		x3.append(line)
f3.close()

shift_x=float(x3[int(sys.argv[2])-1][0])
shift_y=float(x3[int(sys.argv[2])-1][0])


file1="%s%s"%(data_path,sys.argv[1]) #data files are here !!!!!!!!!!!!!  be very careful original data

file_box="%sbox_astro_corr_new/reg%s.txt"%(primary_path,sys.argv[2])
file_box1="%sbox_astro_corr_new/reg%s_check.txt"%(primary_path,sys.argv[2])
hdulist = fits.open(file1)


from astropy.wcs import WCS
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
header =fits.getheader(file1,1)

header1 =fits.getheader(file1,0)
w=WCS(header)

NAXIS1 = header["NAXIS1"]
NAXIS2 = header["NAXIS2"]
NAXIS3 = header["NAXIS3"]
x = np.arange(NAXIS1)
y = np.arange(NAXIS2)

X, Y = np.meshgrid(x, y)
coord = w.wcs_pix2world(X,Y,NAXIS3, 0)

#coord is a list ra and dec and some other shit

scidata=hdulist[1].data
#print scidata.shape
#print scidata[0,1,2]
#wcs = WCS(hdulist)
#print wcs



pad=30
box=15
z_beg=0.81
z_end=1.5
spec_rad=5
Fe_lam=2626.45
all_gal="gal_ha_sel_final.txt"
#all_gal="gal_redshift_cat_new.txt"
x1=[]
with open(all_gal, "r") as f:              #
	for line in f:
		line = line.strip()
                line = line.split()
                x1.append(line)
f.close()
num_rows = sum(1 for line in open(all_gal))
#convert index to wavelength 


j=0
gal_ra=[]
gal_dec=[]
gal_z=[]
gal_id=[]
for i in range(0,num_rows):

	if (float(x1[i][3])>z_beg and float(x1[i][3])<z_end):
		gal_id.append(x1[i][0])		
		gal_ra.append(x1[i][1])
		gal_dec.append(x1[i][2])
		gal_z.append(x1[i][3])		
		j=j+1
tot_gal=j

for j in range(0,tot_gal):
	gal_ra[j]=float(gal_ra[j])
	gal_dec[j]=float(gal_dec[j])
	gal_z[j]=float(gal_z[j])



lam_beg=(1+z_beg)*Fe_lam
lam_end=(1+z_end)*Fe_lam

lam_min=header1["WAVELMIN"]*10 ###wavelength are in nm.
lam_max=header1["WAVELMAX"]*10
total=NAXIS3
del_lam=(lam_max-lam_min)/total

ind_beg=int((lam_beg-lam_min)/del_lam)
ind_ned=int((lam_end-lam_min)/del_lam)

ax1=header["NAXIS1"]
ay1=header["NAXIS2"]
final_out=np.zeros((tot_gal,box,box))

k1=0
#enter the file here
#find image centre and extent in ra and dec
#for i in range(ind_beg:ind_end):

k2=0
am=[]
an=[]
am1=[]
an1=[]
#print NAXIS1,NAXIS2
for x in range(0, NAXIS1):
	for y in range(0,NAXIS2):
		if np.isnan(scidata[0,y,x])==False :
			am.append(x)
			an.append(y)	
			#print x,y
			k2=k2+1
			break

#print k2
k2=0
for x in range(0, NAXIS1):
	for y in range(0,NAXIS2):
		if np.isnan(scidata[0,NAXIS2-y-1,NAXIS1-x-1])==False :
			am1.append(NAXIS1-x)
			an1.append(NAXIS2-y)	
			#print NAXIS1-x,NAXIS2-y
			k2=k2+1
			break


#print k2





arg_max=np.argmax(an1)
arg_min=np.argmin(an)

a4=am1[arg_max]
b4=an1[arg_max]
a2=am[arg_min]
b2=an[arg_min]


a1=am[0]
b1=an[0]

a3=am1[0]
b3=an1[0]

#slope1=(b2-b1)/(a2-a1)
#slope2=(b4-b3)/(a4-a3)

'''
a_min=max((a1+a2)*0.5,(a1+a4)*0.5)
a_max=min((a2+a3)*0.5,(a3+a4)*0.5)
b_min=max((b1+b2)*0.5,(b2+b3)*0.5)
b_max=min((b3+b4)*0.5,(b1+b4)*0.5)
'''

#(x1,y1)=(pad,b1)
#(x3,y3)=(a3-pad,b3)
#print a1,b1,a2,b2,a3,b3,a4,b4
#print a_min,a_max,b_min,b_max

#y=y1+slope1*(x-x1)
'''

nx=int(a_max-a_min)
ny=int(b_max-b_min)
'''
import matplotlib.path as mplPath
poly = [a1, b1, a2,b2,a3,b3,a4,b4]
bbPath = mplPath.Path(np.array([[a1, b1],
                     [a2, b2],
		     [a3, b3],
		     [a4, b4]]))

k1=0

pixx_g=[]
pixy_g=[]
pixx_c=[]
pixy_c=[]
gal_id_sel=[]
gal_z_sel=[]
z_index_g=[]
ind_g=[]
im=np.zeros((box,box))

file_check= 'ADP.2016-06-15T12:54:04.047.fits'

for j in range(0,tot_gal):
	
	r1 = w.wcs_world2pix(float(gal_ra[j]),float(gal_dec[j]),NAXIS3,0)
	#print r1
	pixx = (r1[0])# changed integer thing
	pixy = (r1[1])
	pixx_c=pixx+shift_x
	pixy_c=pixy+shift_y
	z_index = int(((1+gal_z[j])*Fe_lam-lam_min)/del_lam)
	
	#lam[k]=(lam_min+del_lam*k)/(1+float(sys.argv[6]))

	ch1=bbPath.contains_point((pixx-box*0.5, pixy-box*0.5))
	ch2=bbPath.contains_point((pixx+box*0.5, pixy-box*0.5))
	ch3=bbPath.contains_point((pixx+box*0.5, pixy+box*0.5))
	ch4=bbPath.contains_point((pixx-box*0.5, pixy+box*0.5))

	if(ch1==1 and ch2==1 and ch3==1 and ch4==1 and z_index<NAXIS3):
		ind_g.append(j)
		pixx_g.append(pixx_c)
		pixy_g.append(pixy_c)
		gal_id_sel.append(float(gal_id[j]))
		gal_z_sel.append(float(gal_z[j]))
		print gal_id[j]
		z_index_g.append(z_index)		
		k1=k1+1
		r2=w.wcs_pix2world(pixx_c, pixy_c,NAXIS3,0)
		gal_ra[j]=r2[0]
		gal_dec[j]=r2[1]
		#print gal_ra[j],gal_dec[j]
		print k1
		#print k1, gal_id[j], gal_ra[j], gal_dec[j], gal_z[j], z_index, pixx, pixy, file1
		#print k1, gal_id[j], gal_ra[j], gal_dec[j], gal_z[j], sys.argv[1]
		
		for l in range(0,box):
			for m in range(0,box):
				final_out[j,m,l]=final_out[j,m,l]+scidata[z_index,pixy_c-box*0.5+m,pixx_c-box*0.5+l]
				im[m,l]=final_out[j,m,l] 
			
		#		im[m,l]=scidata[z_index,pixy-box*0.5+m,pixx-box*0.5+l]
		#create fits for each galaxy and header file with info

		#print k1, j, gal_id[j], gal_ra[j], gal_dec[j], gal_z[j], z_index, pixx, pixy
		#print z_index
		prihdr = fits.getheader(file1,1)
		prihdr['CRPIX1'] = box/2
		prihdr['CRPIX2'] = box/2
		prihdr['CRVAL1'] = float(gal_ra[j])
		prihdr['CRVAL2'] = float(gal_dec[j])
		filenew="%sgalaxy_stamps_corr/reg%s/gal%s_%s.fits"%(primary_path,sys.argv[2],k1,gal_id[j])
		file_path="%sgalaxy_stamps_corr/reg%s/"%(primary_path,sys.argv[2])
		file_path_spec="%sgalaxy_stamps_corr/reg%s_spec/"%(primary_path,sys.argv[2])

		exist_dir_check=os.path.exists(file_path)
		if exist_dir_check==False:
			os.makedirs(os.path.dirname(file_path))
		exist_dir_check=os.path.exists(file_path_spec)
		if exist_dir_check==False:
			os.makedirs(os.path.dirname(file_path_spec))

		fits.writeto(filenew, data=im, header=header,clobber=True)
		im=np.zeros((box,box))
final_out1=np.zeros((NAXIS3))
lam=np.zeros((NAXIS3))
flux_sum_circ=np.zeros((NAXIS3))
j=0
a1=[]
b1=[]

for l in range(-spec_rad,spec_rad):
	for m in range(-spec_rad,spec_rad):
		if (np.sqrt(l**2+m**2))<=spec_rad:
			a1.append(l)
			b1.append(m)
			j=j+1

print (a1[0],b1[0])
circ_points=j
for i in range(k1):
	for j in range(0,NAXIS3):
		for l in range(0,circ_points):
			om=pixy_g[i]+b1[l]
			on=pixx_g[i]+a1[l]
			flux_sum_circ[j]=flux_sum_circ[j]+scidata[j,om,on]

	

for i in range(k1):
	print "spectrum calculaiton of galaxy"
	print i
	for j in range(NAXIS3):
		for l in range(0,box):
			for m in range(0, box):
				final_out1[j]=final_out1[j]+scidata[j,pixy_g[i]-box*0.5+m, pixx_g[i]-box*0.5+l]
#				if (np.sqrt(l**2+m**2))<=spec_rad:
 #                                       flux_sum_circ[j]=flux_sum_circ[j]+scidata[j,pixy_g[i]+m,pixx_g[i]+l]
#box_val=np.arange(0,box,1.0)
		lam[j]=(lam_min+del_lam*j)/(1+gal_z_sel[i])
	np.savetxt("%s/galaxy_stamps_corr/reg%s_spec/reg%s_gal%s_id%s_spec.txt"%(primary_path,sys.argv[2],sys.argv[2],(i+1),gal_id_sel[i]),\
	np.vstack((lam,final_out1,flux_sum_circ)).T)
        final_out1=np.zeros((NAXIS3))
	lam=np.zeros((NAXIS3))


#		lam[j]=(lam_min+del_lam*j)/(1+gal_z_sel[i])
		




#np.savetxt("new_region.reg", np.vstack((hstphot_goods['ra'],hstphot_goods['dec'])).T, fmt="fk5; circle %.6fd  %.6fd 1.5\" # width=1",header='fk5')#, comment='')
#file_check= 'ADP.2016-06-15T12:54:04.047.fits'
#np.savetxt("new_region_no_shift.reg", np.vstack((pixx_g,pixy_g)).T, fmt="circle( %.6f , %.6f , 7.5)")#, comment='')
np.savetxt("%s/region_files_corr/region_box_reg%s.reg"%(primary_path,sys.argv[2]),  np.vstack((pixx_g,pixy_g,gal_id_sel)).T, fmt="box( %.6f , %.6f , 30.0,30.0,0)# text={%.1f}")
np.savetxt("%s/region_files_corr/region_circ_reg%s.reg"%(primary_path,sys.argv[2]), np.vstack((pixx_g,pixy_g, gal_id_sel)).T, fmt="circle( %.6f , %.6f , 7.5)# text={%.1f}")#, comment='') ###to be used for furtehr analysis






#for l in range(0,box):
#	for m in range(0,box):
			                                #final_out[j,m,l]=final_out[j,m,l]+scidata[z_index,pixy-box*0.5+m,pixx-box*0.5+l]
							                                #im[m,l]=final_out[
#		im[m,l]=scidata[619,pixy-box*0.5+m,pixx-box*0.5+l]

#fits.writeto("check.fits", data=im, header=header, clobber=True)

#ind_g=int(ind_g)

with open(file_box, "w") as out_file:
	for i in range(0,k1):
		d=ind_g[i]
		out_string= ""
		out_string+= str(i+1)
		out_string += " " + str(gal_id[d])
		out_string += " " + str(gal_ra[d])
		out_string += " " + str(gal_dec[d])
		out_string += " " + str(gal_z[d])
		out_string += " " + str(z_index_g[i])
		out_string += " " + str(pixx_g[i])
		out_string += " " + str(pixy_g[i])
		out_string += " " + str(sys.argv[1])
		out_string += "\n"
		out_file.write(out_string)

out_file.close()




#with open(file_box, "w") as out_file:
 #       for i in range(0,k1):
#		d=ind_g[i]
#		out_string+= ""
#		out_string+= str(i+1)
#		out_string += " " + str(gal_id[d])





with open(file_box1, "w") as out_file1:
	for i in range(0,k1):
		d=ind_g[i]
		out_string1= ""
		out_string1+= str(i+1)
		#out_string1 += " " + str(ind_g[i])
		out_string += " " + str(gal_id[d])
		out_string1 += " " + str(gal_ra[d])
		out_string1 += " " + str(gal_dec[d])
		out_string1 += " " + str(gal_z[d])
		out_string1 += " " + str(sys.argv[1])
		out_string1 += "\n"
		out_file1.write(out_string)

out_file1.close()


#do something
print "time taken"
print datetime.now() - startTime
print "start "
print startTime


#im=np.zeros((box,box))
#############################################################im=final_out.sum(axis=0)
#for l in range(0,box):
#	for m in range(0,box):
#		im[m,l]=final_out[4129,m,l]
#fig,ax = plt.subplots(1,1,figsize=(5,5),tight_layout=True)
#show=ax.imshow(im,origin='lower',interpolation='none',norm=LogNorm(1,200))
#plt.colorbar(show)
#plt.title('trial')
#fig.savefig('trial.png', dvi=400)
#plt.show()
#import pandas as pd
#df = pd.DataFrame(scidata)
#print pd.isnull(scidata)
#print np.isnan(scidata)#scidata.isnull().sum().sum()
#x = np.arange(box) 
#y = np.arange(box)
#plt.ploti(final_out())
