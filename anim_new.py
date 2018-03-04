import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import animation as anim
from astropy.io import fits
from astropy.io import fits
from scipy.interpolate import interp1d
from matplotlib.colors import LogNorm
from matplotlib import colors
import matplotlib.patches as patches

from mpl_toolkits.axes_grid1 import make_axes_locatable






file1= "ADP.2016-06-15T12:54:04.046.fits"
file2= "ADP.2016-06-15T12:54:04.047.fits" 
file3="/home/sourabh/muse_analysis/reg1/gal_17328.0.fits"
gal_dat="gal1.txt"
spec_dat="gal1_reg1_id17328.0_spec.txt"

hdulist1=fits.open(file1)
data3d=fits.open(file1)

hdulist2=fits.open(file2)
data2d=hdulist2[1].data

hdulist3=fits.open(file3)
gal1_2d=hdulist3[0].data

x1=[]

with open(gal_dat, "r") as f:              #
    for line in f :
            line = line.strip()
            line = line.split()
            x1.append(line)
f.close()
num_rows = sum(1 for line in open(gal_dat))

fig = plt.figure(tight_layout=True)
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 4)

#line, = ax1.plot([], [])

#def init():
 #   line.set_data([], [])
    
  #  return line

Fe_lam1=2586.65
Fe_lam2=2600.16
Fe_lam1_st=2612.65
Fe_lam2_st=2626.45

box=15
def animate(i):
    #plt.clf()
    
    im1=ax1.imshow(data2d,origin='lower',interpolation='none',cmap='Greys' ,norm=LogNorm(0.017,1.0))
    ax1.set_title("region1")
    ax1.add_patch(
        patches.Rectangle(
            (float(x1[i][6])-box/2, float(x1[i][7])-box/2),   # (x,y)
            box,          # width
            box,          # height
            fill=False,
            edgecolor="red"
        )
    )
    

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im1, cax=cax, orientation='vertical')
    file_2d="/home/sourabh/muse_analysis/reg1/gal_%s.fits"%(x1[i][1])
    hd=fits.open(file_2d)
    gal1_2d=hd[0].data
    
    ax2.imshow(gal1_2d,origin='lower',interpolation='none',cmap='Greys' ,norm=LogNorm(0.017,0.5))
    ax2.set_title("gal%s id %s"%(i+1,x1[i][1]))
    
    
    spec_dat="/home/sourabh/muse_analysis/galaxy_stamps_new/reg1_spec/gal%s_reg1_id%s_spec.txt"%(i+1,x1[i][1])
    x2=[]
    with open(spec_dat, "r") as f1:              #
        for line in f1 :
            line = line.strip()
            line = line.split()
            x2.append(line)
    f1.close()
   
    num_rows1 = sum(1 for line in open(spec_dat))
    x=np.zeros((num_rows1))
    y=np.zeros((num_rows1))
    for j in range(num_rows1):
        x[j]=(float(x2[j][0]))
        y[j]=(float(x2[j][2]))#using third column which sums over circle
    
    
    #cond_fe=(x>2550) & (x<2680)
    #ax3.plot (x[cond_fe],y[cond_fe])
    ax3.plot(x,y)
    ax3.set_xlabel("rest wavelength")
    ax3.set_ylabel("flux ")
    ax3.set_title("spectrum gal%s id %s"%(i+1,x1[i][1]))
    ax3.set_xlim(2550,2680)
    ax3.set_ylim(-500,500)
    ax3.axvline(Fe_lam1, color='r')
    ax3.axvline(Fe_lam2, color='r')
    ax3.axvline(Fe_lam1_st, color='g')
    ax3.axvline(Fe_lam2_st, color='g')
    #plt.clf()
    ax1.hold(False)
    ax2.hold(False)
    ax3.hold(False)
    return im1,

ani=anim.FuncAnimation(fig,animate,interval=1000, frames=num_rows)#init_func=init)

ani.save('animation.mp4', fps=1)

plt.show()
