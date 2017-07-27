import numpy as np
from nexusformat.nexus import *
from matplotlib import pyplot as plt

def static_mask(nx=516,ny=1556,npadsx=2,npadsy=6,bordersize=8,padspacing = 8):
    '''
    Masking the border and the large pixels between pads.

    Input
    -----
    nx,ny: number of pixels along x and y (int)
    npadsx,npadsy: number of pads along x and y (int)
    bordersize: how many pixesl to mask from the edges
    padspacing: how many pixles to mask between pads

    Output
    ------
    mask (2x2 numpy array, boolean)
    '''
    mask = np.ones((nx,ny),dtype=bool)

    #Mask edges
    mask[:bordersize,:], mask[:,:bordersize] = 0,0
    mask[-bordersize:,:], mask[:,-bordersize:] = 0,0

    #mask pixels between pads
    dx = int(nx/npadsx)
    for i in range(0,npadsx):
        mask[dx*i-padspacing/2:dx*i+padspacing/2,:] = 0
    dy = int(ny/npadsy)+1#+padspacing/2
    for j in range(0,npadsy):
        mask[:,dy*j-padspacing/2-2:dy*j+padspacing/2-2] = 0

    print 'static: masked %.2f percent'%(len(mask[mask==0])/float(nx*ny)*100.)
    return mask

def beam_blocker_mask(nx=516,ny=1556,center_x=277,center_y=725,radius=30,slope=0.007,offset=-8,thickness=15):
    '''
    Masking the beam blocker
    '''
    mask = np.ones((nx,ny),dtype=bool)

    #disk in the center
    x,y= np.mgrid[0:nx,0:ny]
    x-=center_x
    y-=center_y
    rho = np.sqrt(x**2+y**2)
    mask[rho<radius]=0

    # line
    line1 = (1./slope)*(x+offset)-y
    line2 = line1+(1./slope)*thickness
    line = np.ones((nx,ny),dtype=int)
    line[:,center_y:] = line1[:,center_y:]*line2[:,center_y:]
    mask[line<0]=0 
    print 'beam blocker : masked %.2f percent'%(len(mask[mask==0])/float(nx*ny)*100.)

    return mask

def mean_std_mask(data,nx=516,ny=1556, min_thr_mean=1e-3, max_thr_mean = 1e5, plot_flag = True):
#npadsx=2,npadsy=6,bordersize=8,padspacing = 8):
    '''
    Masking using thresholds on the mean intensity and standard
    deviation per pixel 

    Input
    -----
    data: data array (numpy nd array: [t,x,y])
    nx,ny: number of pixels along x and y (int)
    min_thr_mean: minimum threshold for the mean cut off
    max_thr_mean: maximum threshold for the mean cut off
    min_thr_sd: minimum threshold for the standard deviation cut off
    max_thr_sd: maximum threshold for the standard deviation cut off

    Output
    ------
    mask (2x2 numpy array, boolean)
    '''

    # -- Parameters
    #min_thr_mean, max_thr_mean = 0.4,1000
    #min_thr_mean, max_thr_mean = 1e-3,1e5
    min_thr_sd,max_thr_sd = np.sqrt(min_thr_mean),np.sqrt(max_thr_mean)

    # -- Definitions
    mean_data = np.average(data,axis=0)
    sd_data = np.std(data,axis=0)
    mask_mean =np.ones((nx,ny),dtype=bool) 
    mask_mean[mean_data>max_thr_mean]=0
    mask_mean[mean_data<min_thr_mean] = 0

    mask_sd =np.ones((nx,ny),dtype=bool) 
    mask_sd[sd_data>max_thr_sd]=0 
    mask_sd[sd_data<min_thr_sd] = 0

    mask = mask_mean*mask_sd
    if plot_flag is True:
        plot_mean_std_mask(mean_data,sd_data,max_thr_mean,min_thr_mean,max_thr_sd,min_thr_sd,mask)

    #print 'mean & sd: masked %.2f percent'%(len(mask[mask==0])/float(nx*ny)*100.)
    print 'mean & sd: masked %d pixels'%(len(mask[mask==0]))
    return mask


def plot_mean_std_mask(mean_data,sd_data,max_thr_mean,min_thr_mean,max_thr_sd,min_thr_sd,mask):
    '''
    Assistive function for choosing thresholds via plotting the mean
    and standard deviation masks
    '''
    plt.figure()
    plt.subplot(3,1,1)
    #plt.imshow(mean_data,vmax= 5,vmin=0.4,interpolation='none')
    plt.imshow(np.log10(mean_data),vmax= 1,vmin=-3,interpolation='none')
    plt.title('mean')
    plt.colorbar()

    plt.subplot(3,1,2)
    plt.imshow(sd_data,vmax=2.5,vmin=0.6,interpolation='none')
    plt.imshow(np.log10(sd_data),vmax= 1,vmin=-1,interpolation='none')
    plt.title('standard deviation')
    plt.colorbar()

    plt.subplot(3,1,3)
    plt.imshow(mask,vmax=1,vmin=0,interpolation='none')
    plt.title('mask')
    plt.colorbar()
    plt.tight_layout()


    plt.figure()
    plt.subplot(2,1,1)
    bi,bf,db = 0.01,1000,0.01
    hy,hx = np.histogram(mean_data,bins = np.arange(bi,bf,db)) 
    plt.plot(hx[:-1],hy,label='raw',color='blue')
    plt.yscale('log',nonposy='clip')
    plt.xscale('log',nonposy='clip')
    plt.title('mean')
    plt.xlabel('Mean intensity (photons)')
    plt.ylabel('Number of pixels')
    plt.ylim([1e-1,1e5])
    #plt.axvline(x=min_thr_mean,color='red',ls='--')
    plt.axvline(x=max_thr_mean,color='red',ls='--')

    plt.subplot(2,1,2)
    bi,bf,db = 0.01,1000,0.01
    hy,hx = np.histogram(sd_data,bins = np.arange(bi,bf,db)) 
    plt.plot(hx[:-1],hy,label='raw',color='green')
    plt.yscale('log',nonposy='clip')
    plt.xscale('log',nonposy='clip')
    plt.title('standard deviation')
    plt.xlabel('Standard deviation (dphotons)')
    plt.ylabel('Number of pixels')
    plt.tight_layout()
    plt.ylim([1e-1,1e5])
    #plt.axvline(x=min_thr_sd,color='red',ls='--')
    plt.axvline(x=max_thr_sd,color='red',ls='--')

    #plt.show()
    return




