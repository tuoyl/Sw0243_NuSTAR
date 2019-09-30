#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from tqdm import tqdm

def fsearch0243(data,fmin,fmax,f1,f2,fstep,bin_cs=20,bin_profile=20,t0=0):

    if t0 == 0:
        t0 =min(data)
    t_0 = t0
    #data = raw_data - min(raw_data);data.sort();
    raw_data = data
    #data = data - data[0]
#    t_0 = raw_data[0]
    N = len(data)
    # bin_cs=20 is DOF for chi_square test
    # bin_profile=20 is profile bin size
    b = N/bin_cs
    f = np.arange(fmin,fmax,fstep)
    #f1 = np.arange(f1min,f1max,f1step)
    #chi_square = [0] * len(f)
    chi_square = np.array([])

    for i in tqdm(range(0, len(f))):
        phi_tmp = np.mod((data-t_0)*f[i] + (1.0/2)*((data-t_0)**2)*f1 + (1.0/6.0)*((data-t_0)**3)*f2,1.0)
        p_num = np.histogram(phi_tmp,bin_cs)[0]
#        chi_square = np.append(chi_square,(np.std(p_num)**2/np.mean(p_num)))
        chi_square = np.append(chi_square, np.sum(((p_num-b)**2)/b))

#    for i in range(0,len(f)):
#        data = data - t_0
#        phi_tmp = np.mod(data*f[i] + (data**2)*f1*0.5 + (data**3)*f2/6,1.0)
#        p_num = np.histogram(phi_tmp,bin_cs)[0]
#        bb = b * np.ones(bin_cs)
#        chi_square[i] = np.sum((p_num-bb)**2)/b
#
#        percent = float(i)*100/len(f)
#        sys.stdout.write(" fsearch complete: %.2f"%percent);
#        sys.stdout.write("%\r");
#        sys.stdout.flush()
    fbest = f[np.where(chi_square==max(chi_square))][0]
    phi = np.mod((data-t_0)*fbest + ((data-t_0)**2)*f1*0.5 + ((data-t_0)**3)*f2/6,1.0)
    p_num = np.histogram(phi,bin_profile)[0]
#    p_num = p_num/np.mean(p_num)
    #p_num = p_num/(max(data)-min(data))
    p_num_x = np.arange(0.,bin_profile,1)/bin_profile

    p_num_x_2 = np.append(p_num_x, p_num_x+1)
    p_num_2 = np.append(p_num, p_num)
#    p_num_x_2_tmp = p_num_x + 1;p_num_x_2_tmp.tolist();
#    p_num_x_2 = p_num_x.tolist();p_num_2 = p_num.tolist();
#    p_num_x_2.extend(p_num_x_2_tmp);p_num_2.extend(p_num_2);

    print("fbest: ",fbest)
    print("T: ",1/fbest)
    print("t_0(the very first arrived photon): ",t_0)
    print("done")

    return p_num_x_2, p_num_2, f, chi_square, fbest

def create_gti(startlist, stoplist, outfile):
    START = startlist
    STOP  = stoplist
    
    # Table
    c1 = fits.Column(name='START', array=START,format='1D')
    c2 = fits.Column(name='STOP',  array=STOP, format='1D')
    tb = fits.BinTableHDU.from_columns([c1,c2])
    # Prmary Header
    header = fits.Header()
    primary_hdr = fits.Header()
    primary_hdr['comments'] = 'FITS (Flexible Image Transport System) format is defined in "Astronomyand Astrophysics", volume 376, page 359; bibcode: 2001A&A...376..359H'

    primary_hdu = fits.PrimaryHDU(header=primary_hdr)
    
    hdul = fits.HDUList([primary_hdu, tb])
    hdul.writeto(outfile,overwrite=True)
    
    # write Keywords
    TRUE = np.bool(True)
    FALSE = np.bool(False)
    hdulist = fits.open(outfile)

    hdulist[1].header['TUNIT1']= 's       '           
    hdulist[1].header['TUNIT2']= 's       '         
    hdulist[1].header['EXTNAME']= 'STDGTI  '           
    hdulist[1].header['MJDREFI']=   '55197 '
    hdulist[1].header['MJDREFF']=   '7.660185200000E-04 '
    hdulist[1].header['TIMEUNIT']= 's       '          
    hdulist[1].header['TIMEREF']= 'LOCAL   '           
    hdulist[1].header['CREATOR']= 'Xselect V2.4g'
    hdulist[1].header['HDUCLASS']= 'OGIP    '          
    hdulist[1].header['HDUCLAS1']= 'GTI     '          
    hdulist[1].header['HDUCLAS2']= 'STANDARD'          
    hdulist[1].header['TIMEZERO']=  '   0.0000000000E+00'
    hdulist[1].header['TSTART']=    str(min(START))
    hdulist[1].header['TSTOP']=     str(max(STOP))
    
    hdulist.writeto(outfile,overwrite=True)


##read events file
def read_file(filename):
    hdulist = fits.open(filename)
    time = hdulist[1].data.field("TIME")
    hdulist.close()
    return time



##calculate phases

##create time intervals for each phase


if __name__ == "__main__":
    filename = "/home/tuoyl/Work/Sw0243/nustar/data/nu90302319004_out/nu90302319004B01_cl.evt"
    time = read_file(filename)

    f = 0.10158457772688247
    f = 0.10158457772688252
    frange = 1.5e-5
    fstep = 1e-6
    pro_x, pro, fre, chi2, fbest = fsearch0243(time, f-frange, f+frange, 0,0, fstep, bin_cs=10, bin_profile=10, t0=0)

    t0 = min(time)
    phases =  np.mod((time-t0)*fbest, 1.0)

    period = 1/fbest

    time = np.sort(time)
    phasenum = 10

    for i in np.arange(0,phasenum):
        phase_i = i
            
        t_tmp = np.arange(min(time), max(time), period)
        tstart = t_tmp + (phase_i/phasenum) * period
        tstop  = tstart + ((phase_i+1)/phasenum) * period
        print(len(tstart), len(tstop))
        for i in range(len(tstart)):
            print(tstart[i], tstop[i])
    
        #create gti file 
        print(">>>>>>>creating GTI file<<<<<<<<")
        create_gti(tstart, tstop, "../data/B01gti_phase_%s-%s.fits"%(str(phase_i/phasenum),str((phase_i+1)/phasenum)))

    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(pro_x, pro)
    plt.subplot(2,1,2)
    plt.plot(fre, chi2)

    plt.show()
