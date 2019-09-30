#!/usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import glob

def create_nu_spec_file(counts, exposure, outfile, instrument="FPMA"):
    Channel = np.arange(0,4096)
    COUNTS = counts
    
    # Table
    c1 = fits.Column(name='Channel', array=Channel,format='1J')
    c2 = fits.Column(name='COUNTS' , array=COUNTS ,format='1J')
    tb = fits.BinTableHDU.from_columns([c1,c2])
    # Prmary Header
    header = fits.Header()
    primary_hdr = fits.Header()
    primary_hdr['comments'] = 'FITS(Flexible Image Transport System)'
    primary_hdu = fits.PrimaryHDU(header=primary_hdr)
    
    hdul = fits.HDUList([primary_hdu, tb])
    hdul.writeto(outfile,overwrite=True)
    
    # write Keywords
    T= np.bool(True)
    F= np.bool(False)
    hdulist = fits.open(outfile)

    hdulist[1].header['EXTEND  '] = T
    hdulist[1].header['CONTENT '] = 'PHA SPECTRUM'       
    hdulist[1].header['XTENSION'] = 'BINTABLE'           
    hdulist[1].header['PCOUNT  '] = 0 
    hdulist[1].header['GCOUNT  '] = 1 
    hdulist[1].header['TFIELDS '] = 3 
    hdulist[1].header['TTYPE1  '] = 'CHANNEL '           
    hdulist[1].header['TFORM1  '] = 'J       '           
    hdulist[1].header['TTYPE2  '] = 'COUNTS  '           
    hdulist[1].header['TFORM2  '] = 'J       '           
    hdulist[1].header['TUNIT2  '] = 'count   '           
    hdulist[1].header['EXTNAME '] = 'SPECTRUM'           
    hdulist[1].header['HDUCLASS'] = 'OGIP    '           
    hdulist[1].header['HDUCLAS1'] = 'SPECTRUM'           
    hdulist[1].header['HDUVERS1'] = '1.2.0   '           
    hdulist[1].header['HDUVERS '] = '1.2.0   '           
    hdulist[1].header['HDUCLAS2'] = 'DERIVED '           
    hdulist[1].header['HDUCLAS3'] = 'COUNT   '           
    hdulist[1].header['TLMIN1  '] = 0 
    hdulist[1].header['TLMAX1  '] = 4095 
    hdulist[1].header['TELESCOP'] = 'NuSTAR  '           
    hdulist[1].header['INSTRUME'] = instrument           
    hdulist[1].header['DETNAM  '] = 'NONE    '           
    hdulist[1].header['FILTER  '] = 'NONE    '           
    hdulist[1].header['EXPOSURE'] = exposure
    hdulist[1].header['AREASCAL'] = 1.000000E+00 
    hdulist[1].header['BACKFILE'] = 'NONE    '
    hdulist[1].header['BACKSCAL'] =  1.0 
    hdulist[1].header['CORRFILE'] = 'NONE    '  
    hdulist[1].header['CORRSCAL'] =  1.000000E+00
    hdulist[1].header['RESPFILE'] = 'NONE    '
    hdulist[1].header['ANCRFILE'] = 'NONE    '           
    hdulist[1].header['PHAVERSN'] = '1992a   '           
    hdulist[1].header['DETCHANS'] = 4096 
    hdulist[1].header['CHANTYPE'] = 'PI      '           
    hdulist[1].header['POISSERR'] = T 
    hdulist[1].header['SYS_ERR '] = 0 
    hdulist[1].header['GROUPING'] = 0 

    hdulist.writeto(outfile,overwrite=True)

if __name__ == "__main__" :
    #get files
#    filelist = glob.glob("/Users/tuoyouli/Work/Sw0243/nustar/data/phase_resolved/A01_products_0*/nu90302319004A01_sr.pha")
    filelist = glob.glob("/Users/tuoyouli/Work/Sw0243/nustar/data/phase_resolved/A01_products_0*/nu90302319004A01_bk.pha")
    counts_all = []
    exposure_all = np.array([])
    for filename in filelist:
        hdulist = fits.open(filename)
        channel = hdulist[1].data.field("CHANNEL")
        counts  = hdulist[1].data.field("COUNTS")
        exposure= np.float(hdulist[1].header["EXPOSURE"])

        counts_all.append(counts)
        exposure_all = np.append(exposure_all, exposure)

    #sum up
    counts_sum = np.zeros(len(channel))
    exposure_sum = 0
    for i in range(len(counts_all)):
        counts_sum += counts_all[i]/exposure_all[i] #NOTE:counts average by expo
        exposure_sum += exposure_all[i]
    counts_sum = counts_sum/len(counts_all) #NOTE:rate for 10 bins
    counts_sum = counts_sum*exposure_sum    #NOTE:calculate counts on sum exposure

    print(counts_sum, exposure_sum)

    #create spec file
#    create_nu_spec_file(counts_sum, exposure_sum, "../data/phase_resolved/merge_spec/my_A01_merge_sr.pha",instrument="FPMA")
    create_nu_spec_file(counts_sum, exposure_sum, "../data/phase_resolved/merge_spec/my_A01_merge_bk.pha",instrument="FPMA")



