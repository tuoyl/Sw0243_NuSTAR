from __future__ import division
from astropy.io import fits
import numpy as np
import glob
import os

def get_time_intervals(pipelineout):

    MJDREFI = 55197
    MJDREFF = 7.6601852000000E-04

    FPMA_evt = glob.glob(pipelineout+"/nu*A01_cl.evt")[0]
    FPMB_evt = glob.glob(pipelineout+"/nu*B01_cl.evt")[0]


    hdulist = fits.open(FPMA_evt)
    time_A = hdulist[1].data.field("TIME")

    hdulist = fits.open(FPMB_evt)
    time_B = hdulist[1].data.field("TIME")

    split_sec = 10

    tmp = np.arange(min(time_A), max(time_A), split_sec)
    tmp = np.append(tmp, max(time_A))
    tstart_A = tmp[:-1]
    tstop_A   = tmp[1:]

    tmp = np.arange(min(time_B), max(time_B), split_sec)
    tmp = np.append(tmp, max(time_B))
    tstart_B = tmp[:-1]
    tstop_B   = tmp[1:]

    gti_list = [tstart_A, tstop_A, tstart_B, tstop_B]
    return gti_list

def get_gti_file(gti_list,gti_template,outdir):
    #create gti file
    hdulist = fits.open(gti_template)

    startA = gti_list[0]
    stopA  = gti_list[1]
    startB = gti_list[2]
    stopB  = gti_list[3]

    #A
    for i in range(len(startA)):
        print("%s interval"%startA[i])
        create_gti(startA[i], stopA[i], outfile=os.path.join(outdir,"A01_gti_%s.fits"%str(i)))
    #B
    for i in range(len(startB)):
        print("%s interval"%startB[i])
        create_gti(startB[i], stopB[i], outfile=os.path.join(outdir,"B01_gti_%s.fits"%str(i)))


def create_gti(start, stop, outfile):
    START = np.array([start])
    STOP  = np.array([stop])


    
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
#    hdulist[1].header['EXTNAME'] = 'SPECTRUM'
#    hdulist[1].header['DETCHANS'] = 256
#    hdulist[1].header['BACKFILE'] = 'NONE'
#    hdulist[1].header['BACKSCAL'] = 1
#    hdulist[1].header['CORRFILE'] = 'NONE'
#    hdulist[1].header['CORRSCAL'] = 0
#    hdulist[1].header['RESPFILE'] = 'NONE'
#    hdulist[1].header['ANCRFILE'] = 'NONE'
#    hdulist[1].header['FILETER']  = 'NONE'
#    hdulist[1].header['PHAVERSN'] = '1992a'
#    hdulist[1].header['STATERR']  = TRUE
#    hdulist[1].header['SYSERR']   = FALSE
#    hdulist[1].header['POISSERR'] = FALSE
#    hdulist[1].header['GROUPING'] = 1
#    hdulist[1].header['QUALITY']  = 1
#    hdulist[1].header['AREASCAL'] = 1
#    hdulist[1].header['EXPOSURE'] = exposure
#    hdulist[1].header['LIVETIME'] = 1
#    hdulist[1].header['DEADC']    = 0
#    hdulist[1].header['DETID']    = 0
#    hdulist[1].header['CHANTYPE'] = 'PI'
#    hdulist[1].header['TLMIN2']   = 0
#    hdulist[1].header['TLMAX2']   = 255
#    hdulist[1].header['TELESCOP'] = 'HXMT'
#    hdulist[1].header['INSTRUME'] = 'HE'


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
    hdulist[1].header['TSTART']=    str(start)         
    hdulist[1].header['TSTOP']=     str(stop)          
    
    hdulist.writeto(outfile,overwrite=True)

def get_nuproducts_cmd(gtidir, outdir, obsid):
    cmd = []
    gtifiles = glob.glob(gtidir+"/A01_gti_*.fits")
    for i in range(len(gtifiles)):
        gtiname   = os.path.join(gtidir , "A01_gti_%s.fits"%str(i))
        print("LLLLLLLLL",outdir)
        outfolder = os.path.join(outdir , "A01_products_%s"%str(i))
        nuproducts_cmd = "nuproducts infile=./nu%s_out/nu%sA01_cl.evt "\
                "srcregionfile=nu%sA01_src.reg bkgregionfile=nu%sA01_bkg.reg "\
                "indir=nu%s_out outdir=%s "\
                "instrument=FPMA steminputs=nu%s bkgextract=yes clobber=yes  "\
                "usrgtifile=%s binsize=0.5"%(obsid, obsid, obsid, obsid, obsid, outfolder, obsid, gtiname)
        cmd.append(nuproducts_cmd)
    gtifiles = glob.glob(gtidir+"/B01_gti_*.fits")
    for i in range(len(gtifiles)):
        gtiname   = os.path.join(gtidir , "B01_gti_%s.fits"%str(i))
        outfolder = os.path.join(outdir , "B01_products_%s"%str(i))
        nuproducts_cmd = "nuproducts infile=./nu%s_out/nu%sB01_cl.evt "\
                "srcregionfile=nu%sB01_src.reg bkgregionfile=nu%sB01_bkg.reg "\
                "indir=nu%s_out outdir=%s "\
                "instrument=FPMB steminputs=nu%s bkgextract=yes clobber=yes  "\
                "usrgtifile=%s binsize=0.5"%(obsid, obsid, obsid, obsid, obsid, outfolder, obsid, gtiname)
        cmd.append(nuproducts_cmd)
    with open("tmp.sh",'w')as fout:
        for i in range(len(cmd)):
            fout.write(cmd[i]+'\n')
    

if __name__ == "__main__":
    obsid = "90302319004"
    pipelineout = "nu90302319004_out"
    gti_template = "gti_template.fits"
    gti_outdir  = "gtifiles"
    products_dir = "./split_products"

    gti_list = get_time_intervals(pipelineout)
    get_gti_file(gti_list, gti_template, gti_outdir)
    get_nuproducts_cmd(gti_outdir, products_dir, obsid)


