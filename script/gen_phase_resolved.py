import glob
import numpy as np











infile = "./nu90302319004_out/nu90302319004B01_cl.evt"
srcreg = "nu90302319004B01_src.reg"
bkreg  = "nu90302319004B01_bkg.reg"
indir  = "nu90302319004_out"
stem   = "nu90302319004"
instru = "FPMB"
phase_num = 10
for i in np.arange(0,10):
    phase_i = i/phase_num
    phase_i1= (i+1)/phase_num
    outdir = "./phase_resolved/%s01_products_%s-%s"%(instru[-1],str(phase_i),str(phase_i1))
    gti    = "%s01gti_phase_%s-%s.fits"%(instru[-1], str(phase_i), str(phase_i1))

    test_cmd = "nuproducts infile=%s "\
            "srcregionfile=%s bkgregionfile=%s "\
            "indir=%s outdir=%s "\
            "instrument=%s steminputs=%s bkgextract=yes clobber=yes "\
            "usrgtifile=%s "%(infile, srcreg, bkreg, indir, outdir, instru, stem, gti)

    print(test_cmd)
