#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import glob
from matplotlib.ticker import MultipleLocator,AutoMinorLocator,FormatStrFormatter
from matplotlib import gridspec

def ratio_err(u,x,y,err_x, err_y):
    tmp = (err_x**2/x**2 + err_y**2/y**2)*(u**2)
    err_u = np.sqrt(tmp)
    return err_u

#phalist = glob.glob("../data/phase_resolved/B01_products_0*/nu90302319004B01_sr.pha")
phalist = sorted(glob.glob("../data/phase_resolved/A01_products_0*/nu90302319004A01_sr.pha"))

counts_sum  = [] 
error_sum   = []
exposure_sum = np.array([])
onsoure_sum  = np.array([])

fig = plt.figure(figsize=(11, 8))
gs = gridspec.GridSpec(2,1,height_ratios=[2,1])
fig.subplots_adjust(hspace=0)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1], sharex=ax1)
ax1.xaxis.set_ticks_position('none') 

for pha in phalist:
    print(">> pha : %s"%pha)

    hdulist = fits.open(pha)
    exp     = np.float(hdulist[1].header["EXPOSURE"])
    ont     = np.float(hdulist[1].header["ONTIME"])
    print(exp,ont)
    channel = hdulist[1].data.field("channel")
    counts  = hdulist[1].data.field("counts")
    error   = np.sqrt(counts)

    counts_sum.append(counts)
    error_sum.append(error)
    exposure_sum = np.append(exposure_sum, exp)
    onsoure_sum  = np.append(onsoure_sum, ont)
    #counts_sum = counts_sum + counts
    #error_sum  = np.sqrt(error_sum**2 + (error)**2)
    #exposure_sum = exposure_sum + exp
    #onsoure_sum  = onsoure_sum + ont

#    ax1.errorbar(channel, counts, yerr=error,color='black',capsize=0)

#normalization on exposure
counts_sum_norm = np.zeros(len(channel))
error_sum_norm  = np.zeros(len(channel))
for i in range(len(counts_sum)):
    counts_sum_norm += counts_sum[i]/exposure_sum[i]
    error_sum_norm  += (error_sum[i]/exposure_sum[i])**2
counts_sum_norm= counts_sum_norm/10
error_sum_norm = np.sqrt(error_sum_norm)/10
ax1.errorbar(channel, counts_sum_norm, xerr=0.5, yerr=error_sum_norm, fmt='o', color='blue',mfc="none",mec='blue', label="sum",capsize=0)
#non-norm
#ax1.errorbar(channel_sum, counts_sum,xerr=0.5, yerr=error_sum, fmt='o', color='blue',mfc="none",mec='blue', label="sum (ontime=%s s, livetime=%s s)"%(str(onsoure_sum),str(exposure_sum)),capsize=0)

##--- sum spec by addspec
pha = glob.glob("/Users/tuoyouli/Work/Sw0243/nustar/data/phase_resolved/merge_spec/A01total_sr.pha")[0]
hdulist = fits.open(pha)
counts_sum_addspec = hdulist[1].data.field("COUNTS")
exposure_sum_addspec = np.float(hdulist[1].header["exposure"])
ax1.errorbar(channel, counts_sum_addspec/exposure_sum_addspec,fmt='-.',capsize=0)

#pha all
#pha = glob.glob("../data/phase_resolved/B01_products_all/nu90302319004B01_sr.pha")[0]
pha = glob.glob("../data/phase_resolved/A01_products_all/nu90302319004A01_sr.pha")[0]
print(">> pha all: %s"%pha)
hdulist = fits.open(pha)
channel_all = hdulist[1].data.field("CHANNEL")
counts_all  = hdulist[1].data.field("COUNTS")
error_all   = np.sqrt(counts_all)
exposure_all= np.float(hdulist[1].header["EXPOSURE"])
onsource_all = np.float(hdulist[1].header["ONTIME"])
counts_all_norm = counts_all/exposure_all
error_all_norm  = error_all/exposure_all
#ax1.errorbar(channel_all, counts_all, xerr=0.5, yerr=error_all, color='red', mfc="none", mec='red',fmt='o', capsize=0,
#        label="total (ontime=%s s, livetime=%s s)"%(str(onsource_all),str(exposure_all)))
ax1.errorbar(channel_all, counts_all_norm, xerr=0.5, yerr=error_all_norm, color='red', mfc="none", mec='red',fmt='o', capsize=0,
        label="total")


#ax2.axhline(y=0,color='black',lw=0.5)
#ax2.errorbar(channel, counts_sum_norm-counts_all_norm,xerr=0.5,
#        yerr=np.sqrt(error_sum_norm**2 + error_all_norm**2) ,color='black',capsize=0,fmt='.',ms=0)
ax2.axhline(y=1, color='black', lw=0.5)
error_ratio_norm =  ratio_err(counts_all_norm/counts_sum_norm, counts_all_norm, counts_sum_norm, error_all_norm, error_sum_norm)
ax2.errorbar(channel_all, counts_all_norm/counts_sum_norm, xerr=0.5, yerr=error_ratio_norm, capsize=0,color='black',ms=0,fmt='.')


print(">>> cmp expo (sum, all): ", exposure_sum, exposure_all)
print(">>> cmp onsource (sum, all): ", onsoure_sum, onsource_all)

ax1.tick_params(axis='both',direction='in',which='both',top=True, bottom=True, left=True, right=True, labelbottom=False)
plt.setp(ax1.get_yticklabels()[0],visible=False) #hide the last tick marker of y-axis to avoid overlaps
plt.setp(ax1.get_yticklabels()[-1],visible=False) #hide the last tick marker of y-axis to avoid overlaps
ax1.set_ylabel("Counts", fontsize=15)
ax2.set_ylabel("$\delta$", fontsize=15)
ax2.set_xlabel("Channel",fontsize=15)
#ax1.set_yscale('log',nonposy='clip')
ax1.set_xscale('log',nonposy='clip')
#ax2.set_yscale('log',nonposy='clip')
ax2.set_xscale('log',nonposy='clip')

ax1.set_xlim([.9,4100])
ax2.set_xlim([.9,4100])
#ax1.set_ylim([9,3000])
#ax1.set_ylim([-20,500])
#ax1.set_ylim([-20,500])
#ax1.set_ylim([-0.5,3.5])
#ax1.set_ylim([-0.5,3.5])


ax1.legend(numpoints=1)
#plt.savefig("../results/B01phase_10_spec_cmp.png")
plt.savefig("../results/A01phase_10_spec_cmp.png")
plt.show()
