#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys
from os import listdir
from os.path import isfile, join
from itertools import islice

#use this line for LS period search. # of cols to use depends on how many periods
#ls_period_array = np.genfromtxt('./LS_periods/LS_periods.data', usecols=(1,4,7,10,13))
#ls_period = [ f for f in listdir('/Volumes/ext_drive/work/KELT/jon_lcs/raw/ls_periods_month/') if isfile(join('/Volumes/ext_drive/work/KELT/jon_lcs/raw/ls_periods_month/',f))][1:]
#ls_periodogram = [ f for f in listdir('/Volumes/ext_drive/work/KELT/jon_lcs/raw/period_month/') if isfile(join('/Volumes/ext_drive/work/KELT/jon_lcs/raw/period_month/',f))]

trial = 'f'
subsample = '10.0'
sigma = '5'
minperiod = '0.01'
maxperiod = '10.0'

ls_periodogram = [ f for f in listdir('/media/sf_Astro/KELT/jon_lcs/tfa/LS/periodograms/{0}'.format(trial)) if isfile(join('/media/sf_Astro/KELT/jon_lcs/tfa/LS/periodograms/{0}'.format(trial),f))] #[1:]

ls_period = [ f for f in listdir('/media/sf_Astro/KELT/jon_lcs/tfa/LS/peaks/{0}'.format(trial)) if isfile(join('/media/sf_Astro/KELT/jon_lcs/tfa/LS/peaks/{0}'.format(trial),f))] #[1:]



#ls_period = ['/ABE-050/best_per_day_{0}.txt'.format(trial), '/ABE-075/best_per_day_{0}.txt'.format(trial)]
#ls_periodogram = 'w_N03_lc_41046.dat.ls' #ABE-075
#ls_periodogram = ['w_N11_lc_122551_{0}.dat.ls'.format(trial),'w_N03_lc_41046_{0}.dat.ls'.format(trial)] #ABE-050
#kelt_id = ['w_N11_lc_122551.dat','w_N03_lc_41046.dat']
#abenums = ['ABE-050','ABE-075']
#hd_num = ['HD 345439','HD 23478']

for (i,f) in enumerate(ls_periodogram):
#for (i,fname) in enumerate(ls_period[:4]):
    fname = f[:-3]
    data1 = np.genfromtxt('/media/sf_Astro/KELT/jon_lcs/tfa/LS/peaks/{0}/{1}.txt'.format(trial,fname),usecols=2) #data for phase folded plot with the best periods
    data1 = data1[1:]
    p_1 = data1[0]
    p_2 = data1[3]
    p_3 = data1[6]
    p_4 = data1[9]
    p_5 = data1[12]

    data2 = np.genfromtxt('/media/sf_Astro/KELT/jon_lcs/tfa/LS/periodograms/{0}/{1}.ls'.format(trial,fname)) #data for periodogram
    freq = data2[:,0]
    power = data2[:,1]
    period = 1/freq

    lc_data_comb = np.genfromtxt('/media/sf_Astro/KELT/jon_lcs/tfa/combined/{0}'.format(fname)) #data for lightcurve (time and magnitude)


    t0 = 1500.
    time_ = lc_data_comb[:,0]
    mag = lc_data_comb[:,1]
#   time_e=east_data[:,0]
#   mag_e = east_data[:,1]
#   time_w=west_data[:,0]
#   mag_w = west_data[:,1]
    time = time_ - t0
    ph1 = time / p_1 - np.floor(time/p_1)
    ph2 = time / p_2 - np.floor(time/p_2)
    ph3 = time / p_3 - np.floor(time/p_3)
    ph4 = time / p_4 - np.floor(time/p_4)
    ph5 = time / p_5 - np.floor(time/p_5)


#------------------begin averaging stuff--------------------
    #break up the ph1, ph2,... into 50 or so bins, and get the averages to plot on the same figure as the data points
    bin_size = 0.02
    Nbins = 50
    (bin_mag1,bin_mag2,bin_mag3,bin_mag4) = ([],[],[],[]) #initializes lists to hold mag info for each bin
    w_bin_mag = []
    w_med_list = []
    med_list = []
    middle_phase_list = []
    for m in range(1,5):
        exec( '(bin_mag{0},med_list{0},middle_phase_list{0}) = ([],[],[])'.format(m) )
        for j in range(Nbins):
            mid_phase =(j*bin_size + (j+1)*bin_size)/2
            for k in range(len(eval('ph{0}'.format(m)))):
                if (j*bin_size) < eval('ph{0}[k]'.format(m)) < ((j+1)*bin_size):
                    eval('bin_mag{0}'.format(m)).append(mag[k])
            if len(eval('bin_mag{0}'.format(m))) > 30: # and len(w_bin_mag) > 30:
                exec('bin_med{0} = np.median(bin_mag{0})'.format(m))
                exec('med_list{0}.append(bin_med{0})'.format(m))
                eval('middle_phase_list{0}'.format(m)).append(mid_phase)
            exec( 'bin_mag{0} = []'.format(m) )

        exec('middle_phase_array{0} = np.array(middle_phase_list{0})'.format(m))
    
#------------------end averaging stuff--------------------

    f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3,2,sharex=False, sharey=False)
    f.subplots_adjust(top = 0.9)
    f.suptitle('  LS Periodogram' + '  subsample= ' + subsample + ' non-whitened with ' + sigma + ' sigma clipping' \
               +  '   Trial ' + trial + '\n' + fname + '    Min Period = ' + minperiod + '   Max Period = ' + maxperiod + '\n')

    ax1.scatter(time/1000., mag, marker=',',s=0.05)
    ax1.set_title(' \n_')
    ax1.grid(True)
    ax1.set_xlabel('HJD / 1000',fontsize = 10)
    ax1.set_ylabel('Approx. V-mag')
    plt.xlim([time.min() - 20, time.max() + 20])
    plt.ylim([mag.min() - 0.01, mag.max() + 0.01])
    ax1.invert_yaxis()
#    plt.ylim([mag.max() + 0.01, mag.min() - 0.01])

    ax2.vlines(period,0,power,colors = 'k', linestyles = 'solid')
#    ax6.title('LS Periodogram for ' + str( ls_periodogram[i][:-3]) + ' ' + abenums[i] \
#              + '\n' + '  subsample = '+subsample + ' whitened with '+sigma+  ' sigma clipping' + '  Trial ' + trial + '\n')

    plt.ylim([0, power.max()])
    ax2.set_xlabel('Period (days)')
    ax2.set_ylabel('LS Power')
    ax2.set_title(' \n')
    #ax6.xscale('log')
    plt.xlim([0, period.max()])
    fig = plt.gcf()
#    ax = fig.gca()
    ax2.xaxis.grid(True)
    ax2.grid(True,which='both')

    ax3.plot(ph1, mag, 'k,',markersize=1.0)
    ax3.plot(ph1-1, mag, 'k,',markersize=2.0)
    ax3.plot(middle_phase_array1,med_list1,'r-')
    ax3.plot(middle_phase_array1-1,med_list1,'r-')
    ax3.set_title('Phase 1: ' + str(p_1) + ' days')
    ax3.ticklabel_format(axis='y', style='sci', scilimits=(-2,2), useOffset=False)
    ax3.grid(True)
    plt.xlim([-1.,1.])
#    plt.ylim([mag.max() + 0.01, mag.min() - 0.01])
    plt.ylim([mag.min() - 0.01, mag.max() + 0.01])
    ax3.invert_yaxis()
    

    ax4.plot(ph2, mag, 'k,',markersize=1.0)
    ax4.plot(ph2-1, mag, 'k,',markersize=2.0)
    ax4.plot(middle_phase_array2,med_list2,'r-')
    ax4.plot(middle_phase_array2-1,med_list2,'r-')
    ax4.set_title('Phase 2: ' + str(p_2) + ' days')
    ax4.ticklabel_format(axis='y', style='sci', scilimits=(-2,2), useOffset=False)
    ax4.grid(True)
    plt.xlim([-1.,1.])
    plt.ylim([mag.min() - 0.01, mag.max() + 0.01])
    ax4.invert_yaxis()
#    plt.ylim([mag.max() + 0.01, mag.min() - 0.01])

    ax5.plot(ph3, mag, 'k,',markersize=1.0)
    ax5.plot(ph3-1, mag, 'k,',markersize=2.0)
    ax5.set_title('Phase 3: ' + str(p_3) + ' days')
    ax5.plot(middle_phase_array3,med_list3,'r-')
    ax5.plot(middle_phase_array3-1,med_list3,'r-')
    ax5.ticklabel_format(axis='y', style='sci', scilimits=(-2,2), useOffset=False)
    ax5.grid(True)
    plt.xlim([-1.,1.])
    plt.ylim([mag.min() - 0.01, mag.max() + 0.01])
    ax5.invert_yaxis()
#    plt.ylim([mag.max() + 0.01, mag.min() - 0.01])

    ax6.plot(ph4, mag, 'k,',markersize=1.0)
    ax6.plot(ph4-1, mag, 'k,',markersize=2.0)
    ax6.set_title('Phase 4: ' + str(p_4) + ' days')
    ax6.plot(middle_phase_array4,med_list4,'r-')
    ax6.plot(middle_phase_array4-1,med_list4,'r-')
    ax6.ticklabel_format(axis='y', style='sci', scilimits=(-2,2), useOffset=False)
    ax6.grid(True)
    ax6.set_xlabel('Phase')
    plt.xlim([-1.,1.])
    plt.ylim([mag.min() - 0.01, mag.max() + 0.01])
    ax6.invert_yaxis()
#    plt.ylim([mag.max() + 0.01, mag.min() - 0.01])

    fig.set_size_inches(10.0,10.0)


    plt.tight_layout(pad=3.0, w_pad=1.5, h_pad = 1.5)

    plt.savefig('/media/sf_Astro/KELT/jon_lcs/tfa/LS/figs/{0}/{1}.png'.format(trial,fname))
    plt.close()


"""


    f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, sharex=False, sharey=True, figsize=(8,17))

    plt.grid(True)
    f.subplots_adjust(top = 0.82)
#   f.suptitle(str(filename), horizontalalignment='right')
    f.text(.1, .9, str( kelt_id[i]) + '   min per = 0.01d   ' + 'max per = 15.0d' + '\n'  \
           + ' subsample = ' + subsample + ' whitened with ' + sigma + ' sigma clipping' + '\n' \
           + 'HD ' + hd_num[i] + '    ' + abenums[i] + 'Trial ' + trial)
    #f.text(.1, .86, 'HD number, RA, Dec (SIMBAD): \n'+ simbad_hd_coords[i])
    plt.ylabel('Approx V-mag (KELT_mag - 3.9)')

#    ax1.plot(time_e, mag_e, 'r.', markersize=2.0)
#    ax1.plot(time_w, mag_w, 'b.',markersize=2.0)
    ax1.scatter(time, mag, s=2.0)
    ax1.set_title(ls_periodogram[i][:-3])
    ax1.grid(True)
    plt.xlim([time.min() - 20, time.max() + 20])
    plt.ylim([mag.max() + 0.1, mag.min() - 0.1])

    ax2.plot(ph1, mag, 'k.',markersize=2.0)
    ax2.plot(ph1-1, mag, 'k.',markersize=2.0)
    ax2.set_title('Phase: ' + str(p_1) + ' days')
    ax2.grid(True)
    plt.xlim([-1.,1.])
    plt.ylim([mag.max() + 0.1, mag.min() - 0.1])    

    ax3.plot(ph2, mag, 'k.',markersize=2.0)
    ax3.plot(ph2-1, mag, 'k.',markersize=2.0)
    ax3.set_title('Phase: ' + str(p_2) + ' days')
    ax3.grid(True)
    plt.xlim([-1.,1.])
    plt.ylim([mag.max() + 0.1, mag.min() - 0.1])

    ax4.plot(ph3, mag, 'k.',markersize=2.0)
    ax4.plot(ph3-1, mag, 'k.',markersize=2.0)
    ax4.set_title('Phase: ' + str(p_3) + ' days')
    ax4.grid(True)
    plt.xlim([-1.,1.])
    plt.ylim([mag.max() + 0.1, mag.min() - 0.1])

    ax5.plot(ph4, mag, 'k.',markersize=2.0)
    ax5.plot(ph4-1, mag, 'k.',markersize=2.0)
    ax5.set_title('Phase: ' + str(p_4) + ' days')
    ax5.grid(True)
    plt.xlim([-1.,1.])
    plt.ylim([mag.max() + 0.1, mag.min() - 0.1])

    plt.tight_layout()

    plt.savefig('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/Sigma_Ori_E/{0}/phased_lcs_ls/{0}_{1}.png'.format(abenums[i], trial))
    plt.cla()
    plt.clf()

    plt.vlines(period,0,power,colors = 'k', linestyles = 'solid')
    plt.title('LS Periodogram for ' + str( ls_periodogram[i][:-3]) + ' ' + abenums[i] \
              + '\n' + '  subsample = '+subsample + ' whitened with '+sigma+  ' sigma clipping' + '  Trial ' + trial + '\n')

    plt.ylim([0, power.max()])
    plt.xlabel('Period (days)')
    plt.ylabel('LS Power')
    plt.xscale('log')
    plt.xlim([0, period.max()])
    fig = plt.gcf()
    ax = fig.gca()
    ax.xaxis.grid(True)
    plt.grid(True,which='both')
    fig.set_size_inches(8.0,8.0)
    plt.tight_layout()
    plt.savefig('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/Sigma_Ori_E/{0}/periodograms_ls/{0}_{1}.png'.format(abenums[i], trial) )
    plt.clf()
    plt.cla()


    plt.vlines(period_clipped,0,power_clipped,colors = 'k', linestyles = 'solid')
    plt.title('Clipped LS Periodogram for ' + str( ls_periodogram[i][:-3]) + ' ' + abenums[i] \
              + '\n' + '  subsample = '+ subsample  + ' whitened with ' + sigma + ' sigma clipping' + '  Trial ' + trial + '\n')

    plt.ylim([0, power.max()])
    plt.xlabel('Period (days)')
    plt.ylabel('LS Power')
    plt.xscale('log')
    plt.xlim([0, period.max()])
    fig = plt.gcf()
    ax = fig.gca()
    ax.xaxis.grid(True)
    plt.grid(True,which='both')
    plt.savefig('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/Sigma_Ori_E/{0}/periodograms_ls_clipped/{0}_{1}.png'.format(abenums[i], trial) )
    plt.clf()
    plt.cla()
    plt.close()
"""
