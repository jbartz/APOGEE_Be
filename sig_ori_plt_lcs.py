#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys
from os import listdir
from os.path import isfile, join
from itertools import islice

kelt_id = ['w_N11_lc_122551.dat','w_N03_lc_41046.dat']
abenums = ['ABE-050','ABE-075']
hd_num = ['HD 345439','HD 23478']
ra_ = ['299.70087', '56.66858']
dec_ = ['23.09003','32.28999']
vmag_simbad = ['11.11', '6.688']
x_coord=['3202.197','1668.470']
y_coord=['744.317','2144.353']


for (i,fname) in enumerate(abenums):
    #lc_data = np.genfromtxt('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/west/{0}'.format(kelt_id[i])) #data for lightcurve (time and magnitude)
    lc_data = np.genfromtxt('/media/sf_Astro/KELT/jon_lcs/tfa/west/{0}'.format(kelt_id[i])) #data for lightcurve (time and magnitude)
    time_ = lc_data[:,0]
    mag = lc_data[:,1]
    
    fig = plt.figure()
    fig.set_size_inches(8.0,5.0)
    plt.scatter(time_, mag, marker=',', s=0.5)
    plt.gca().invert_yaxis()
    
    plt.xlabel('HJD')
    plt.ylabel('Approx. V magnitude')
    plt.grid(True)
    plt.xlim((time_.min()-50., time_.max() + 50))
    plt.ylim((mag.min(), mag.max()))
    fig.patch.set_facecolor('white')
    #fig.suptitle(fname + '     ' + kelt_id[i] + '     ' + hd_num[i] + '     x: ' + x_coord[i] + ' \nRA: ' + ra_simbad[i] + '    Dec: ' + \
    #                   dec_simbad[i] + '   V_mag: ' + vmag_simbad[i] + '     y: ' + y_coord[i])
    
    fig.suptitle(fname + '     ' + kelt_id[i][:-4] + '     RA: '+ ra_[i] + '     x: ' + x_coord[i] + ' \n' + hd_num[i] + '    V_mag: ' + \
                       vmag_simbad[i] + '     Dec: ' + dec_[i] +  '     y: ' + y_coord[i])
    
    
    #plt.tight_layout(pad=2.0, w_pad=0.5, h_pad=2.0)
    #plt.show()
    plt.savefig('/media/sf_Astro/KELT/jon_lcs/tfa/Sigma_Ori_E/{0}/{0}_lc.png'.format(abenums[i]))
    plt.cla()
    plt.clf()
    