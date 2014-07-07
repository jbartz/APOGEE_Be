import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join

#lcs_e = [ f for f in listdir('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/east') if isfile(join('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/east',f))]
#lcs_w = [ f for f in listdir('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/west') if isfile(join('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/west',f))]
lcs_combined = [ f for f in listdir('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/combined/') if isfile(join('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/combined/',f))]
vartools = '/Volumes/Payne/Users/jon/VARTOOLS1.202/vartools'
#used the combined east+west lightcurves and run a LS period search on each lc

trial = 'f'
min_p = '0.01'
max_p = '10.0'
subsample = '0.1'
sigma = '5'

input_dir = '/Volumes/ext_drive/work/KELT/jon_lcs/tfa/combined/'
output_dir = '/Volumes/ext_drive/work/KELT/jon_lcs/tfa/LS/'
output_dir_peaks = '/Volumes/ext_drive/work/KELT/jon_lcs/tfa/LS/peaks/'

ls_outfile = '/Volumes/ext_drive/work/Scripts/ls_{0}.run'.format(trial)

with open(ls_outfile, 'a') as outputfile:
    for obj in lcs_combined:
        line = vartools + ' -i ' + input_dir + obj + ' -oneline -ascii -LS ' + min_p + ' ' + max_p + ' ' + subsample + ' 15 1 ' + output_dir + 'periodograms/' \
               + trial +'/ whiten clip ' + sigma + ' 1 > ' + output_dir_peaks + trial + '/' + str(obj) + '.txt' + '\n'
        outputfile.write(line)
        
outputfile.close()
