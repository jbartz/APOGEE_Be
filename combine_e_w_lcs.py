import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
#if an object has both east and west lightcurves, this script will combine the two, hopefully aligning the medians or something correctly

lc_files = [ f for f in listdir('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/all/') if isfile(join('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/all/',f))]
match_key = np.genfromtxt('/Volumes/ext_drive/work/a_k_match_trim.txt',dtype=str)
id_e_list = match_key[:,0]
id_w_list = match_key[:,1]
id_a_list = match_key[:,2]
#west_new = np.genfromtxt('./be_matches_kelt_id.list',dtype=str)

west_files=[]
east_files=[]
fields=[]
for i in range(len(match_key)):
#for i in range(2):
    e_id = match_key[i][0]
    w_id = match_key[i][1]
    a_id = match_key[i][2]
    if (e_id != 'NULL'): #if this object has an east lightcurve
        field = e_id[2:5]
        lc_num_e = e_id[6:-4] + '_mag.data'
        #e_data = np.genfromtxt('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/all/{0}'.format(lc_num_e))
        try:
            e_data = np.genfromtxt('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/all/{0}'.format(lc_num_e))
        except IOError:
            e_data = 'NULL'
            lc_num_e = 'NULL'
            e_time = 'NULL'
            e_mag = 'NULL'
        else:
            e_time = e_data[:,0]
            e_mag = e_data[:,1]
        fields.append(field)
    else: #if this object doesn't have an east lightcurve
        e_data = 'NULL'
        lc_num_e = 'NULL'
        e_time = 'NULL'
        e_mag = 'NULL'
        
    if (w_id != 'NULL'): #if this object has a west lightcurve
        field = w_id[2:5]
        lc_num_w = w_id[6:-4] + '_mag.data'
        #w_data = np.genfromtxt('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/all/{0}'.format(lc_num_w))  
        #new code here!!!-----------
        try:
            w_data = np.genfromtxt('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/all/{0}'.format(lc_num_w))  
            #break
        except IOError:
            w_data = 'NULL'
            lc_num_w = 'NULL'
            w_time = 'NULL'
            w_mag = 'NULL'
        else:
        #end new code--------

            w_time = w_data[:,0]
            w_mag = w_data[:,1]
    else: #if this object doesn't have a west lightcurve
        w_data = 'NULL'
        lc_num_w = 'NULL'
        w_time = 'NULL'
        w_mag = 'NULL'
    ##so far, we've loaded the data for the east and/or west lightcurve for this object. Need to plot east and west data points on one plot
    ##if both != 'NULL': proceed as below
    if (e_data != 'NULL') & (w_data != 'NULL'):
        if e_time.min() < w_time.min():
            t0 = e_time.min()
        else:
            t0 = w_time.min()
        if e_time.max() < w_time.max():
            tmax = w_time.max()
        else:
            tmax = e_time.max()
        bin_size = 10
        Nbins = (int(tmax + 1) - int(t0))/bin_size
        e_bin_mag = [] #initializes lists to hold mag info for each bin
        w_bin_mag = []
        w_med_list = []
        e_med_list = []
        middle_time_list = []
        for j in range(Nbins + 1):
            for k in range(e_data.shape[0]):
                if (t0 + j*bin_size) < e_data[k,0] < (t0 + (j+1)*bin_size):
                    e_bin_mag.append(e_data[k,1])
            for k in range(w_data.shape[0]):
                if (t0 + j*bin_size) < w_data[k,0] < (t0 + (j+1)*bin_size):
                    w_bin_mag.append(w_data[k,1])
            if len(e_bin_mag) > 30 and len(w_bin_mag) > 30:
                w_bin_med = np.median(w_bin_mag)
                e_bin_med = np.median(e_bin_mag)
                w_med_list.append(w_bin_med)
                e_med_list.append(e_bin_med)
                middle_time = t0 + (j*bin_size + (j+1)*bin_size)/2
                middle_time_list.append(middle_time)
            e_bin_mag = []
            w_bin_mag = []
        w_allbins_med = np.median(w_med_list) - 3.9
        e_allbins_med = np.median(e_med_list) - 3.9
        offset = w_allbins_med - e_allbins_med
        w_data[:,1] -= offset
    
        e_data[:,1] -= 3.9 #puts mag on roughly v-mag system
        w_data[:,1] -= 3.9
        #outfile = east_files[i][:2] + east_files[i][3:]
    
        #outputs a lc with both east and west data in one file
        lc_data = np.vstack([e_data,w_data])
        np.savetxt('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/combined/{0}'.format(e_id),lc_data)
    #    np.savetxt('./east_lcs/'.format(i),east_data)
    #    np.savetxt('./west_lcs/04w_lc_023161.raw.mag',west_data)
    
        e_med_array=np.array(e_med_list)- 3.9
        w_med_array = np.array(w_med_list) - 3.9 - offset
        #plots mag vs. time and saves each lc
        plt.scatter(e_data[:,0],e_data[:,1],s=0.5,marker='.',color='red',label='east')
        plt.scatter(w_data[:,0],w_data[:,1],s=0.5,marker='.',color='blue',label='west')
    #    plt.scatter(middle_time_list,e_med_array,s=10.0, marker='o',color='green')
    #    plt.scatter(middle_time_list,w_med_array,s=10.0, marker='o',color='black')
    #    plt.axhline(y=(w_allbins_med - offset), color='blue')
    #    plt.axhline(y=e_allbins_med, color='red', ls='dashed')
        plt.gca().invert_yaxis()
        plt.legend(loc='upper left',numpoints=4)
        plt.ylabel('Approx. V-Magnitude')
        plt.xlabel('HJD')
        plt.grid(True)
        plt.xlim((t0-20,tmax + 20))
        plt.title('field ' + field + ' ' + e_id + ' ' + w_id + ' ' + a_id)
        #plt.show()
        plt.savefig('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/lc_figures/{0}.png'.format(e_id[2:-4]),bbox_inches=0)
        plt.close()
    
    if (e_data != 'NULL') & (w_data == 'NULL'): #if there is only an east lightcurve...
        t0 = e_time.min()
        tmax = e_time.max()
        e_data[:,1] -= 3.9 #puts mag on roughly v-mag system
        np.savetxt('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/east/{0}'.format(e_id),e_data)

        #plots mag vs. time and saves each lc
        plt.scatter(e_data[:,0],e_data[:,1],s=0.5,marker='.',color='red',label='east')
        #plt.scatter(w_data[:,0],w_data[:,1],s=0.5,marker='.',color='blue',label='west')
    #    plt.scatter(middle_time_list,e_med_array,s=10.0, marker='o',color='green')
    #    plt.scatter(middle_time_list,w_med_array,s=10.0, marker='o',color='black')
    #    plt.axhline(y=(w_allbins_med - offset), color='blue')
    #    plt.axhline(y=e_allbins_med, color='red', ls='dashed')
        plt.gca().invert_yaxis()
        plt.legend(loc='upper left',numpoints=4)
        plt.ylabel('Approx. V-Magnitude')
        plt.xlabel('HJD')
        plt.grid(True)
        plt.xlim((t0-20,tmax + 20))
        plt.title('field ' + field + '   ' + e_id + '   No w_id   ' + 'ABE-' + a_id)
        #plt.show()
        plt.savefig('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/lc_figures/{0}.png'.format(e_id[2:-4]),bbox_inches=0)
        plt.close()
    if (w_data != 'NULL') & (e_data == 'NULL'): #if there is only a west lightcurve...
        t0 = w_time.min()
        tmax = w_time.max()
        w_data[:,1] -= 3.9
        np.savetxt('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/west/{0}'.format(w_id),w_data)
        plt.scatter(w_data[:,0],w_data[:,1],s=0.5,marker='.',color='blue',label='west')
        plt.gca().invert_yaxis()
        plt.legend(loc='upper left',numpoints=4)
        plt.ylabel('Approx. V-Magnitude')
        plt.xlabel('HJD')
        plt.grid(True)
        plt.xlim((t0-20,tmax + 20))
        plt.title('field ' + field + '   No e_id   ' + w_id + '  ' + 'ABE-' + a_id)
        #plt.show()
        plt.savefig('/Volumes/ext_drive/work/KELT/jon_lcs/tfa/lc_figures/{0}.png'.format(w_id[2:-4]),bbox_inches=0)
        plt.close()
  
#for i in lc_files:
#    if i[5] == 'e' and i[-1] == 'g':
#        east_files.append(i)
#    if i[5] == 'w' and i[-1] == 'g':
#        west_files.append(i)
#    west_files.append(i[1]+'w_'+i[2])
#    east_files.append(i[1]+'e_lc_'+i[0]+'_mag.data')


"""
for i,filename in enumerate(east_files):
    east_data = np.genfromtxt('./be_stars_v2/{0}'.format(east_files[i]))
    west_data = np.genfromtxt('./be_stars_v2/{0}'.format(west_files[i]))
    
    east_time = east_data[:,0]
    west_time = west_data[:,0]
    if east_time.min() < west_time.min():
        t0 = east_time.min()
    else:
        t0 = west_time.min()
    if east_time.max() < west_time.max():
        tmax = west_time.max()
    else:
        tmax = east_time.max()
#makes bins of size 10 days
    bin_size = 10
    Nbins = (int(tmax + 1) - int(t0))/bin_size
    e_bin_mag = [] #initializes lists to hold mag info for each bin
    w_bin_mag = []
    w_med_list = []
    e_med_list = []
    middle_time_list = []
    for j in range(Nbins + 1):
        for k in range(east_data.shape[0]):
            if (t0 + j*bin_size) < east_data[k,0] < (t0 + (j+1)*bin_size):
                e_bin_mag.append(east_data[k,1])
        for k in range(west_data.shape[0]):
            if (t0 + j*bin_size) < west_data[k,0] < (t0 + (j+1)*bin_size):
                w_bin_mag.append(west_data[k,1])
        if len(e_bin_mag) > 30 and len(w_bin_mag) > 30:
            w_bin_med = np.median(w_bin_mag)
            e_bin_med = np.median(e_bin_mag)
            w_med_list.append(w_bin_med)
            e_med_list.append(e_bin_med)
            middle_time = t0 + (j*bin_size + (j+1)*bin_size)/2
            middle_time_list.append(middle_time)
        e_bin_mag = []
        w_bin_mag = []
    w_allbins_med = np.median(w_med_list) - 3.9
    e_allbins_med = np.median(e_med_list) - 3.9
    offset = w_allbins_med - e_allbins_med
    west_data[:,1] -= offset

    east_data[:,1] -= 3.9 #puts mag on roughly v-mag system
    west_data[:,1] -= 3.9
    outfile = east_files[i][:2] + east_files[i][3:]

    #outputs a lc with both east and west data in one file
    lc_data = np.vstack([east_data,west_data])
    np.savetxt('./combined/{0}'.format(filename),lc_data)
#    np.savetxt('./east_lcs/'.format(i),east_data)
#    np.savetxt('./west_lcs/04w_lc_023161.raw.mag',west_data)

    e_med_array=np.array(e_med_list)- 3.9
    w_med_array = np.array(w_med_list) - 3.9 - offset
    #plots mag vs. time and saves each lc
    plt.scatter(east_data[:,0],east_data[:,1],s=0.5,marker='.',color='red',label='east')
    plt.scatter(west_data[:,0],west_data[:,1],s=0.5,marker='.',color='blue',label='west')
#    plt.scatter(middle_time_list,e_med_array,s=10.0, marker='o',color='green')
#    plt.scatter(middle_time_list,w_med_array,s=10.0, marker='o',color='black')
#    plt.axhline(y=(w_allbins_med - offset), color='blue')
#    plt.axhline(y=e_allbins_med, color='red', ls='dashed')
    plt.gca().invert_yaxis()
    plt.legend(loc='upper left',numpoints=4)
    plt.ylabel('Approx. V-Magnitude')
    plt.xlabel('HJD')
    plt.grid(True)
    plt.xlim((t0-20,tmax + 20))
    plt.title('field ' + field + ' ' + e_id + ' ' + w_id + ' ' + a_id)
#    plt.show()
    plt.savefig('./lc_figures/{0}.png'.format(filename),bbox_inches=0)
    plt.close()

"""
