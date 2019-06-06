#author: Tiffany Zhao
#version: 7/1/17
import numpy as np
from matplotlib import pyplot as plt



import os
import random
from loadspec_spSpec import spSpec_sdss, spec_norm, spec_rebin

from sklearn.decomposition import NMF, PCA, IncrementalPCA
import h5py

from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=8, usetex=False)

#import qso catalog_path


filename='PCA_spectra_data.hdf5'
with h5py.File(filename,'r') as hf:
    data_spectra=hf["spectra"][:,:]
    wave_value=hf["wave"][:]

sub_data=data_spectra
wave = wave_value
print data_spectra.shape
print wave_value


# PCA


#n_components<no_spectra
n_components = 30
#first attempted NMF, too slow, so use IPCA
#nmf = NMF(n_components,max_iter=100)
#nmf.fit(sub_data)
#nmf_comp = nmf.components_

sub_data=sub_data-np.mean(sub_data)
#n_components= kept number of components
#n = number of components used out of n_components
ipca=IncrementalPCA(n_components=n_components)
ipca.fit(sub_data)
ipca_comp=ipca.components_
ipca_comp=ipca_comp.T

#n=30
#coeff = np.dot(sub_data[0,:],ipca_comp[:,0:n])
#new=np.dot(ipca_comp[:,0:n], coeff[0:n])
#plt.plot(sub_data[0,:],color='black')
#plt.plot(new,color='red')

#plt.show()
#print coeff.shape

#coeff=ipca.transform(sub_data)
#mean_spectra=np.mean(sub_data,axis=0)
#coeff = np.dot(nmf_comp, sub_data[1,:])
#define the fitted continuum:
x=0
n=18
coeff = np.dot(sub_data[x,:],ipca_comp[:,0:n])
flux = sub_data[x,:]
fit = np.dot(ipca_comp[:,0:n], coeff[0:n])
norm_flux=flux/fit
print flux
print fit
print flux/fit
plt.plot(wave,norm_flux)
plt.plot(wave,wave*0.0+1.0,color='red')
plt.xlim(np.min(wave),np.max(wave))
plt.ylim(-10,10)
plt.show()

#------------------------------------------------------------
for x in range(10):
    fig = plt.figure(figsize=(5, 5))
    fig.subplots_adjust(hspace=0, top=0.95, bottom=0.1, left=0.12, right=0.93)
    fig.suptitle('IPCA Reconstruction of Spectrum ' + str(x+1))
    for i, n in enumerate([ 4,12,18,20]):
        print x,n
        ax = fig.add_subplot(411 + i)

        coeff = np.dot(sub_data[x,:],ipca_comp[:,0:n])
        flux = sub_data[x,:]
        fit=np.dot(ipca_comp[:,0:n], coeff[0:n])
        #norm_flux=flux/fit
        print flux
        print fit
        print flux/fit


        ax.plot(wave, flux, '-', c='black')
        ax.plot(wave, fit, c='red')
        #ax.plot(wave,norm_flux,c='red')

        #ax.set_ylim([-5,5])




        if i < 3:
            ax.xaxis.set_major_formatter(plt.NullFormatter())

        ax.set_ylabel('flux')
        ax.set_ylim([-2,2])

        if n == 0:
            text = "mean"
        elif n == 1:
            text = "1 component\n"
        else:
            text = "%i components\n" % n

        ax.text(0.02, 0.93, text, ha='left', va='top', transform=ax.transAxes)

    fig.axes[-1].set_xlabel(r'${\rm wavelength\ (\AA)}$')
    #plt.pause(5)
    plt.show()
