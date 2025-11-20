import numpy as np
from scipy.signal import find_peaks
import PyDSTool as dst
from PyDSTool import *
import pickle as pkl
from motif2_final import ffi



#Fundamental functions for the analysis
#Firing rate fuctions
def firing_rate_cal_func(pts, threshold=5): #Make changes to the prominence and threshold. Thresh was 10 now 5
    peak_indx_1 = find_peaks(pts['Vs'], height=threshold, prominence=60) # was 70, prominences makes sure, spurious spikes below 0mV or extended depolarization blocks are not taken as spikes.
    peak_indx_2 = find_peaks(pts['Vpost'], height=threshold, prominence=60)
    firing_freq_pre_val = len(peak_indx_1[0])/max(pts['t'])
    firing_freq_post_val = len(peak_indx_2[0])/max(pts['t'])
    return firing_freq_pre_val, firing_freq_post_val


#----------------------------------------------------------------**-----------------------------------------------------------
#FIRING RATE ANALYSIS No stim
#Glutamate excitotoxicity - dependence of glutamate with firing rate. - firing rate
def glutamate_analysis(run):
    extrasyn_glu = np.arange(0,0.61, 0.01)
    freq_firing_pre = [] #presynaptic neuron
    freq_firing_post = [] #post_synaptic neuron firing
    for g in extrasyn_glu:
        print(g)
        _, x_, pts = ffi(Istim=2, extrasyn_glu=g, forward_syn_weight=15)
        fr_pre, fr_post = firing_rate_cal_func(pts)
        freq_firing_pre.append(fr_pre)
        freq_firing_post.append(fr_post)
        del _,x_, pts
        gc.collect()

    fr_dict = {
        'fr_pre':freq_firing_pre, 
        'fr_post':freq_firing_post
    }
    with open('glutamate_analysis_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(fr_dict, f)
    del freq_firing_pre, freq_firing_post



#atp dependency firing rate
def atp_dependency(run):
    atp_levels1 = np.arange(0.01, 15, 0.5)
    atp_levels2 = np.arange(atp_levels1[-1]+0.5 , 35, 1)

    freq_firing_pre = []
    freq_firing_post = []
    for atp in atp_levels1:
        print(atp)
        _, x_, pts = ffi(Istim=2, atp=atp, forward_syn_weight=15)
        fr_pre, fr_post = firing_rate_cal_func(pts)
        freq_firing_pre.append(fr_pre)
        freq_firing_post.append(fr_post)
        del _, x_, pts
        gc.collect()

    for atp in atp_levels2:
        print(atp)
        _, x_, pts = ffi(Istim=2, atp=atp, forward_syn_weight=15)
        fr_pre, fr_post = firing_rate_cal_func(pts)
        freq_firing_pre.append(fr_pre)
        freq_firing_post.append(fr_post)
        del _, x_, pts
        gc.collect()

    fr_dict = {
        'fr_pre':freq_firing_pre, 
        'fr_post':freq_firing_post
    }
    with open('atp_dependency_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(fr_dict, f)
    del freq_firing_pre, freq_firing_post




#Mitochondrial dependency firing rate
def mito_dependency(run):
    mito_levels = np.arange(0, 180, 5)
    freq_firing_pre = []
    freq_firing_post = []
    for mito in mito_levels:
        print(mito)
        _, x_, pts = ffi(Istim=2, del_phi_mito=mito, forward_syn_weight=15)
        fr_pre, fr_post = firing_rate_cal_func(pts)
        freq_firing_pre.append(fr_pre)
        freq_firing_post.append(fr_post)
        del _,x_, pts
        gc.collect()

    fr_dict = {
        'fr_pre':freq_firing_pre, 
        'fr_post':freq_firing_post
    }
    with open('mito_dependency_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(fr_dict, f)
    del freq_firing_pre, freq_firing_post 



#Extracellular potassium
def Ko_dependency(run):
    Ko_levels = np.arange(2,12, 0.25)
    freq_firing_pre = []
    freq_firing_post = []
    for Ko in Ko_levels:
        print(Ko)
        _, x_, pts = ffi(Istim=2, Ko_rest = Ko, forward_syn_weight=15)
        fr_pre, fr_post = firing_rate_cal_func(pts)
        freq_firing_pre.append(fr_pre)
        freq_firing_post.append(fr_post)
        del _,x_, pts
        gc.collect()

    fr_dict = {
        'fr_pre':freq_firing_pre, 
        'fr_post':freq_firing_post
    }
    with open('Ko_dependency_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(fr_dict, f)      
    del freq_firing_pre, freq_firing_post







#synaptic connectivity
def synweight_dependency(run):
    syn_weight = np.arange(0,80,2) #Dimensionless
    freq_firing_pre = []
    freq_firing_post = []
    for syn in syn_weight:
        print(syn)
        _,x_, pts = ffi(Istim=2, forward_syn_weight=syn)
        fr_pre, fr_post = firing_rate_cal_func(pts)
        freq_firing_pre.append(fr_pre)
        freq_firing_post.append(fr_post)
        del _,x_, pts
        gc.collect()

    fr_dict = {
        'fr_pre':freq_firing_pre, 
        'fr_post':freq_firing_post
    }
    with open('synweight_dependency_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(fr_dict, f)
    del freq_firing_pre, freq_firing_post


