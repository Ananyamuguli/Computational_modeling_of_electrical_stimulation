from motif2_final import ffi
import numpy as np
from scipy.signal import find_peaks
import PyDSTool as dst
from PyDSTool import *
import pickle as pkl
import gc



#Firing rate fuction
# %%
def firing_rate_cal_func(pts, threshold=5):
    peak_indx_1 = find_peaks(pts['Vs'], height=threshold, prominence=60) #prominences makes sure, spurious spikes below 0mV or extended depolarization blocks are not taken as spikes.
    peak_indx_2 = find_peaks(pts['Vpost'], height=threshold, prominence=60)
    firing_freq_pre_val = len(peak_indx_1[0])/max(pts['t'])
    firing_freq_post_val = len(peak_indx_2[0])/max(pts['t'])
    return firing_freq_pre_val, firing_freq_post_val



#Single variable AC stim
def glu_ac_bn(run):
    extrasyn_glu = [0, 0.06, 0.2, 0.29, 0.45]
    ac_amp = np.arange(0.5,9)
    ac_freq = [0.5,1,2,5,10,15,20,25,30,50,60,100,120,140, 150]
    glu_acstim_dict = {}
    print('AC stimulation')

    for g in extrasyn_glu:
        print('Glutamate: ', g)
        freq_firing_pre_stim = [] #presynaptic neuron
        freq_firing_post_stim = [] #post_synaptic neuron firing

        for ac in ac_amp:
            print('AC amplitude: ', ac)
            freq_firing_pre = []
            freq_firing_post = []
            for f in ac_freq:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, extrasyn_glu=g, forward_syn_weight=15, neuromod='VE_sin(t)', freq_sin = f, Amp_sin = ac, stim_Vs=1, stim_Vpost=1)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        glu_acstim_dict[g] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()

    with open('glutamate_acstim_analysis_bn_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(glu_acstim_dict, f)




def glu_ac_pyr(run):
    extrasyn_glu = [0, 0.06, 0.2, 0.29, 0.45]
    ac_amp = np.arange(0.5,9)
    ac_freq = [0.5,1,2,5,10,15,20,25,30,50,60,100,120,140, 150]
    glu_acstim_dict = {}
    print('AC stimulation')

    for g in extrasyn_glu:
        print('Glutamate: ', g)
        freq_firing_pre_stim = [] #presynaptic neuron
        freq_firing_post_stim = [] #post_synaptic neuron firing

        for ac in ac_amp:
            print('AC amplitude: ', ac)
            freq_firing_pre = []
            freq_firing_post = []
            for f in ac_freq:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, extrasyn_glu=g, forward_syn_weight=15, neuromod='VE_sin(t)', freq_sin = f, Amp_sin = ac, stim_Vs=1, stim_Vpost=0)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        glu_acstim_dict[g] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()

    with open('glutamate_acstim_analysis_pyr_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(glu_acstim_dict, f)





def glu_ac_inh(run):
    extrasyn_glu = [0, 0.06, 0.2, 0.29, 0.45]
    ac_amp = np.arange(0.5,9)
    ac_freq = [0.5,1,2,5,10,15,20,25,30,50,60,100,120,140, 150]
    glu_acstim_dict = {}
    print('AC stimulation')

    for g in extrasyn_glu:
        print('Glutamate: ', g)
        freq_firing_pre_stim = [] #presynaptic neuron
        freq_firing_post_stim = [] #post_synaptic neuron firing

        for ac in ac_amp:
            print('AC amplitude: ', ac)
            freq_firing_pre = []
            freq_firing_post = []
            for f in ac_freq:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, extrasyn_glu=g, forward_syn_weight=15, neuromod='VE_sin(t)', freq_sin = f, Amp_sin = ac, stim_Vs=0, stim_Vpost=1)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        glu_acstim_dict[g] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()

    with open('glutamate_acstim_analysis_inh_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(glu_acstim_dict, f)





def Ko_ac_bn(run):
    Ko_vals = [2,4,5.8, 8, 9.9, 11]
    ac_amp = np.arange(0.5,9)
    ac_freq = [0.5,1,2,5,10,15,20,25,30,50,60,100,120,140,150]
    Ko_acstim_dict = {}
    print('AC stimulation')

    for Ko in Ko_vals:
        print('Ko: ', Ko)
        freq_firing_pre_stim = []
        freq_firing_post_stim = []

        for ac in ac_amp:
            print('AC amplitude: ', ac)
            freq_firing_pre = []
            freq_firing_post = []
            for f in ac_freq:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, Ko_rest=Ko, forward_syn_weight=15, neuromod='VE_sin(t)', freq_sin = f, Amp_sin = ac, stim_Vs=1, stim_Vpost=1)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        Ko_acstim_dict[Ko] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()

    with open('Ko_acstim_analysis_bn_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(Ko_acstim_dict, f)





def Ko_ac_pyr(run):
    Ko_vals = [2,4,5.8, 8, 9.9, 11]
    ac_amp = np.arange(0.5,9)
    ac_freq = [0.5,1,2,5,10,15,20,25,30,50,60,100,120,140,150]
    Ko_acstim_dict = {}
    print('AC stimulation')

    for Ko in Ko_vals:
        print('Ko: ', Ko)
        freq_firing_pre_stim = []
        freq_firing_post_stim = []

        for ac in ac_amp:
            print('AC amplitude: ', ac)
            freq_firing_pre = []
            freq_firing_post = []
            for f in ac_freq:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, Ko_rest=Ko, forward_syn_weight=15, neuromod='VE_sin(t)', freq_sin = f, Amp_sin = ac, stim_Vs=1, stim_Vpost=0)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        Ko_acstim_dict[Ko] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()


    with open('Ko_acstim_analysis_pyr_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(Ko_acstim_dict, f)





def Ko_ac_inh(run):
    Ko_vals = [2,4,5.8, 8, 9.9, 11]
    ac_amp = np.arange(0.5,9)
    ac_freq = [0.5,1,2,5,10,15,20,25,30,50,60,100,120,140,150]
    Ko_acstim_dict = {}
    print('AC stimulation')

    for Ko in Ko_vals:
        print('Ko: ', Ko)
        freq_firing_pre_stim = []
        freq_firing_post_stim = []

        for ac in ac_amp:
            print('AC amplitude: ', ac)
            freq_firing_pre = []
            freq_firing_post = []
            for f in ac_freq:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, Ko_rest=Ko, forward_syn_weight=15, neuromod='VE_sin(t)', freq_sin = f, Amp_sin = ac, stim_Vs=0, stim_Vpost=1)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        Ko_acstim_dict[Ko] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()


    with open('Ko_acstim_analysis_inh_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(Ko_acstim_dict, f)





def atp_ac_bn(run):
    atp_vals = [0.01, 0.51, 1.01, 1.51, 2, 5, 10]
    ac_amp = np.arange(0.5,9)
    ac_freq = [0.5,1,2,5,10,15,20,25,30,50,60,100,120,140,150]
    atp_acstim_dict = {}
    print('AC stimulation')

    for atp_v in atp_vals:
        print('ATP: ', atp_v)
        freq_firing_pre_stim = []
        freq_firing_post_stim = []

        for ac in ac_amp:
            print('AC amplitude: ', ac)
            freq_firing_pre = []
            freq_firing_post = []
            for f in ac_freq:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, atp=atp_v, forward_syn_weight=15, neuromod='VE_sin(t)', freq_sin = f, Amp_sin = ac, stim_Vs=1, stim_Vpost=1)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        atp_acstim_dict[atp_v] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()

    with open('atp_acstim_analysis_bn_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(atp_acstim_dict, f)





def atp_ac_pyr(run):
    atp_vals = [0.01, 0.51, 1.01, 1.51, 2, 5, 10]
    ac_amp = np.arange(0.5,9)
    ac_freq = [0.5,1,2,5,10,15,20,25,30,50,60,100,120,140,150]
    atp_acstim_dict = {}
    print('AC stimulation')
    for atp_v in atp_vals:
        print('ATP: ', atp_v)
        freq_firing_pre_stim = []
        freq_firing_post_stim = []
        for ac in ac_amp:
            print('AC amplitude: ', ac)
            freq_firing_pre = []
            freq_firing_post = []
            for f in ac_freq:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, atp=atp_v, forward_syn_weight=15, neuromod='VE_sin(t)', freq_sin = f, Amp_sin = ac, stim_Vs=1, stim_Vpost=0)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        atp_acstim_dict[atp_v] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()

    with open('atp_acstim_analysis_pyr_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(atp_acstim_dict, f)





def atp_ac_inh(run):
    atp_vals = [0.01, 0.51, 1.01, 1.51, 2, 5, 10]
    ac_amp = np.arange(0.5,9)
    ac_freq = [0.5,1,2,5,10,15,20,25,30,50,60,100,120,140,150]
    atp_acstim_dict = {}
    print('AC stimulation')
    for atp_v in atp_vals:
        print('ATP: ', atp_v)
        freq_firing_pre_stim = []
        freq_firing_post_stim = []
        for ac in ac_amp:
            print('AC amplitude: ', ac)
            freq_firing_pre = []
            freq_firing_post = []
            for f in ac_freq:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, atp=atp_v, forward_syn_weight=15, neuromod='VE_sin(t)', freq_sin = f, Amp_sin = ac, stim_Vs=0, stim_Vpost=1)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        atp_acstim_dict[atp_v] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()

    with open('atp_acstim_analysis_inh_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(atp_acstim_dict, f)




def syn_ac_bn(run):
    syn_vals = [0, 5, 15, 20, 30, 40, 50]
    ac_amp = np.arange(0.5,9)
    ac_freq = [0.5,1,2,5,10,15,20,25,30,50,60,100,120,140,150]
    syn_acstim_dict = {}
    print('AC stimulation')

    for syn in syn_vals:
        print('Syn: ', syn)
        freq_firing_pre_stim = [] 
        freq_firing_post_stim = [] 

        for ac in ac_amp:
            print('AC amplitude: ', ac)
            freq_firing_pre = []
            freq_firing_post = []
            for f in ac_freq:
                print('Frequency (Hz): ', f)
                _, x_, pts = ffi(Istim=2, forward_syn_weight=syn, neuromod='VE_sin(t)',  freq_sin = f, Amp_sin = ac, stim_Vs=1, stim_Vpost=1)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        syn_acstim_dict[syn] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()


    with open('syn_acstim_analysis_bn_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(syn_acstim_dict, f)




def syn_ac_pyr(run):
    syn_vals = [0, 5, 15, 20, 30, 40, 50]
    ac_amp = np.arange(0.5,9)
    ac_freq = [0.5,1,2,5,10,15,20,25,30,50,60,100,120,140,150]
    syn_acstim_dict = {}
    print('AC stimulation')

    for syn in syn_vals:
        print('Syn: ', syn)
        freq_firing_pre_stim = [] 
        freq_firing_post_stim = [] 

        for ac in ac_amp:
            print('AC amplitude: ', ac)
            freq_firing_pre = []
            freq_firing_post = []
            for f in ac_freq:
                print('Frequency (Hz): ', f)
                _, x_, pts = ffi(Istim=2, forward_syn_weight=syn, neuromod='VE_sin(t)',  freq_sin = f, Amp_sin = ac, stim_Vs=1, stim_Vpost=0)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        syn_acstim_dict[syn] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()


    with open('syn_acstim_analysis_pyr_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(syn_acstim_dict, f)





def syn_ac_inh(run):
    syn_vals = [0, 5, 15, 20, 30, 40, 50]
    ac_amp = np.arange(0.5,9)
    ac_freq = [0.5,1,2,5,10,15,20,25,30,50,60,100,120,140,150]
    syn_acstim_dict = {}
    print('AC stimulation')

    for syn in syn_vals:
        print('Syn: ', syn)
        freq_firing_pre_stim = [] 
        freq_firing_post_stim = [] 

        for ac in ac_amp:
            print('AC amplitude: ', ac)
            freq_firing_pre = []
            freq_firing_post = []
            for f in ac_freq:
                print('Frequency (Hz): ', f)
                _, x_, pts = ffi(Istim=2, forward_syn_weight=syn, neuromod='VE_sin(t)',  freq_sin = f, Amp_sin = ac, stim_Vs=0, stim_Vpost=1)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        syn_acstim_dict[syn] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()


    with open('syn_acstim_analysis_inh_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(syn_acstim_dict, f)




def mito_ac_bn(run):
    mito_vals = [25,45,50,55, 60, 75, 100, 120]
    ac_amp = np.arange(0.5,9)
    ac_freq = [0.5,1,2,5,10,15,20,25,30,50,60,100,120,140,150]
    mito_acstim_dict = {}
    print('AC stimulation')

    for mito in mito_vals:
        print('Mito: ', mito)
        freq_firing_pre_stim = [] 
        freq_firing_post_stim = [] 

        for ac in ac_amp:
            print('AC amplitude: ', ac)
            freq_firing_pre = []
            freq_firing_post = []
            for f in ac_freq:
                print('Frequency (Hz): ', f)
                _, x_, pts = ffi(Istim=2, forward_syn_weight=15, neuromod='VE_sin(t)',freq_sin = f, Amp_sin = ac, stim_Vs=1, stim_Vpost=1, del_phi_mito=mito)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        mito_acstim_dict[mito] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()


    with open('mito_acstim_analysis_bn_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(mito_acstim_dict, f)




def mito_ac_pyr(run):
    mito_vals = [25,45,50,55, 60, 75, 100, 120]
    ac_amp = np.arange(0.5,9)
    ac_freq = [0.5,1,2,5,10,15,20,25,30,50,60,100,120,140,150]
    mito_acstim_dict = {}
    print('AC stimulation')

    for mito in mito_vals:
        print('Mito: ', mito)
        freq_firing_pre_stim = [] 
        freq_firing_post_stim = [] 

        for ac in ac_amp:
            print('AC amplitude: ', ac)
            freq_firing_pre = []
            freq_firing_post = []
            for f in ac_freq:
                print('Frequency (Hz): ', f)
                _, x_, pts = ffi(Istim=2, forward_syn_weight=15, neuromod='VE_sin(t)',freq_sin = f, Amp_sin = ac, stim_Vs=1, stim_Vpost=0, del_phi_mito=mito)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        mito_acstim_dict[mito] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()


    with open('mito_acstim_analysis_pyr_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(mito_acstim_dict, f)




def mito_ac_inh(run):
    mito_vals = [25,45,50,55, 60, 75, 100, 120]
    ac_amp = np.arange(0.5,9)
    ac_freq = [0.5,1,2,5,10,15,20,25,30,50,60,100,120,140,150]
    mito_acstim_dict = {}
    print('AC stimulation')

    for mito in mito_vals:
        print('Mito: ', mito)
        freq_firing_pre_stim = [] 
        freq_firing_post_stim = [] 

        for ac in ac_amp:
            print('AC amplitude: ', ac)
            freq_firing_pre = []
            freq_firing_post = []
            for f in ac_freq:
                print('Frequency (Hz): ', f)
                _, x_, pts = ffi(Istim=2, forward_syn_weight=15, neuromod='VE_sin(t)',freq_sin = f, Amp_sin = ac, stim_Vs=0, stim_Vpost=1, del_phi_mito=mito)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        mito_acstim_dict[mito] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()


    with open('mito_acstim_analysis_inh_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(mito_acstim_dict, f)





#-----------------------------------------Single variable DBS---------------------------------
def glu_dbs_bn(run):
    extrasyn_glu = [0, 0.06, 0.2, 0.29, 0.45]
    dbs_amp = np.arange(1.5,10.5)
    dbs_f = [50,60,75,100,120,140,150]
    stim_dict = {}
    print('DBS stimulation')

    for g in extrasyn_glu:
        print('Glutamate: ', g)
        freq_firing_pre_stim = [] #presynaptic neuron
        freq_firing_post_stim = [] #post_synaptic neuron firing

        for amp in dbs_amp:
            print('DBS amplitude: ', amp)
            freq_firing_pre = []
            freq_firing_post = []
            for f in dbs_f:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, extrasyn_glu=g, forward_syn_weight=15, neuromod='dbs_monophase(t)', dbs_freq = f, Amp_dbs = amp, stim_Vs=1, stim_Vpost=1)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        stim_dict[g] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()

    with open('glutamate_dbsstim_analysis_bn_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(stim_dict, f)



def glu_dbs_pyr(run):
    extrasyn_glu = [0, 0.06, 0.2, 0.29, 0.45]
    dbs_amp = np.arange(1.5,10.5)
    dbs_f = [50,60,75,100,120,140,150]
    stim_dict = {}
    print('DBS stimulation')

    for g in extrasyn_glu:
        print('Glutamate: ', g)
        freq_firing_pre_stim = [] #presynaptic neuron
        freq_firing_post_stim = [] #post_synaptic neuron firing

        for amp in dbs_amp:
            print('DBS amplitude: ', amp)
            freq_firing_pre = []
            freq_firing_post = []
            for f in dbs_f:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, extrasyn_glu=g, forward_syn_weight=15, neuromod='dbs_monophase(t)', dbs_freq = f, Amp_dbs = amp, stim_Vs=1, stim_Vpost=0)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        stim_dict[g] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()

    with open('glutamate_dbsstim_analysis_pyr_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(stim_dict, f)



def glu_dbs_inh(run):
    extrasyn_glu = [0, 0.06, 0.2, 0.29, 0.45]
    dbs_amp = np.arange(1.5,10.5)
    dbs_f = [50,60,75,100,120,140,150]
    stim_dict = {}
    print('DBS stimulation')

    for g in extrasyn_glu:
        print('Glutamate: ', g)
        freq_firing_pre_stim = [] #presynaptic neuron
        freq_firing_post_stim = [] #post_synaptic neuron firing

        for amp in dbs_amp:
            print('DBS amplitude: ', amp)
            freq_firing_pre = []
            freq_firing_post = []
            for f in dbs_f:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, extrasyn_glu=g, forward_syn_weight=15, neuromod='dbs_monophase(t)', dbs_freq = f, Amp_dbs = amp, stim_Vs=0, stim_Vpost=1)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        stim_dict[g] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()

    with open('glutamate_dbsstim_analysis_inh_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(stim_dict, f)



def Ko_dbs_bn(run):
    Ko_vals = [2,4,5.8, 8, 9.9, 11]
    dbs_amp = np.arange(1.5,10.5)
    dbs_f = [50,60,75,100,120,140,150]
    stim_dict = {}
    print('DBS stimulation')

    for Ko in Ko_vals:
        print('Ko: ', Ko)
        freq_firing_pre_stim = [] #presynaptic neuron
        freq_firing_post_stim = [] #post_synaptic neuron firing

        for amp in dbs_amp:
            print('DBS amplitude: ', amp)
            freq_firing_pre = []
            freq_firing_post = []
            for f in dbs_f:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, Ko_rest=Ko, forward_syn_weight=15, neuromod='dbs_monophase(t)', dbs_freq = f, Amp_dbs = amp, stim_Vs=1, stim_Vpost=1)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        stim_dict[Ko] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()

    with open('Ko_dbsstim_analysis_bn_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(stim_dict, f)



def Ko_dbs_pyr(run):
    Ko_vals = [2,4,5.8, 8, 9.9, 11]
    dbs_amp = np.arange(1.5,10.5)
    dbs_f = [50,60,75,100,120,140,150]
    stim_dict = {}
    print('DBS stimulation')

    for Ko in Ko_vals:
        print('Ko: ', Ko)
        freq_firing_pre_stim = [] #presynaptic neuron
        freq_firing_post_stim = [] #post_synaptic neuron firing

        for amp in dbs_amp:
            print('DBS amplitude: ', amp)
            freq_firing_pre = []
            freq_firing_post = []
            for f in dbs_f:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, Ko_rest=Ko, forward_syn_weight=15, neuromod='dbs_monophase(t)', dbs_freq = f, Amp_dbs = amp, stim_Vs=1, stim_Vpost=0)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        stim_dict[Ko] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()

    with open('Ko_dbsstim_analysis_pyr_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(stim_dict, f)



def Ko_dbs_inh(run):
    Ko_vals = [2, 4, 5.8, 8, 9.9, 11]
    dbs_amp = np.arange(1.5,10.5)
    dbs_f = [50,60,75,100,120,140,150]
    stim_dict = {}
    print('DBS stimulation')

    for Ko in Ko_vals:
        print('Ko: ', Ko)
        freq_firing_pre_stim = [] #presynaptic neuron
        freq_firing_post_stim = [] #post_synaptic neuron firing

        for amp in dbs_amp:
            print('DBS amplitude: ', amp)
            freq_firing_pre = []
            freq_firing_post = []
            for f in dbs_f:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, Ko_rest=Ko, forward_syn_weight=15, neuromod='dbs_monophase(t)', dbs_freq = f, Amp_dbs = amp, stim_Vs=0, stim_Vpost=1)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        stim_dict[Ko] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()

    with open('Ko_dbsstim_analysis_bn_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(stim_dict, f)


def atp_dbs_bn(run):
    atp_vals = [0.01, 0.51, 1.01, 1.51, 2, 5, 10]
    dbs_amp = np.arange(1.5,10.5)
    dbs_f = [50,60,75,100,120,140,150]
    stim_dict = {}
    print('DBS stimulation')

    for atp in atp_vals:
        print('ATP: ', atp)
        freq_firing_pre_stim = [] #presynaptic neuron
        freq_firing_post_stim = [] #post_synaptic neuron firing

        for amp in dbs_amp:
            print('DBS amplitude: ', amp)
            freq_firing_pre = []
            freq_firing_post = []
            for f in dbs_f:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, atp=atp, forward_syn_weight=15, neuromod='dbs_monophase(t)', dbs_freq = f, Amp_dbs = amp, stim_Vs=1, stim_Vpost=1)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        stim_dict[atp] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()

    with open('atp_dbsstim_analysis_bn_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(stim_dict, f)



def atp_dbs_pyr(run):
    atp_vals = [0.01, 0.51, 1.01, 1.51, 2, 5, 10]
    dbs_amp = np.arange(1.5,10.5)
    dbs_f = [50,60,75,100,120,140,150]
    stim_dict = {}
    print('DBS stimulation')

    for atp in atp_vals:
        print('ATP: ', atp)
        freq_firing_pre_stim = [] #presynaptic neuron
        freq_firing_post_stim = [] #post_synaptic neuron firing

        for amp in dbs_amp:
            print('DBS amplitude: ', amp)
            freq_firing_pre = []
            freq_firing_post = []
            for f in dbs_f:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, atp=atp, forward_syn_weight=15, neuromod='dbs_monophase(t)', dbs_freq = f, Amp_dbs = amp, stim_Vs=1, stim_Vpost=0)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        stim_dict[atp] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()

    with open('atp_dbsstim_analysis_pyr_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(stim_dict, f)



def atp_dbs_inh(run):
    atp_vals = [0.01, 0.51, 1.01, 1.51, 2, 5, 10]
    dbs_amp = np.arange(1.5,10.5)
    dbs_f = [50,60,75,100,120,140,150]
    stim_dict = {}
    print('DBS stimulation')

    for atp in atp_vals:
        print('ATP: ', atp)
        freq_firing_pre_stim = [] #presynaptic neuron
        freq_firing_post_stim = [] #post_synaptic neuron firing

        for amp in dbs_amp:
            print('DBS amplitude: ', amp)
            freq_firing_pre = []
            freq_firing_post = []
            for f in dbs_f:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, atp=atp, forward_syn_weight=15, neuromod='dbs_monophase(t)', dbs_freq = f, Amp_dbs = amp, stim_Vs=0, stim_Vpost=1)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        stim_dict[atp] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()

    with open('atp_dbsstim_analysis_inh_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(stim_dict, f)



def syn_dbs_bn(run):
    syn_vals = [0, 5, 15, 20, 30, 40, 50]
    dbs_amp = np.arange(1.5,10.5)
    dbs_f = [50,60,75,100,120,140,150]
    stim_dict = {}
    print('DBS stimulation')

    for syn in syn_vals:
        print('Syn: ', syn)
        freq_firing_pre_stim = [] #presynaptic neuron
        freq_firing_post_stim = [] #post_synaptic neuron firing

        for amp in dbs_amp:
            print('DBS amplitude: ', amp)
            freq_firing_pre = []
            freq_firing_post = []
            for f in dbs_f:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, forward_syn_weight=syn, neuromod='dbs_monophase(t)', dbs_freq = f, Amp_dbs = amp, stim_Vs=1, stim_Vpost=1)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        stim_dict[syn] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()

    with open('syn_dbsstim_analysis_bn_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(stim_dict, f)



def syn_dbs_pyr(run):
    syn_vals = [0, 5, 15, 20, 30, 40, 50]
    dbs_amp = np.arange(1.5,10.5)
    dbs_f = [50,60,75,100,120,140,150]
    stim_dict = {}
    print('DBS stimulation')

    for syn in syn_vals:
        print('Syn: ', syn)
        freq_firing_pre_stim = [] #presynaptic neuron
        freq_firing_post_stim = [] #post_synaptic neuron firing

        for amp in dbs_amp:
            print('DBS amplitude: ', amp)
            freq_firing_pre = []
            freq_firing_post = []
            for f in dbs_f:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, forward_syn_weight=syn, neuromod='dbs_monophase(t)', dbs_freq = f, Amp_dbs = amp, stim_Vs=1, stim_Vpost=0)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        stim_dict[syn] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()

    with open('syn_dbsstim_analysis_pyr_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(stim_dict, f)



def syn_dbs_inh(run):
    syn_vals = [0, 5, 15, 20, 30, 40, 50]
    dbs_amp = np.arange(1.5,10.5)
    dbs_f = [50,60,75,100,120,140,150]
    stim_dict = {}
    print('DBS stimulation')

    for syn in syn_vals:
        print('Syn: ', syn)
        freq_firing_pre_stim = [] #presynaptic neuron
        freq_firing_post_stim = [] #post_synaptic neuron firing

        for amp in dbs_amp:
            print('DBS amplitude: ', amp)
            freq_firing_pre = []
            freq_firing_post = []
            for f in dbs_f:
                print('Frequency (Hz)', f)
                _, x_, pts = ffi(Istim=2, forward_syn_weight=syn, neuromod='dbs_monophase(t)', dbs_freq = f, Amp_dbs = amp, stim_Vs=0, stim_Vpost=1)
                fr_pre, fr_post = firing_rate_cal_func(pts)
                freq_firing_pre.append(fr_pre)
                freq_firing_post.append(fr_post)
                del _, x_, pts
                gc.collect()
            freq_firing_pre_stim.append(freq_firing_pre)
            freq_firing_post_stim.append(freq_firing_post)
            del freq_firing_pre, freq_firing_post
            gc.collect()

        stim_dict[syn] = (freq_firing_pre_stim, freq_firing_post_stim)
        del freq_firing_pre_stim, freq_firing_post_stim
        gc.collect()

    with open('syn_dbsstim_analysis_inh_{}.pkl'.format(run), 'wb+') as f:
        pkl.dump(stim_dict, f)


