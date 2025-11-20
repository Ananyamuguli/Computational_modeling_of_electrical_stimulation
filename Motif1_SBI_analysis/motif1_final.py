from sre_parse import Verbose
from PyDSTool.Toolbox.phaseplane import *
import PyDSTool as dst
from PyDSTool import *
import numpy as np

def ffe(Istim1 = 2,  Ko_rest = 4, forward_syn_weight = 15, extrasyn_glu = 0.0, atp = 20, del_phi_mito=180, 
            Amp_dc=0, freq_sin = 10, Amp_sin = 0, Amp_dbs = 0, pho_dbs = 20, neuromod='0', time=10000):

    ngv_model = dst.args(name='twoneuron_system_ffe')
    ngv_model.pars = { #Fixed nernst potentials mV
                        'ECan':-20,
                        'Cm':1,

                    #Conductances in mS/cm2
                        'gNa':55,
                        'gNaP':0.15,
                        'gK':15,
                        'gCa_s':1.5,
                        'g_Can':0.025,
                        'gA':1,
                        'gKs':2,
                        'gCl':0.05,
                        'tau_Ca_s':0.24, #ms
                        

                    #a, b parameters
                        'a_can':0.0056, #([s(mM)]^-1)
                        'b_can':0.002, #/s 500s
            

                    #Ionic concentrations equations (mM)
                        'Nao_rest':120,
                        'Nai_rest':18,
                        'Ko_rest': Ko_rest,
                        'Ki_rest':110,
                        'Cao_rest':3,
                        'Cai_rest':80e-6,
                        'Cli_rest':6,
                        'Clo_rest':65,
                        'NaiA':18, #mM
            
                    #Conversion parameters
                        'beta':3,#Intracellular Volume 
                        'F':96485, #C/mol
                        'Am':922, #um2 neurons surface area
                        'vi_neuron':2160, #um3
                        'vo_neuron':720, #um3
                        'gamma_in_neuron': 0.04424, #[mM/s/uA/cm2]
                        'beta' : 3,
                        'R':8.314,
                        'T':310,
                        
                    #ATPases
                        'rho_max' : 15, #uA/cm2
                        'epsilon_kmax' : 0.25, #/s ~ 4s
                        'Ukcc2' : 0.3, #mM/s
                        'Unkcc1':0.1, #mM/s
                        'glia_uptake' : 12, #uA/cm2

                    #Receptors, synaptic, stroke
                        'Mg2':0.01,
                        'Enmda':0, 
                        'Eampa' : 0,  #reversal voltage
                        'gnmdar': 0.02, #0.03, changed for exp2
                        'gampa' : 0.01, #mS/cm2
                        'alpha_ampar' : 1.1e3,  #synaptic
                        'beta_ampar': 190,  #synaptic
                        'alpha_nmdar' : 72, #mM/s #synaptic
                        'beta_nmdar' : 6.6, #s #synaptic
                        'T_max': 20,  #mM maximum glutamate relase from presynaptic terminal. Liable for adjusment. Can range from 1-15mM
                        'V_p': 0, #mV
                        'K_p': 5,  #mV
                        'tau_g':0.025,

                    #neuromodulation, 
                        'Amp_dc':Amp_dc, 
                        'freq_sin':freq_sin, 
                        'Amp_sin':Amp_sin, 
                        'pho_dbs': pho_dbs, #ms, time period
                        'delta_D':0.6, #ms, pulse width
                        'i_dbs':Amp_dbs, 



                    #network parameters
                        'glu_extrasyn':extrasyn_glu, #mM
                        'glu_syn_max':0.5,
                        'Istim1':Istim1,
                        'forward_syn_weight':forward_syn_weight,
                        'atp':atp, #uM #saturating state atp
                        'min_atp':0.1, #uM
                        'atp_k':2, #uM
                        'del_phi_shift_mf': 85, #mV,
                        'del_phi_mito' : del_phi_mito, # change this instead #mV mitochondria inner membrane potential. #Reduces with increased ROS -- Mitochondrial dysfunction.
                        'fca' : 0.01,
    }



    ##Neuromodulation paradigms
    neuromodulation_paradigms = {   'VE_dc': (['t'], 'Amp_dc'),
                                    'VE_sin': (['t'], 'Amp_sin * sin(2*3.14*freq_sin*t/1000)'), 
                                    'dbs':(['t'], 'i_dbs*heav(sin(2*3.14*t/pho_dbs))*(1-heav(sin(2*3.14*(t+delta_D)/pho_dbs)))'),
                                    'Iext': (['t'], neuromod)
        }

    #ATP factor
    metabolic_equations = { 'mito_dyfn_factor' :([''], '1/(1+exp(-(del_phi_mito - del_phi_shift_mf)/12))') ,
                            'energy_factor': ([''], '(1/(1+ atp_k /(min_atp + atp * mito_dyfn_factor())))'), 
                             
    }

    #Regular spiking neurons. RS.
    Na_functions_RS = { 'alpha_m_Na_RS': (['Vs'], '-0.1*(Vs+31)/(exp(-(Vs+31)/10)-1)'),
                    'beta_m_Na_RS' : (['Vs'], '4*exp(-(Vs+56)/18)'),
                    'alpha_h_Na_RS': (['Vs'], '0.07*exp(-(Vs+47)/20)'),
                    'beta_h_Na_RS' : (['Vs'], '1/(1+exp(-(Vs+17)/10))')}

    Kdr_functions_RS = {'alpha_n_Kdr_RS': (['Vs'], '0.01*(Vs+34)/(1-exp(-(Vs+34)/10))'),
                    'beta_n_Kdr_RS': (['Vs'], '0.125*exp(-(Vs+44)/80)')}


    CaL_s_functions_RS = {'m_inf_CaL_RS': (['V'], '1/(1+exp(-(V+20)/9))')}


    CaN_functions_RS = {'m_inf_CaN_RS': (['Cai'], 'a_can*Cai/(a_can*Cai + b_can)'),
                    'tau_m_CaN_RS': (['Cai'], '1/(a_can*Cai + b_can)')}



    #Ionic relations.
    Ionic_equations = {'Nao': (['Nain'], 'Nao_rest - beta*(Nain-Nai_rest)'),
                    'Ki': (['Clin', 'Nain'], 'Ki_rest + (Nai_rest - Nain) - (Cli_rest - Clin)'),
                    'Cao': (['Cai'], 'Cao_rest - beta * (Cai - Cai_rest)'),
                    'Clo': (['Clin'], 'Clo_rest - beta * (Clin - Cli_rest)')}


    #Nernst potentials
    Nernst_eqns = { 'ENa': (['Nain', 'Naout'], '26.69 * log(Naout/Nain)'),
                    'EK': (['Kin', 'Kout'], '26.69 * log(Kout/Kin)'),
                    'ECl': (['Clin', 'Clout'], '26.69 * log(Clin/Clout)'),
                    'ECa': (['Cain', 'Caout'], '26.69 * 0.5* log(Caout/Cain)'),
                }


    #Ion pump transporters.
    pumps_transporters = { 'INaKpump' : ( ['Nain', 'Kout'], 'energy_factor() * rho_max*((1 + exp((28-Nain)/3))**-1) *((1 + exp(5.5-Kout))**-1)'),
            'Iglialpump' : (['Kout'], 'energy_factor() * rho_max*((1 + exp((30-NaiA)/4))**-1)*((1 + exp(5.5-Kout)/3)**-1)*(3**-1)'),
            'Iglia' : (['Kout'], 'energy_factor() * glia_uptake * (1 + exp((25-Kout)/2.5))**-1'),
            'Idiff' : (['Kout'],  'energy_factor()  * epsilon_kmax * (Kout - Ko_rest)'),
            'Ikcc2' : (['Kin', 'Kout', 'Clin','Clout'],  'Ukcc2 * log( (Kin * Clin)/(Kout * Clout) )'),
            'Inkcc1': (['Kin', 'Kout', 'Clin', 'Clout', 'Nain', 'Naout'], 'Unkcc1 * (1/(1 + exp(16-Kout))) * (log((Kin * Clin)/(Kout * Clout)) + log((Nain * Clin)/(Naout * Clout)) )')
        }



    receptor_functions = { 'Mg2_block':(['V'],'(1/(1 + exp(-0.062*(V))* Mg2/3.57))' ),
                        'glutamate_rel_presyn': (['V'], 'T_max/(1+exp(-(V - V_p)/K_p))' ), 
    }

    stimulus_current = {'Isyn_exc_forward': (['nmdar', 'ampar', 'Vpost'], 'gnmdar * nmdar * Mg2_block(Vpost) * (Vpost - Enmda) + gampa * ampar * (Vpost - Eampa)'),
                        #'Is1':(['t'], 'if(t<5000, Istim1 ,0)'),
                        'Is1':(['t'], 'Istim1')
    }




    #Main ODE Functions
    membrane_voltage_equations = {  #Presynaptic neuron
                                    'Vs': '(Is1(t) + Iext(t) -\
                                        INaKpump(Nai, Ko) - \
                                        gCa_s*(pow(m_inf_CaL_RS(Vs),2)) * (Vs  - ECa(Cai_s, Cao(Cai_s))) - \
                                        gNa*(pow(m,3))*h*(Vs  - ENa(Nai, Nao(Nai))) - \
                                        gK*pow(n, 4)*(Vs  - EK(Ki(Cli, Nai), Ko)) - \
                                        g_Can*pow(m_CaN,2)*(Vs  - ECan) - \
                                        gnmdar * nmdar_extrasyn * Mg2_block(Vs) *(Vs  - Enmda ) -\
                                        gCl*(Vs  - ECl(Cli, Clo(Cli) ) ))/Cm', 


                                #Post synaptic neuron
                                'Vpost': '(Iext(t) - forward_syn_weight * Isyn_exc_forward(nmdar_syn, ampar_syn, Vpost)-\
                                        INaKpump(Nai_post, Ko_post) - \
                                        gCa_s*(pow(m_inf_CaL_RS(Vpost),2)) * (Vpost  - ECa(Cai_post, Cao(Cai_post))) - \
                                        gNa*pow(m_post,3)*h_post *(Vpost  - ENa(Nai_post, Nao(Nai_post))) - \
                                        gK*pow(n_post,4)*(Vpost  - EK(Ki(Cli_post, Nai_post), Ko_post)) - \
                                        gnmdar * nmdar_extrasyn * Mg2_block(Vpost) *(Vpost  - Enmda ) -\
                                        g_Can*pow(m_CaN,2)*(Vpost  - ECan) - \
                                        gCl*(Vpost  - ECl(Cli_post, Clo(Cli_post))))/Cm'
                                        }



    gate_ionic_equations = {#Presynaptic neuron
                        'Cai_s': '-fca* gamma_in_neuron * (gCa_s * pow(m_inf_CaL_RS(Vs),2)* (Vs  - ECa(Cai_s, Cao(Cai_s))) + gnmdar * nmdar_extrasyn * Mg2_block(Vs) *(Vs  - Enmda )) -\
                            energy_factor() * Cai_s/tau_Ca_s',

                        'm': 'alpha_m_Na_RS(Vs) * (1-m) - beta_m_Na_RS(Vs) * m',

                        'h': 'alpha_h_Na_RS(Vs)*(1-h) - beta_h_Na_RS(Vs)*h',

                        'n': 'alpha_n_Kdr_RS(Vs)*(1-n) - beta_n_Kdr_RS(Vs)*n',

                        'm_CaN': '(m_inf_CaN_RS(Cai_s)-m_CaN)/tau_m_CaN_RS(Cai_s)',

                        'Nai': '-gamma_in_neuron * gNa*pow(m,3)*h*(Vs -ENa(Nai, Nao(Nai))) - \
                            3 * INaKpump(Nai, Ko)* gamma_in_neuron  -\
                            Inkcc1(Ki(Cli, Nai), Ko, Cli, Clo(Cli), Nai, Nao(Nai))-\
                            gnmdar * nmdar_extrasyn * Mg2_block(Vs) *(Vs - Enmda) * gamma_in_neuron',

                        'Ko': 'gamma_in_neuron * gK * pow(n,4) * (Vs - EK(Ki(Cli, Nai), Ko)) - \
                            2 * beta * gamma_in_neuron * INaKpump(Nai, Ko) - \
                            Idiff(Ko) -\
                            2 * gamma_in_neuron * Iglialpump(Ko) - \
                            gamma_in_neuron * Iglia(Ko) + beta * (Ikcc2(Ki(Cli, Nai), Ko, Cli, Clo(Cli)) + Inkcc1(Ki(Cli, Nai), Ko, Cli, Clo(Cli), Nai, Nao(Nai)))',

                        'Cli': 'gamma_in_neuron * gCl * (Vs - ECl(Cli,Clo(Cli)))- Ikcc2( Ki(Cli, Nai), Ko, Cli, Clo(Cli)) - 2*Inkcc1(Ki(Cli, Nai), Ko, Cli, Clo(Cli), Nai, Nao(Nai))',

                        'glu_syn' :  '-glu_syn/tau_g +  glutamate_rel_presyn(Vs)',

                        'nmdar_extrasyn': 'glu_extrasyn * alpha_nmdar * (1- nmdar_extrasyn) - beta_nmdar * nmdar_extrasyn', 

                        'nmdar_syn':  'glu_syn * alpha_nmdar*(1- nmdar_syn)/glu_syn_max - beta_nmdar * nmdar_syn',

                        'ampar_syn': 'glu_syn * alpha_ampar*(1-ampar_syn)/glu_syn_max - beta_ampar * ampar_syn', 

                        'Cai_post':    '-fca* gamma_in_neuron * (gCa_s * pow(m_inf_CaL_RS(Vpost),2)* (Vpost  - ECa(Cai_post, Cao(Cai_post))) + gnmdar * nmdar_extrasyn * Mg2_block(Vpost) *(Vpost  - Enmda )) -\
                                        energy_factor() * Cai_post/tau_Ca_s', 

                        'm_post': 'alpha_m_Na_RS(Vpost)*(1-m_post) - beta_m_Na_RS(Vpost) * m_post',

                        'h_post': 'alpha_h_Na_RS(Vpost)*(1-h_post) - beta_h_Na_RS(Vpost) * h_post',

                        'n_post': 'alpha_n_Kdr_RS(Vpost)*(1-n_post) - beta_n_Kdr_RS(Vpost) * n_post',

                        'm_CaN_post': '(m_inf_CaN_RS(Cai_post)-m_CaN_post)/tau_m_CaN_RS(Cai_post)',

                        'Nai_post': '-gamma_in_neuron * gNa* pow(m_post, 3)* h_post * (Vpost - ENa(Nai_post, Nao(Nai_post))) -\
                                    3 * INaKpump(Nai_post, Ko_post)*gamma_in_neuron - \
                                    Inkcc1(Ki(Cli_post, Nai_post), Ko_post, Cli_post, Clo(Cli_post), Nai_post, Nao(Nai_post)) -\
                                    gnmdar * 3 * nmdar_extrasyn * Mg2_block(Vpost) *(Vpost - Enmda ) * gamma_in_neuron',

                        'Ko_post': 'gamma_in_neuron * gK * pow(n_post,4) * (Vpost - EK(Ki(Cli_post, Nai_post), Ko_post)) - \
                                    2 * beta * gamma_in_neuron * INaKpump(Nai_post, Ko_post) -\
                                    Idiff(Ko_post) -\
                                    2 * gamma_in_neuron * Iglialpump(Ko_post) - gamma_in_neuron * Iglia(Ko_post) + beta * (Ikcc2(Ki(Cli_post, Nai_post), Ko_post, Cli_post, Clo(Cli_post)) + \
                                    Inkcc1(Ki(Cli_post, Nai_post), Ko_post, Cli_post, Clo(Cli_post), Nai_post, Nao(Nai_post)))',


                        'Cli_post': 'gamma_in_neuron * gCl * (Vpost - ECl(Cli_post , Clo(Cli_post))) - Ikcc2( Ki(Cli_post, Nai_post), Ko_post, Cli_post, Clo(Cli_post)) - 2*Inkcc1(Ki(Cli_post, Nai_post), Ko_post, Cli_post, Clo(Cli_post), Nai_post, Nao(Nai_post))',

                        }




    ngv_model.fnspecs = { **metabolic_equations, **Na_functions_RS,  **Kdr_functions_RS, **CaL_s_functions_RS, **CaN_functions_RS, **Ionic_equations, **pumps_transporters, 
                        **Nernst_eqns, **receptor_functions,  **neuromodulation_paradigms, **stimulus_current}

    ngv_model.varspecs = {**membrane_voltage_equations,  **gate_ionic_equations}


    ngv_model.ics =  {'Vs':-75,'Cai_s': 80e-5,  'm':0.005, 'h':0.99, 'n':0.005, 'm_CaN':0.1, 
                    'Nai':18, 'Ko':4, 'Cli':6, 'glu_syn':0.0001,  'nmdar_syn':0.01, 'ampar_syn':0.01,'Vpost': -75,  'm_post':0.005, 'h_post':0.15, 'n_post':0.5,
                    'Nai_post': 18, 'Ko_post':4, 'Cli_post':6,'nmdar_extrasyn':0.00001, 'Cai_post': 80e-5, 'm_CaN_post':0.1}
    ngv_model.tdata = [0, time]
    ngv_model.algparams= {'max_pts':100000000}



    ode1 = dst.Generator.Radau_ODEsystem(ngv_model)
    ode1.xdomain =  {'Vs':[-100,80],'Cai_s': [1e-5, 4],  'm':[0,1], 'h':[0,1], 'n':[0,1], 'm_CaN': [0,1], 
                    'Nai':[10,40], 'Ko':[2,35], 'Cli':[2,20],  'glu_syn':[1e-4,10], 
                     'nmdar_syn':[0,1], 'ampar_syn':[0,1], 'Vpost': [-100,80],  'm_post':[0,1], 'h_post':[0,1], 'n_post':[0,1],
                    'Nai_post': [10,40], 'Ko_post':[2,35], 'Cli_post':[2,20], 'nmdar_extrasyn':[0.0,1.0],  'Cai_post':[1e-5, 4], 'm_CaN_post': [0,1]}


    traj = ode1.compute('twoneuron_system_ffe')
    pts  = traj.sample(dt=0.01)
    return ode1, traj, pts