import sys
from motif2_stim_analysis import *


if len(sys.argv)==3:
    if int(sys.argv[1]) == 1:
        print('AC stim analysis with glutamate: PYR-INH')
        glu_ac_bn(sys.argv[2])
        print('-------------Done---------------')
    if int(sys.argv[1]) == 2:
        print('AC stim analysis with glutamate: Stimulating only PYR')
        glu_ac_pyr(sys.argv[2])
        print('-------------Done-----------')
    if int(sys.argv[1]) == 3:
        print('AC stim analysis with glutamate: Stimulating only INH')
        glu_ac_inh(sys.argv[2])
        print('-------------Done-----------')
    if int(sys.argv[1]) == 4:
        print('AC stim analysis with Ko: Stimulating both neurons')
        Ko_ac_bn(sys.argv[2])
        print('-------------Done---------------')
    if int(sys.argv[1]) == 5:
        print('AC stim analysis with Ko: Stimulating only PYR')
        Ko_ac_pyr(sys.argv[2])
        print('-------------Done-----------')
    if int(sys.argv[1]) == 6:
        print('AC stim analysis with Ko: Stimulating only INH')
        Ko_ac_inh(sys.argv[2])
        print('-------------Done-----------')
    if int(sys.argv[1]) == 7:
        print('AC stim analysis with atp: Stimulating both neurons')
        atp_ac_bn(sys.argv[2])
        print('-------------Done---------------')
    if int(sys.argv[1]) == 8:
        print('AC stim analysis with atp: Stimulating only PYR')
        atp_ac_pyr(sys.argv[2])
        print('-------------Done-----------')
    if int(sys.argv[1]) == 9:
        print('AC stim analysis with atp: Stimulating only INH')
        atp_ac_inh(sys.argv[2])
        print('-------------Done-----------')  
    if int(sys.argv[1]) == 10:
        print('AC stim analysis with syn: Stimulating both neurons')
        syn_ac_bn(sys.argv[2])
        print('-------------Done---------------')
    if int(sys.argv[1]) == 11:
        print('AC stim analysis with syn: Stimulating only PYR')
        syn_ac_pyr(sys.argv[2])
        print('-------------Done-----------')
    if int(sys.argv[1]) == 12:
        print('AC stim analysis with syn: Stimulating only INH')
        syn_ac_inh(sys.argv[2])
        print('-------------Done-----------')
    if int(sys.argv[1]) == 13:
        print('AC stim analysis with mito: Stimulating both neurons')
        mito_ac_bn(sys.argv[2])
        print('-------------Done---------------')
    if int(sys.argv[1]) == 14:
        print('AC stim analysis with mito: Stimulating only PYR')
        mito_ac_pyr(sys.argv[2])
        print('-------------Done-----------')
    if int(sys.argv[1]) == 15:
        print('AC stim analysis with mito: Stimulating only INH')
        mito_ac_inh(sys.argv[2])
        print('-------------Done-----------')
    if int(sys.argv[1]) == 16:
        print('DBS stim analysis with glutamate: Stimulating both neurons')
        glu_dbs_bn(sys.argv[2])
        print('-------------Done---------------')
    if int(sys.argv[1]) == 17:
        print('DBS stim analysis with glutamate: Stimulating only PYR')
        glu_dbs_pyr(sys.argv[2])
        print('-------------Done-----------')
    if int(sys.argv[1]) == 18:
        print('DBS stim analysis with glutamate: Stimulating only INH')
        glu_dbs_inh(sys.argv[2])
        print('-------------Done-----------')
    if int(sys.argv[1]) == 19:
        print('DBS stim analysis with Ko: Stimulating both neurons')
        Ko_dbs_bn(sys.argv[2])
        print('-------------Done---------------')
    if int(sys.argv[1]) == 20:
        print('DBS stim analysis with Ko: Stimulating only PYR')
        Ko_dbs_pyr(sys.argv[2])
        print('-------------Done-----------')
    if int(sys.argv[1]) == 21:
        print('DBS stim analysis with Ko: Stimulating only INH')
        Ko_dbs_inh(sys.argv[2])
        print('-------------Done-----------') 
    if int(sys.argv[1]) == 22:
        print('DBS stim analysis with atp: Stimulating both neurons')
        atp_dbs_bn(sys.argv[2])
        print('-------------Done---------------')
    if int(sys.argv[1]) == 23:
        print('DBS stim analysis with atp: Stimulating only PYR')
        atp_dbs_pyr(sys.argv[2])
        print('-------------Done-----------')
    if int(sys.argv[1]) == 24:
        print('DBS stim analysis with atp: Stimulating only INH')
        atp_dbs_inh(sys.argv[2])
        print('-------------Done-----------') 
else:
    print('Wrong number of inputs')
