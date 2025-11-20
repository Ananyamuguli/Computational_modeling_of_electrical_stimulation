import sys
from motif1_SBIanalysis import *


if len(sys.argv)==3:
    if int(sys.argv[1]) == 1:
        print('Glutamate analysis')
        glutamate_analysis(sys.argv[2])
        print('-------------Done---------------')
    if int(sys.argv[1]) == 2:
        print('ATP')
        atp_dependency(sys.argv[2])
        print('-------------Done---------------')
    if int(sys.argv[1]) == 3:
        print('Mito')
        mito_dependency(sys.argv[2])
        print('-------------Done---------------')
    if int(sys.argv[1]) == 4:
        print('Ko')
        Ko_dependency(sys.argv[2])
        print('-------------Done---------------')
    if int(sys.argv[1]) == 5:
        print('Synaptic weight')
        synweight_dependency(sys.argv[2])
        print('-------------Done---------------')        
else:
    print('Wrong set of inputs in the bash file')
