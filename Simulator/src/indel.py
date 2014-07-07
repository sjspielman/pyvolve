# File contains all kindsa indel stuff which has yet to be formally placed.
import misc
import numpy as np
from scipy import special


def indelDist_RZ():
    ''' Use Riemann Zeta function to get indel length distributions for a given set of indel model params.
        THIS SHOULD ULTIMATELY BE DONE IN SOME KIND OF SETUP FILE, but for now here is ok.
    '''
    max_i = model.indelParams['maxIns'] # default should be like 50
    max_d = model.indelParams['maxDel'] # default should be like 50
    alpha_i = model.indelParams['insAlpha'] # default 1.75
    alpha_d = model.indelParams['delAlpha'] # default 1.75
    
    model.indelParams['insDist'] = np.zeros(max_i)
    model.indelParams['delDist'] = np.zeros(max_d)
    
    for x in range(1, max_i):
        model.indelParams['insDist'][x] = (float(x)**(alpha_i*-1.)) / special.zeta(x, 1)
    for x in range (1, max_d):
        model.indelParams['delDist'][x] = (float(x)**(alpha_d*-1.)) / special.zeta(x, 1)

