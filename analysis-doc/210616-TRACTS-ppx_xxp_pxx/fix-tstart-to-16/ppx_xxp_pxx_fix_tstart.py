import numpy
import scipy


def ppx_xxp_pxx(*params):
    """a simple model in which populations Eu and NAT arrive discretely at first generation, Af at a subsequent generation, 
    and Eu discretely yet later. 
    . If a time is not integer, the migration 
    is divided between neighboring times proportional to the non-integer time fraction. 
    We'll assume population 3 
    replaces after the replacement from population 1 and 2 if they arrive at same generation. 
    The first generation is equally composed of pops 1 and 2
    order is ['CEU','NAT','YRI']"""
    (init_Eu, tstart, afam_prop, afam_time, nuEu_prop, nuEu_time) = params[0]

    tstart *= 100
    afam_time *= 100
    nuEu_time *= 100
    # print "times ",tstart,afam_time,nuEu_time

    if afam_time > tstart or nuEu_time > tstart or nuEu_time < 0 or afam_time < 0:
        # that should be caught by constraint. Return empty matrix
        gen = int(numpy.ceil(max(tstart, 0)))+1
        mig = numpy.zeros((gen+1, 3))
        return mig

    #	print "error, afam_time>tstart:", afam_time,">", tstart
    #	raise()

        #print (a,b)

    gen = int(numpy.ceil(tstart))+1
    frac = gen-tstart-1
    mig = numpy.zeros((gen+1, 3))

    init_Nat = 1-init_Eu

    # replace a fraction at second generation to ensure a continuous model distribution with generation
    mig[-1, :] = numpy.array([init_Eu, init_Nat, 0])
    # contents of the second generation

    interEu = frac*init_Eu
    interNat = frac*init_Nat
    mig[-2, :] = numpy.array([interEu, interNat, 0])

    afgen = int(numpy.ceil(afam_time))+1
    fracAf = afgen-afam_time-1

    # we want the total african proportiuon replaced to be afam_prop. We therefore add a fraction
    # f at generation gen-1, and (afam_prop-f)/(1-f) at generation gen.

    mig[afgen-1, 2] = fracAf*afam_prop
    mig[afgen, 2] = (afam_prop-fracAf*afam_prop)/(1-fracAf*afam_prop)

    # finally add the European continuous migration

    gen = int(numpy.ceil(nuEu_time))
    frac = gen-nuEu_time
    # replacement from pop 3 occurs after replacement from pop 1,2, and may span 2 generations
    mig[gen-1, 0] = frac*nuEu_prop
    mig[gen, 0] = (nuEu_prop-frac*nuEu_prop)/(1-frac*nuEu_prop)

    return mig


def outofbounds_ppx_xxp_pxx(*params):
    # constraint function evaluating below zero when constraints not satisfied
    # print "in constraint  outofbounds_211_cont_unif_params"
    ret = 1
    (init_Eu, tstart, afam_prop, afam_time, nuEu_prop, nuEu_time) = params[0]
    ret = min(1-init_Eu, 1-afam_prop, 1-nuEu_prop)
    ret = min(ret, init_Eu, afam_prop, nuEu_prop)

    # print "ret0",  ret

    # pedestrian way of testing for all possible issues
    func = ppx_xxp_pxx
    mig = func(params[0])
    totmig = mig.sum(axis=1)
    # print  "ret1=",ret

    if init_Eu > 1 or afam_prop > 1 or nuEu_prop > 1:
        print("Pulse greater than 1")
    if init_Eu < 0 or afam_prop < 0 or nuEu_prop < 0:
        print("Pulse less than 0")
    # print  "ret2 ",ret
    ret = min(ret, -abs(totmig[-1]-1)+1e-8)
    ret = min(ret, -totmig[0], -totmig[1])
    # print "ret3 " , ret
    ret = min(ret, min(1-totmig), min(totmig))
    # print "ret4 " , ret

    # print "times ",afam_time,tstart
    ret = min(ret, tstart-afam_time)
    ret = min(ret, tstart-nuEu_time)
    ret = min(ret, afam_time)
    ret = min(ret, nuEu_time)
    ret = min(ret, tstart)

    # print "ret5 " , ret
    if abs(totmig[-1]-1) > 1e-8:
        print(mig)
        print("founding migration should sum up to 1. Now:")
        # print mig[-1,:],"sum up to ",self.totmig[-1])

    if totmig[0] > 1e-10:
        print("migrants at last generation should be removed from sample!")
        #print("currently", self.totmig[0])

    if totmig[1] > 1e-10:
        print("migrants at penultimate generation should be removed from sample!")
        #print("currently", self.totmig[1])

    if ((totmig > 1).any() or (mig < 0).any()):
        print("migration rates should be between 0 and 1")
    # print "constraint ",ret
    return ret

# We don't have to calculate all the tract length distributions to have the
# global ancestry proportion right. Here we define a function that
# automatically adjusts the migration rates to have the proper global ancestry
# proportion, saving a lot of optimization time!
# We first define a function that calculates the final proportions of ancestry
# based on the migration matrix


def propfrommig(mig):
    curr = mig[-1, :]
    for row in mig[-2::-1, :]:
        curr = curr*(1-numpy.sum(row))+row
    return curr


# Fix parameter tstart to 16 GA

def ppx_xxp_pxx_fix_tstart(params, fracs):
    """African and last European proportions fixed by ancestry proportions 
    a simple model in which populations Eu and NAT arrive discretely at first generation, Af at a subsequent generation, 
    and Eu discretely yet later. 
    . If a time is not integer, the migration 
    is divided between neighboring times proportional to the non-integer time fraction. 
    We'll assume population 3 
    replaces after the replacement from population 1 and 2 if they arrive at same generation. 
    The first generation is equally composed of pops 1 and 2
    order is ['CEU','NAT','YRI']"""

    (init_Eu, afam_time, nuEu_time) = params
    tstart = 0.15  # Here, I fix this parameters to 16 GA

    # print "times ",tstart,afam_time,nuEu_time

    def fun(props):
        (afam_prop, nuEu_prop) = props
        return propfrommig(ppx_xxp_pxx((init_Eu, tstart, afam_prop, afam_time, nuEu_prop, nuEu_time)))[0:2]-fracs[0:2]

    (afam_prop, nuEu_prop) = scipy.optimize.fsolve(fun, (.2, .2))
    return ppx_xxp_pxx((init_Eu, tstart, afam_prop, afam_time, nuEu_prop, nuEu_time))


def outofbounds_ppx_xxp_pxx_fix_start(params, fracs):
    # constraint function evaluating below zero when constraints not satisfied
    # print "in constraint  outofbounds_211_cont_unif_params"
    ret = 1
    # (init_Eu,tstart,afam_time,nuEu_time)=params[0]
    (init_Eu, afam_time, nuEu_time) = params
    tstart = 0.15

    def fun(props):
        (afam_prop, nuEu_prop) = props
        return propfrommig(ppx_xxp_pxx((init_Eu, tstart, afam_prop, afam_time, nuEu_prop, nuEu_time)))[0:2]-fracs[0:2]
    (afam_prop, nuEu_prop) = scipy.optimize.fsolve(fun, (.2, .2))
    return outofbounds_ppx_xxp_pxx((init_Eu, tstart, afam_prop, afam_time, nuEu_prop, nuEu_time))
