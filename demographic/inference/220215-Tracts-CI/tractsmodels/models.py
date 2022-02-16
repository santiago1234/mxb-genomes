import numpy
import scipy

"""
A simple model in which populations 1 and 2 arrive discretely at first
generation, 3 at a subsequent generation. If a time is not integer, the
migration is divided between neighboring times proportional to the
non-integer time fraction.  We'll assume population 3 still replaces
migrants from 1 and 2 after the replacement from population 1 and 2 if they
arrive at same generation.
"""


def ppx_xxp(*params):
    """ A simple model in which populations 1 and 2 arrive discretely at first
        generation, 3 at a subsequent generation. If a time is not integer,
        the migration is divided between neighboring times proportional to the
        non-integer time fraction.  We'll assume population 3 still replaces
        migrants from 1 and 2 after the replacement from population 1 and 2 if
        they arrive at same generation.

        Parameters: (prop1, tstart, prop3, t3)
        In this prop1 is the initial proportion from population 1,
        prop2=1-prop1, tstart is the arrival times of pops (1,2) t3 is the
        arrival time of pop 3

        The two times are measured in units of 100 generations, because some
        python optimizers work better when all parameters have the same scale.
        """
    (prop1, tstart, prop3, t3) = params[0]

    tstart *= 100
    t3 *= 100

    # some sanity checks
    if t3 > tstart or t3 < 0 or tstart < 0:
        # This will be caught by "outofbounds" function. Return empty matrix
        gen = int(numpy.ceil(max(tstart, 0))) + 1
        mig = numpy.zeros((gen + 1, 3))
        return mig

    # How many generations we'll need to accomodate all the migrations. The +1
    # is not really necessary here.
    gen = int(numpy.ceil(tstart)) + 1

    # How far off the continuous time is from its discrete optimization
    timefrac = gen - tstart - 1

    # Build an empty matrix with enough generations to handle all migrations
    mig = numpy.zeros((gen + 1, 3))

    # replace a fraction at first and second generation to ensure a continuous
    # model
    prop2 = 1 - prop1
    mig[-1, :] = numpy.array([prop1, prop2, 0])

    interEu = prop1 * timefrac
    interNat = prop2 * timefrac
    mig[-2, :] = numpy.array([interEu, interNat, 0])

    # Which integer generation to add the migrants from pop 3
    gen3 = int(numpy.ceil(t3)) + 1
    timefrac3 = gen3 - t3 - 1

    # we want the total proportion replaced  by 3 to be prop3. We therefore add
    # a fraction f at generation gen-1, and (prop3-f)/(1-f) at generation gen.

    mig[gen3 - 1, 2] = timefrac3 * prop3
    mig[gen3, 2] = (prop3 - timefrac3 * prop3) / (1 - timefrac3 * prop3)

    return mig


def outofbounds_ppx_xxp(*params):
    """ Constraint function evaluating below zero when constraints are not
        satisfied. """
    ret = 1

    (prop1, tstart, prop3, t3) = params[0]

    ret = min(1-prop1, 1-prop3)
    ret = min(ret, prop1, prop3)

    # Pedestrian way of testing for all possible issues
    func = ppx_xxp
    mig = func(params[0])
    totmig = mig.sum(axis=1)
    # print  "ret1=",ret

    if prop1 > 1 or prop3 > 1:
        print("Pulse greater than 1")
    if prop1 < 0 or prop3 < 0:
        print("Pulse less than 0")
    # print  "ret2 ",ret

    ret = min(ret, -abs(totmig[-1]-1) + 1e-8)
    ret = min(ret, -totmig[0], -totmig[1])
    # print "ret3 " , ret

    ret = min(ret, min(1-totmig), min(totmig))
    # print "ret4 " , ret

    # print "times ",t3,tstart
    ret = min(ret, tstart-t3)

    ret = min(ret, t3)
    ret = min(ret, tstart)
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


def ppx_xxp_fix(times, fracs):
    (tstart, t3) = times
    # An auxiliary function that, given the "missing" parameters, returns the
    # full migration matrix

    def fun(props):
        (prop3, prop1) = props
        return propfrommig(ppx_xxp((prop1, tstart, prop3, t3)))[0:2] - fracs[0:2]

    # Find the best-fitting parameters
    # (.2,.2) is just the starting point for the optimization function, it
    # should not be sensitive to this, but it's better to start with reasonable
    # parameter values.
    (prop3, prop1) = scipy.optimize.fsolve(fun, (.2, .2))
    # make sure funding generation sum to 1
    mig = ppx_xxp((prop1, tstart, prop3, t3))
    mig[-1, :] = mig[-1, :] / mig[-1, :].sum()
    return mig


def outofbounds_ppx_xxp_fix(params, fracs):
    # constraint function evaluating below zero when constraints not satisfied
    # print "in constraint  outofbounds_211_cont_unif_params"

    ret = 1

    (tstart, t3) = params
    if tstart > 1:
        print("time above 500 generations!")
        return (1 - tstart)

    def fun(props):
        (prop3, prop1) = props
        return propfrommig(ppx_xxp((prop1, tstart, prop3, t3)))[0:2] - fracs[0:2]
    # (.2,.2) is just the starting point for the optimization function, it should not be sensitive to this, but it's better to start with reasonable parameter values.
    (prop3, prop1) = scipy.optimize.fsolve(fun, (.2, .2))
    return outofbounds_ppx_xxp((prop1, tstart, prop3, t3))


startparams_ppx_xxp = [
    0.0667324402332306,  # t1
    0.055764797329902666  # t2
]
# => =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>


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


def ppx_xxp_pxx_fix(params, fracs):
    """African and last European proportions fixed by ancestry proportions 
    a simple model in which populations Eu and NAT arrive discretely at first generation, Af at a subsequent generation, 
    and Eu discretely yet later. 
    . If a time is not integer, the migration 
    is divided between neighboring times proportional to the non-integer time fraction. 
    We'll assume population 3 
    replaces after the replacement from population 1 and 2 if they arrive at same generation. 
    The first generation is equally composed of pops 1 and 2
    order is ['CEU','NAT','YRI']"""

    (init_Eu, tstart, afam_time, nuEu_time) = params

    # print "times ",tstart,afam_time,nuEu_time

    def fun(props):
        (afam_prop, nuEu_prop) = props
        return propfrommig(ppx_xxp_pxx((init_Eu, tstart, afam_prop, afam_time, nuEu_prop, nuEu_time)))[0:2]-fracs[0:2]

    (afam_prop, nuEu_prop) = scipy.optimize.fsolve(fun, (.2, .2))
    mig = ppx_xxp_pxx((init_Eu, tstart, afam_prop,
                      afam_time, nuEu_prop, nuEu_time))
    mig[-1, :] = mig[-1, :] / mig[-1, :].sum()
    return mig


def outofbounds_ppx_xxp_pxx_fix(params, fracs):
    # constraint function evaluating below zero when constraints not satisfied
    # print "in constraint  outofbounds_211_cont_unif_params"
    ret = 1
    # (init_Eu,tstart,afam_time,nuEu_time)=params[0]
    (init_Eu, tstart, afam_time, nuEu_time) = params

    def fun(props):
        (afam_prop, nuEu_prop) = props
        return propfrommig(ppx_xxp_pxx((init_Eu, tstart, afam_prop, afam_time, nuEu_prop, nuEu_time)))[0:2]-fracs[0:2]
    (afam_prop, nuEu_prop) = scipy.optimize.fsolve(fun, (.2, .2))
    return outofbounds_ppx_xxp_pxx((init_Eu, tstart, afam_prop, afam_time, nuEu_prop, nuEu_time))


startparams_ppx_xxp_pxx = [
    0.16347885,
    0.08063866,
    0.0799999,
    0.06329345
]
# => =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>


def ccx_xxp(*params):
    # a simple model in which populations 1 and 2 arrive continuously
    # starting at time t1 and population 3 discretely starting at
    # time t2, at rates a and b, respectively. Note that if t2>t1 is
    # equivalent to t2=t1 (since there is initially a single population).
    # It's therefore preferable to impose t1<t2. If a time is not integer,
    # the migration is divided between neighboring times proportional to the ancestry.
    # we'll assume population 3 replaces after the replacement from population 1 and 2 if
    # they arrive at same generation. The first generation is equally composed of pops 1 and 2
    (frac1, frac2, t1, frac3, t2) = params[0]
    t1 *= 100
    t2 *= 100

    if t2 > t1:
        print("t2>t1:", t2, ">", t1, "should be caught by outofbounds_func")
        gen = max(int(numpy.ceil(t1)), 2)+1
        mig = numpy.zeros((gen+1, 3))
        return mig

        #print (a,b)

    gen = int(numpy.ceil(t1))+1
    frac = gen-t1-1
    mig = numpy.zeros((gen, 3))
    # replace a fraction at second generation to ensure a continuous model distribution with generation time. The first generation is always in complete replacement. The second generation must mostly replace to ensure continuity.
    # frac=0 is exact at gen, so no correction needed at generation gen-1.
    # gen-1 interpolates between the continuous fractions and the full replacement.

    mig[-1, :] = numpy.array([frac1, frac2, 0])

    # contents of the second generation
    inter1 = frac*mig[-1, 0]+(1-frac)*frac1
    inter2 = frac*mig[-1, 1]+(1-frac)*frac2

    mig[-2, :] = numpy.array([inter1, inter2, 0])
    mig[2:-2, 0] = frac1
    mig[2:-2, 1] = frac2

    # print a,inter,mig

    # now proceed with the pulse population
    gen = int(numpy.ceil(t2))
    frac = gen-t2
    # replacement from pop 3 occurs after replacement from pop 1,2, and may span 2 generations
    # mig[gen]*=(1-frac3*(1-frac))
    mig[gen, 2] = frac3*(1-frac)
    mig[gen-1, 2] = frac3*frac

    mig[gen, 2] = (frac3*(1-frac))/(1-frac*frac3)
    mig[gen-1, 2] = frac*frac3
    # normalize arrival generation
    mig[-1, :] /= mig[-1, :].sum()

    return mig

    # print mig
    # return mig


def outofbounds_ccx_xxp(*params):
    # constraint function evaluating below zero when constraints not satisfied
    # print "in constraint  outofbounds_211_cont_unif_params"

    (frac1, frac2, t1, frac3, t2) = params[0]
    ret = 1-(frac3+frac2+frac1)
    ret = min(ret, frac1, frac2, frac3)

    # print "ret0",  ret

    # pedestrian way of testing for all possible issues
    func = ccx_xxp
    mig = func(params[0])
    totmig = mig.sum(axis=1)
    # print  "ret1=",ret

    # print  "ret2 ",ret
    ret = min(ret, -abs(totmig[-1]-1)+1e-8)
    ret = min(ret, -totmig[0], -totmig[1])
    # print "ret3 " , ret
    ret = min(ret, min(1-totmig), min(totmig))
    # print "ret4 " , ret

    # print "times ",afam_time,tstart
    ret = min(ret, t1-t2)
    ret = min(ret, t1)
    ret = min(ret, t2)

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


startparams_ccx_xxp = [
    0.0533,
    0.0533,
    0.1536,
    0.0913,
    0.11673
]
# => =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>


def ppx_ccx_xxp(*params):
    # a simple model in which populations 1 and 2 arrive continuously
    # starting at time t1 and population 3 discretely starting at
    # time t2, at rates a and b, respectively. Note that if t2>t1 is
    # equivalent to t2=t1 (since there is initially a single population).
    # It's therefore preferable to impose t1<t2. If a time is not integer,
    # the migration is divided between neighboring times proportional to the ancestry.
    # we'll assume population 3 replaces after the replacement from population 1 and 2 if
    # they arrive at same generation. The first generation is equally composed of pops 1 and 2

    # The initEu_frac represents the propotion of European migrants in the very first generation.

    (initEu_frac, frac1, frac2, t1, frac3, t2) = params[0]
    t1 *= 100
    t2 *= 100

    if t2 > t1:
        print("t2>t1:", t2, ">", t1, "should be caught by outofbounds_func")
        gen = max(int(numpy.ceil(t1)), 2)+1
        mig = numpy.zeros((gen+1, 3))
        return mig

        #print (a,b)

    gen = int(numpy.ceil(t1))+1
    frac = gen-t1-1
    mig = numpy.zeros((gen, 3))
    # replace a fraction at second generation to ensure a continuous model distribution with generation time. The first generation is always in complete replacement. The second generation must mostly replace to ensure continuity.
    # frac=0 is exact at gen, so no correction needed at generation gen-1.
    # gen-1 interpolates between the continuous fractions and the full replacement.

    mig[-1, :] = numpy.array([frac1, frac2, 0])

    # contents of the second generation
    inter1 = frac*mig[-1, 0]+(1-frac)*frac1
    inter2 = frac*mig[-1, 1]+(1-frac)*frac2

    mig[-2, :] = numpy.array([inter1, inter2, 0])
    mig[2:-2, 0] = frac1
    mig[2:-2, 1] = frac2

    # print a,inter,mig

    # now proceed with the pulse population
    gen = int(numpy.ceil(t2))
    frac = gen-t2
    # replacement from pop 3 occurs after replacement from pop 1,2, and may span 2 generations
    # mig[gen]*=(1-frac3*(1-frac))
    mig[gen, 2] = frac3*(1-frac)
    mig[gen-1, 2] = frac3*frac

    mig[gen, 2] = (frac3*(1-frac))/(1-frac*frac3)
    mig[gen-1, 2] = frac*frac3

    # Set the proportions for the first arriving generation
    mig[-1, :] = numpy.array([initEu_frac, 1 - initEu_frac, 0])

    return mig

    # print mig
    # return mig


def outofbounds_ppx_ccx_xxp(*params):
    # constraint function evaluating below zero when constraints not satisfied
    # print "in constraint  outofbounds_211_cont_unif_params"

    (initEu_frac, frac1, frac2, t1, frac3, t2) = params[0]
    ret = 1-(frac3+frac2+frac1)
    ret = min(ret, frac1, frac2, frac3)

    # print "ret0",  ret

    # pedestrian way of testing for all possible issues
    func = ppx_ccx_xxp
    mig = func(params[0])
    totmig = mig.sum(axis=1)
    # print  "ret1=",ret

    # print  "ret2 ",ret
    ret = min(ret, -abs(totmig[-1]-1)+1e-8)
    ret = min(ret, -totmig[0], -totmig[1])
    # print "ret3 " , ret
    ret = min(ret, min(1-totmig), min(totmig))
    # print "ret4 " , ret

    # print "times ",afam_time,tstart
    ret = min(ret, t1-t2)
    ret = min(ret, t1)
    ret = min(ret, t2)

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


startparams_ppx_ccx_xxp = [
    0.2140,
    0.0533,
    0.0533,
    0.1536,
    0.0913,
    0.11673
]
# => =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>
# => =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>=> =>


MODELS = {
    'ppx_xxp': (ppx_xxp_fix, outofbounds_ppx_xxp_fix, startparams_ppx_xxp),
    'ppx_xxp_pxx': (ppx_xxp_pxx_fix, outofbounds_ppx_xxp_pxx_fix, startparams_ppx_xxp_pxx),
    'ccx_xxp': (ccx_xxp, outofbounds_ccx_xxp, startparams_ccx_xxp),
    'ppx_ccx_xxp': (ppx_ccx_xxp, outofbounds_ppx_ccx_xxp, startparams_ppx_ccx_xxp)
}
