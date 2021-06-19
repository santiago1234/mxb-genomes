import numpy


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
