import numpy as np
from numpy import genfromtxt
data = genfromtxt('datN26.csv', delimiter='\t') # delimiter=' '

results=[0,0,1,2,3]
for k in range(22,29): # range(0,100):
    # clauses (All the triples (3j,3j+1,3j+2), selecting 3 out of n spins, form a clause)
    n=26
    clauses=data[k,:]
    m=int(len(clauses)/3)
    clauses=clauses-1

    # couplings
    h = {}
    J = {}
    for i in range(0,n):
        h[i]=0
        for j in range(0,n):
            J[i,j]=0
    for i in range(0,m):
        J[clauses[3*i],clauses[3*i+1]]=J[clauses[3*i],clauses[3*i+1]]+2
        J[clauses[3*i],clauses[3*i+2]]=J[clauses[3*i],clauses[3*i+2]]+2
        J[clauses[3*i+1],clauses[3*i+2]]=J[clauses[3*i+1],clauses[3*i+2]]+2
    for i in range(0,len(clauses)):
        h[clauses[i]]=h[clauses[i]]-2

    # Run on classical solver
    from dimod.reference.samplers import ExactSolver
    sampler = ExactSolver()
    sampleset = sampler.sample_ising(h,J)
    print(sampleset.lowest(atol=.001))
    oguess=sampleset.first.sample  # get best guess
    s1guess=dict(oguess)
    s1guess[2]=-1*oguess[2]  # sub-optimal guess 1: flip one spin
    s2guess=dict(s1guess)
    s2guess[4]=-1*oguess[4]  # sub-optimal guess 2: flip two spins
    s3guess=dict(s2guess)
    s3guess[5]=-1*oguess[5]  # sub-optimal guess 3: flip three spins
    targetE=sampleset.first.energy

    print(oguess)
    print(s1guess)
    print(s2guess)
    print(s3guess)
    print(targetE)

    # Run on quantum solver
    from dwave.system import EmbeddingComposite, DWaveSampler
    # Define the sampler that will be used to run the problem
    # qpu_advantage = DWaveSampler(solver={'topology__type': 'pegasus'})
    qpu_2000q = DWaveSampler(solver={'topology__type': 'chimera'})
    sampler = EmbeddingComposite(qpu_2000q)
    nr=30 # num_reads
    at=1 # annealing time

    # Sampling without bias
    sampleset = sampler.sample_ising(h,J,
                            num_reads = nr,
                            annealing_time = at,
                            label='EC3_nobias')
    print(sampleset)

    #Biased sample
    vbias=1
    #h_schedule = [[0.0, 1], [1, 0]]
    hbias={}

    for i in range(0,n):
        hbias[i]=h[i]-vbias*oguess[i]
    sampleset0 = sampler.sample_ising(hbias,J,
                                num_reads = nr,
                                annealing_time = at,
                                label='EC3_oguess')
    print(sampleset0)

    for i in range(0,n):
        hbias[i]=h[i]-vbias*s1guess[i]
    sampleset1 = sampler.sample_ising(hbias,J,
                                num_reads = nr,
                                annealing_time = at,
                                label='EC3_s1guess')
    print(sampleset1)

    for i in range(0,n):
        hbias[i]=h[i]-vbias*s2guess[i]
    sampleset2 = sampler.sample_ising(hbias,J,
                                num_reads = nr,
                                annealing_time = at,
                                label='EC3_s2guess')
    print(sampleset2)

    for i in range(0,n):
        hbias[i]=h[i]-vbias*s3guess[i]
    sampleset3 = sampler.sample_ising(hbias,J,
                                num_reads = nr,
                                annealing_time = at,
                                label='EC3_s3guess')
    print(sampleset3)
    print(oguess)

    # average energy of found solution:
    def figomer(sampleset,num_reads,targetE):
        etot=0
        array=sampleset.record
        print(array)
        sucs=0
        for k in range(0,len(array)):
            en=0
            for i in range(0,n):
                en=en+h[i]*array[k][0][i]
                for j in range(0,n):
                    en=en+J[i,j]*array[k][0][i]*array[k][0][j]
            if abs(en-targetE)<0.001:
                sucs=sucs+array[k][2]
            etot=etot+(en-targetE)*array[k][2]
        output=[etot/num_reads,sucs/num_reads]    
        return output
    A=figomer(sampleset,nr,targetE)
    A0=figomer(sampleset0,nr,targetE)
    A1=figomer(sampleset1,nr,targetE)
    A2=figomer(sampleset2,nr,targetE)
    A3=figomer(sampleset3,nr,targetE)

    eav=np.array([A[0],A0[0],A1[0],A2[0],A3[0]])
    sav=np.array([A[1],A0[1],A1[1],A2[1],A3[1]])
    print(eav)
    print(sav)

    # Save results:
    f = open("eresult_N26m26_100.txt", "a")
    np.savetxt(f,[eav])
    f.close()
    f = open("sresult_N26m26_100.txt", "a")
    np.savetxt(f,[sav])
    f.close()


