# clauses (All the triples (3j,3j+1,3j+2), selecting 3 out of n spins, form a clause)
n=8
#clauses=[1,2,3,2,4,8,2,5,7,6,7,8,4,6,7]
clauses=[0,1,2,1,3,7,1,4,6,5,6,7,3,5,6]
m=int(len(clauses)/3)

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

# Run on quantum solver
from dwave.system import EmbeddingComposite, DWaveSampler
# Define the sampler that will be used to run the problem
# qpu_advantage = DWaveSampler(solver={'topology__type': 'pegasus'})
qpu_2000q = DWaveSampler(solver={'topology__type': 'chimera'})
sampler = EmbeddingComposite(qpu_2000q)

# Run the problem on the sampler and print the results
sampleset = sampler.sample_ising(h,J,
                           num_reads = 1,
                           annealing_time = 1,
                           label='EC3')
print(sampleset)

guess=sampleset.first.sample
print(guess)
