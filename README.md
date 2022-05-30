# dwave_bias
Codes to explore quantum annealing on the D-wave in the presence of bias fields

Annealing with bias has been proposed in: Phys. Rev. Lett. 123, 120501 (2019).
To test the algorithm, instances of exact cover 3 are used.

- ec3.py runs D-wave sampler for 1 EC3 problem defined in the code in standard way
- ec3_itbias.py runs D-wave sampler for 1 EC3 problem using as bias outcome from earlier anneal 
- ec3many.py runs D-wave sampler for several EC3 instances imported from file with instance definitions.

The instance definition files are called datN?M?.cvs for instances with N qubits and M clauses.
All these instances have unique satisfying assignments.
The cvs files contain 3M numbers per line; 
every (of the 100) lines representing one problem instance;
numbers {3i-2,3i-1,3i} of each line select the spins of a clause.
