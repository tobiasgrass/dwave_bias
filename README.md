# dwave_bias
Codes to explore quantum annealing on the D-wave in the presence of bias fields

Annealing with bias has been proposed in: Phys. Rev. Lett. 123, 120501 (2019).
To test the algorithm, instances of exact cover 3 are used.

- ec3.py runs D-wave sampler for 1 EC3 problem defined in the code in standard way
- ec3_itbias.py runs D-wave sampler for 1 EC3 problem using as bias outcome from earlier anneal 

