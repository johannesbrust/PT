# PT
**Algorithms for COVID-19 Pooling Tests (PT)**

MATLAB implementations of algorithms from the article

"COVID-19 Pooling Matrix Designs", J.J. Brust (2022), IMAGE (68)

Content:
  * ALGS/
    - FGD_PT.m (Finite Geometry Design (FGD): {Sparse algorithm to construct finite affine and projective planes (prime powers)})
    - SDD_PT.m (Shifted Diagonal Design (SDD): {Algorithm to construct pooling matrix M (prime number pools)})    
    - SCOMP_PT.m (Sequential Combinatorial Orthogonal Matching Pursuit (SCOMP): {Decoding Algorithm})
    - DD_PT.m (Definite Defectives (DD): {Decoding Algorithm}
    - COMP_PT.m (Combinatorial Orthogonal Matching Pursuit (COMP): {Decoding Algorithm})
  * TESTS/
    - example_Pooling.m (Example pooling test with m=4, k=4 using FGD and COMP)
    - example_Large.m (Large test on identifying 100 positives from 100000 samples (matrix about 3 seconds) )
    - test_COMP_PT_1.m (Test of COMP on an example matrix from Wikipedia)
  * AUXILIARY/
  * EXTERNAL/

## Example
You can run an example test from within the /TESTS folder. On
the MATLAB command prompt:

`> example_Pooling`

```
----------- FGD Algorithm ----------- 
Problem Size 
Pools (m):                     4 
Disjunct (k):                  4 
Expected tests:                20 
OUTPUT(S):############## 
Time (constr):        3.536667e-02 


----------- COMP Algorithm ----------- 
Problem Size 
Tests:                 20 
Samples:               16 
Input Expct. Positive: 4 
OUTPUTS:################ 
Time (search):         9.010810e-04 
Positive items:        4 
######################## 
Identified Indices:    2 
Identified Indices:    3 
Identified Indices:    7 
Identified Indices:    13 
