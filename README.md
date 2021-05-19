# Lipschitz-Hankel Integrals
 
A set of Lipschitz-Hankel integrals in MATLAB.

From a fundamental integral (J_000),
other integrals were obtained through differentiation and recursion relationships 
(see [Eason, Noble & Sneddon](https://royalsocietypublishing.org/doi/abs/10.1098/rsta.1955.0005)),
via a Mathematica script.

The script was then converted to MATLAB. Tests were conducted to make sure
that the results were the same for both of them.

Key integrals with simple analytical solutions were also tested.

```matlab
Jmnp = LH_INTEGRALS(r,z,alpha,K,E,cas);
``` 
calculates the Jmnp LH integral at the point (r,z) and at a cut angle 'alpha' (see [Paynter et al., 2009](https://www.sciencedirect.com/science/article/pii/S0020768308003648)).

The integral is characterised by the case 'cas', passed as a string of the integral index 'mnp'. For example, J<sub>0,0,1</sub> at a point (0,1) is:

```matlab
J001 = LH_INTEGRALS(0,1,0,K,E,'001');
```

'K' and 'E' are the complete elliptic integrals of the first and second kind for the given point (r,z). It is calculated outside of the fucntion to improve computational performance, since multiple integrals can use the same 'K' and 'E' as long as they are being calculated for the same (r,z) point. It also provides opportunity for optimisation of 'K' and 'E' calculations outside of the LH_INTEGRALS script.

Some integrals require also the calculation of the complete elliptic integral of the third kind. Since not all of them require this, it is calculated when needed.

The MATLAB script can handle the following integrals:


* J<sub>0,0,0</sub>, J<sub>0,0,1</sub>, J<sub>0,0,2</sub>, J<sub>0,0,3</sub>
* J<sub>0,1,-1</sub>, J<sub>0,1,0</sub>, J<sub>0,1,1</sub>, J<sub>0,1,2</sub>, J<sub>0,1,3</sub>
* J<sub>0,2,0</sub>, J<sub>0,2,1</sub>
* J<sub>1,0,-1</sub>, J<sub>1,0,0</sub>, J<sub>1,0,1</sub>, J<sub>1,0,2</sub>, J<sub>1,0,3</sub>
* J<sub>1,1,-1</sub>, J<sub>1,1,0</sub>, J<sub>1,1,1</sub>, J<sub>1,1,2</sub>, J<sub>1,1,3</sub>
* J<sub>1,2,-1</sub>, J<sub>1,2,0</sub>, J<sub>1,2,1</sub>, J<sub>1,2,2</sub>
* J<sub>2,0,-1</sub>, J<sub>2,0,0</sub>, J<sub>2,0,1</sub>, J<sub>2,0,2</sub>, J<sub>2,0,3</sub>
* J<sub>2,1,0</sub>, J<sub>2,1,1</sub>, J<sub>2,1,2</sub>, J<sub>2,1,3</sub>
* J<sub>2,2,0</sub>, J<sub>2,2,1</sub>

In cases where there is a need of calculating the ratio between the integral and the r coordinate (i.e. Jmnp/r), special care must be taken.
We must take the limits properly as r approaches zero. For this reason, for the integrals that might be singular when r goes to zero, an alternative form is available as 'Jmnpbyr'.
For example, the integral J<sub>010</sub> cannot be divided by r without proper care. In this case, J010/r at a point (r,z) is given by:

```matlab
LH_INTEGRALS(r,z,alpha,K,E,'010byr');
```

