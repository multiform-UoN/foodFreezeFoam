Axialsymmetric heat exchanger with meat
---------------------------------------

- The 0/ field files contain nonsense patchFields.
  All interesting work is done using the changeDictionaryDicts.

- Base Transient case with:
	- air at T= 225 K (-48 degC) with constant properties. Set up with k-epsilon model and 0.1 m/s velocity jet. 
	- the case is initialised using hydraulic diameter and turbulence intensity of 5%  (standard)
	- full setup of post-processing, apart from local Nusselt, Reynolds and Peclet numbers 
        - fully resolved meat advecting through the domain at approx 0.01 m/s and T = 303 K (30 degC).
	- Note: 
	  - Domain discretised to the second order schemes. 
	  - Mesh indendence was not run.
	  - Overall - standard k-omegaSST results in this case should be taken with caution due to the mesh resolution.



