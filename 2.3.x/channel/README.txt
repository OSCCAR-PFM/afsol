Tue Sep 2 09:18:31 CEST 2014

GAW: Changed mesh thickness from 0.35 m to 0.035 m, so that mesh is 2D slab and
     cell aspect ratio close to 1.

GAW: Changed frontBack patch name to more conventional frontAndBack.

GAW: Cleaned up formatting of constant/polyMesh/blockMeshDict.

GAW: Cleaned up formatting of dictionary files in constant and system.

GAW: Adjusted probes function object to new (thinner) domain thickness.

GAW: Adjusted time step size to be fixed at 5e-5 s. This runs with a Courant
     number of around 0.57.
