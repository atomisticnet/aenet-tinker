### Input Files for Molecular Dynamics (MD) Simulation using the aenet-Tinker Interface

* `Li.ann`, `Si.ann`

   ANN potential files. The original files from the ANN potential training `Li.15t-15t.nn-04153` and `Si.15t-15t.nn-04153` were renamed.

* `aenet.prm`

   File with atom type parameters for each atomic species.

* `md.key`

   Tinker key file with parameters of the MD simulation. Some notes:
   
   `OPENMP-THREADS 64` requests an MD simulation using 64 parallel processes.
   
   The following two lines activate the use of aenet's ANN potentials:
   
   `EXTRATERM ONLY`<br/>
   `PARAMETERS aenet.prm`
   
   The parameter file can also be directly specified when `dynamic.x` is run on the command line.  
   This technique is used in the script `mdrun-NVT.sh`.
   
   Tinker's optimizer does, per default, not expect large negative energies.  
   The following statement sets the minimal expected energy to a very small value:
   
   `FCTMIN -50000000`
   
* `md.xyz`

   Initial amorphous LiSi atomic structure with 1024 atoms in Tinker's XYZ format.
   
* `mdrun-NVT.sh`

   Bash script that runs the actual MD simulation in the NVT ensemble.

