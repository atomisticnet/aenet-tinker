=====================
ænet-Tinker Interface
=====================

If you make use of the aenet-Tinker interface, please cite the following reference:

M.S. Chen, T. Morawietz, T.E. Markland, N. Artrith, AENET-LAMMPS and AENET-TINKER: Interfaces for Accurate and Efficient Molecular Dynamics Simulations with Machine Learning Potentials, J. Chem. Phys. 155, 074801 (2021). doi: https://doi.org/10.1063/5.0063880

Preprint (Preprint where the data is discussed)
arXiv (2021) doi: https://arxiv.org/abs/2107.11311

The database can be obtained from the Materials Cloud repository, DOI: https://doi.org/10.24435/materialscloud:dx-ct

ænet-Tinker
-------------

This directory provides all files needed to link the Tinker_ molecular
dynamics code with ``aenetLib``.

The present interface code has been tested with Tinker version 8.2.1
(from June 2017).

Note
-------------

The interface supports OpenMP_ shared memory parallelism, but only using
the `GNU Fortran`_ compiler.  The OpenMP implementation provided by
`Intel compilers`_ (as of version 15) is not compatible.  However, even
using the GNU compiler ``aenetLib`` a significant speed-up can be
achieved by linking with the `Intel MKL`_ library when available.

Presently, only energies and forces are provided by the interface.  This
means Tinker functionality that requires derivatives of higher orders
are not supported.

Prerequisites
-------------

To be linked with Tinker, ænet has to be compiled as library.  Build the
library by running

.. code-block:: sh

  make -f ./makefiles/Makefile.gfortran_serial lib

in the ænet ``src`` directory.  Alternatively, if the `Intel MKL`_
library is available, ``Makefile.gfortran_mkl_serial`` can be used.

Building Tinker with ænet support
---------------------------------

Copy the following files to the Tinker ``source`` directory:

- aenettinker.f90
- extra.f
- extra1.f
- Makefile.aenetlib_gfortran

and adjust the path to aenet (``AENETDIR``) in the Makefile.  The
Makefile contains two alternative blocks for numerical libraries.  To
use the system ``LAPACK`` and ``BLAS`` libraries uncomment the line

.. code-block:: makefile

  NUMLIB = -llapack -lblas

for use with the `Intel MKL`_ library, use the following instead

.. code-block:: makefile

  NUMLIB = -Wl,--start-group \
           $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
           $(MKLROOT)/lib/intel64/libmkl_core.a \
           $(MKLROOT)/lib/intel64/libmkl_sequential.a \
           -Wl,--end-group -lpthread -lm -ldl

Tinker can then be compiled with

.. code-block:: sh

  make -f Makefile.aenetlib_gfortran

from within the ``source`` directory.

Running Tinker simulations using ænet potentials
------------------------------------------------

To run simulations using ænet artificial neural network (ANN)
potentials, the parameter file ``aenet.prm`` can be used.  It only
contains the masses of all chemical species, as no further information
is required.

All ANN potential files have to follow the naming convention

::

   <species>.ann

where ``<species>`` is the chemical symbol (``H``, ``He``, ``Li``, etc.)
and have to be present in the working directory.

In the Tinker ``key`` principal input file, ænet potentials are
activated with the keyword

::

  EXTRATERM only

The number of threads used in parallel runs can be controlled with the
keyword ``OPENMP-THREADS``.

See also the example Tinker input file ``tinker.key``.

.. _Tinker: https://dasher.wustl.edu/tinker/
.. _OpenMP: http://www.openmp.org
.. _GNU Fortran: https://gcc.gnu.org/fortran/
.. _Intel compilers: https://software.intel.com/en-us/intel-compilers
.. _Intel MKL: https://software.intel.com/en-us/mkl
