Files needed to run Tinker with aenet potentials
================================================

To use Tinker with aenet potentials, a parameter file with type definitions is needed.  Such a file is provided here named :code:`aenet.prm`.

Use this parameter file in your Tinker `key` file as follows:

.. code::

   # use AENET to compute energies & forces
   EXTRATERM ONLY
   PARAMETERS aenet.prm
   
   # Parallel run using 16 threads
   OPENMP-THREADS 16


See also the example provided in the parent directory.
