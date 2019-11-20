Compiling
=========

Typically, you need to use a ``module`` command to get access to a compiler
and MPI library.  This could be, for example::

   module load intel openmpi/intel

To comile, go to one of the expriments/ directories, and use ``make info``
to see macro names and values that would be used in ``make``:::

   cd experiments/turbulence
   make info

If the automatically selected make options shown are OK, just do:::

   make -j

If the compiler chosen by default is not appropriate, do instead, 
for example:::

   make COMPILER=iotran -j

.. toctree::
   :maxdepth: 3

   option_groups
   options

