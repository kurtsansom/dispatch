Compiling
=========

To comile, go to the expriments/example directory, and use ``make info``
to see macro names and values that would be used in ``make``:::

   cd experiments/turbulence
   make info

If the automatically selected make options shown are OK, just do:::

   make -j

If the compiler chosen by default is not appropriate, do instead, 
for example:::

   make COMPILER=gfortran -j

.. toctree::
   :maxdepth: 3

   option_groups
   options

