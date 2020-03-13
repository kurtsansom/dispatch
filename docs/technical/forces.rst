
Forces
------

Forces are now part of the extras_mod, so to add forcing take a copy
of ``extras/extras_mod.f90`` into the ``experiment/xxx/`` dir, and uncomment
the lines marked `... ! forces`

The new ``forces_mod`` does NOT collide with the old ``force_mod`` modules, so
for compatibility one can construct new forces_mod modules by using ``force%selected``
to load values into ``patch%force_per_unit_mass``.

However, in the future ``forces_mod.f90`` files should be written directly,
since this alleviates the need to stick to a definitive (and long) argument 
list for the ``force%selected`` functions.

.. toctree::
   :maxdepth: 3

 
