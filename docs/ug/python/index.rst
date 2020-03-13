.. _python:

Python support
==============

To import Python modules, add the full path of the DISPATCH directory ``utilities/python/``
to the Python path (see ref:`environment`), and do::

   import dispatch
   import dispatch.select
   import dispatch.graphics
   import dispatch.yt

The modules are self-documented, using the standard Python documentation
mechanisms, cf.
::

   help(dispatch)
   help(dispatch.select)
   dispatch.snapshot(<TAB>

The source code may be found in ``utilities/python/dispatch/``, and is
also available via the :ref:`auto-generated documentation <doxy>` (the
"Files" menu contains source code).

.. toctree::
   :maxdepth: 4

   dispatch
   dispatch_select
   dispatch_graphics
   dispatch_yt
   notebooks
   aux_data
