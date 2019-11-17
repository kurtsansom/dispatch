Python support
==============

To load Python modulesi, add the full path of the DISPATCH directory ``utilities/python/``
to the Python path (environment variable ``PYTHONPATH``), and do::

    import dispatch
    import dispatch.yt
    import dispatch.select
    import dispatch.graphics

The modules are self-documented, cf:

.. toctree::
   :maxdepth: 4

   dispatch
   dispatch_yt
   dispatch_select
   dispatch_graphics
   notebooks
   aux_data

The source code may be found in ``utilities/python/dispatch/``, and is
also available via the :ref:`auto-generated documentation <doxy>` (the
"Files" menu contains source code).
