Derefine event sequence
------------------------

* In ``remove_patch()``

  * remove the task from the task list (it will remain in garbage bin while needed)
  * use the existing nbor list to set bits%init_nbors on nbors
  * now that the task is removed, check if nbors have become ready
  * move the task to the garbage bin

* In ``task_mesg_t%npack()``

  * detect bits%init_nbors, inside detection of suicide note
  * set bits%init_nbors on the nbors (under lock!), but only for ``task%step <= 1``
  * continue the usual business

This is actually no different than the handling of ``bits%init_nbors`` when
refining, and can be handled inside ``list_t%init_nbors``


.. toctree::
   :maxdepth: 4


