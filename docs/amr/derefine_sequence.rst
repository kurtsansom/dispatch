Derefine event sequence
------------------------

* Call ``list_t%remove_remove_and_reset()`` to

  * remove the task from the task list (it will remain in garbage bin while
    needed)
  * use the existing nbor list to set ``bits%init_nbors`` on nbors
  * now that the task is removed, check if nbors have become ready
  * move the task to the garbage bin
  * when the task is no longer needed (no longer is a member of an nbor list)
    the task is deallocated and deleted by ``list_t%garbage_collect()``

* In ``task_mesg_t%unpack()``

  * detect ``bits%remove`` and initiate the removal
  * remove the virtual task from the task list add it to the garbage bin
  * call ``set_init_nbors()``, which sets ``bits%init_nbors`` in the tasks
    in the nbor list of the task
  * continue the usual business
  * when the nbor lists of all of the tasks nbors have been remade, the
    number of references to the task reaches zero, and the task is deallocated
    and deleted by the ``garbage_collect()`` procedure

Part of this is identical to what is done on the owner rank, and is taken
care of by the ``list_t%remove_and_reset()`` procedure.


.. toctree::
   :maxdepth: 4


