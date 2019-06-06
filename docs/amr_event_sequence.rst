AMR event sequences
-------------------

If threads only ever change the nbor list of their active task, and leave changing
the nbor lists of other tasks to the threads updating those tasks, then much of 
the nbor list locking and copying can be avoided.  The choice is made by not 
setting ``omp_lock%links`` in task%refine

When the active task accessed the nbor list of another task -- this happens 
primarily in ``check_nbors()``, it needs to be sure that the another thread
doesn't change the nbor list in the mean time, so it needs to lock it.  However,
it does not need to lock it's own nbor list during that time.

Each thread needs to lock its own active tasks nbor list in ``init_nbors()``, while
it is being changed, and then needs to lock the nbor's nbor lists when they are
being used (``in check_nbors()``).

Each thread needs to lock it's own task memory while it is being changed (inside
the hd_t%update or in timestep_t%update), and then needs to lock the memory of
nbors as it accesses it (in ``download_t%download_link()``).  

The lines with task memory locking should be marked with "TASK LOCK", while the 
nbor list locking should be marked with "LINK LOCK":

Refining
........

* in make_child_patch()

  * create the new patch
  * add its parent patch as the only nbor, and prolong
  * set bits%init_nbors
  * add it to the task list
  * add it to the ready queue
  * switch perspective to the first time it enters check_current()

* in check_current()

  * detect bits%init_nbors
  * clear bits%init_nbors
  * generate a new nbor list
  * set bits%init_nbors on the nbors (under lock!)
  * continue the usual business

The weak points here are:

* how to trigger check_nbors(), and get nbors into the ready queue

  * this will be ok, since they will get "ready" w/o knowing about the new task
  * it will happen by the normal check_nbors(), done after updates

* how to prevent the bits%init_nbors from spreading beyond the first nbors

  * just look at the task%istep, and do the %set only for tasks with istep<=1

* in task_mesg_t%unpack for virtual tasks
  * must do whatever an active task has to do
  * virtual task automatically do check_nbors(), as part of task_mesg_t%unpack
  * that procedure must also detect bits%init_nbors, and act accordingly

Derefining
...........

* in ``remove_patch()``

  * remove the task from the task list (it will remain in garbage bin while needed)
  * use the existing nbor list to set bits%init_nbors on nbors
  * now that the task is removed, check if nbors have become ready
  * move the task to the garbage bin

* in ``task_mesg_t%npack()``

  * detect bits%init_nbors, inside detection of suicide note
  * set bits%init_nbors on the nbors (under lock!), but only for ``task%step <= 1``
  * continue the usual business

This is actually no different than the handling of ``bits%init_nbors`` when refining,
and can be handled inside ``list_t%init_nbors``

.. toctree::
   :maxdepth: 3

