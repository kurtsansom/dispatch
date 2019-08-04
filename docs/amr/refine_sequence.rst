Refine event sequence
--------------------

* In make_child_patch()

  * create the new patch
  * add its parent patch as the only nbor, and prolong
  * set bits%init_nbors
  * add it to the task list
  * add it to the ready queue
  * switch perspective to the first time it enters check_current()

* In check_current()

  * detect bits%init_nbors
  * clear bits%init_nbors
  * generate a new nbor list
  * set bits%init_nbors on the nbors (under lock!)
  * continue the usual business

The weak points here are:

* how to trigger ``check_nbors()``, and get nbors into the ready queue

  * this will be ok, since they will get "ready" w/o knowing about the new task
  * it will happen by the normal check_nbors(), done after updates

* how to prevent the bits%init_nbors from spreading beyond the first nbors

  * just look at the task%istep, and do the %set only for tasks with istep<=1

* in ``task_mesg_t%unpack()`` for virtual tasks

  * must do whatever an active task has to do
  * virtual task automatically do check_nbors(), as part of task_mesg_t%unpack
  * that procedure must also detect bits%init_nbors, and act accordingly


.. toctree::
   :maxdepth: 4


