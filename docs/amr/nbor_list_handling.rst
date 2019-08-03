nbor list handling
-------------------

Read-write access
..................

If threads were to be allowed to update the nbor lists of tasks that another
thread might be working on, the situation would become complicated. The thread
would need to lock the nbor list of the other task while updating it, but
it cannot be allowed to lock it's own nbor list during that time, since the 
thread working on that task might happen to be running the same procedure at
the same time, and a deadlock could then occur.

Read-only access
..................

If threads only ever change the nbor list of their active task, and leave changing
the nbor lists of other tasks to the threads updating those tasks, then much of 
the nbor list locking and copying can be avoided.  

When the active task accesses the nbor list of another task -- this happens 
primarily in ``check_nbors()``, it needs to be sure that the another thread
doesn't change the nbor list in the mean time, so it needs to lock it.  That
locking is only effective if a thread that updates its own nbor list actually
locks it when it does so.  Hence even with read-only access, the task owning
thread must lock the nbor list link while changing the nbor list.  Currently
this is primarily done by ``list_t%init_nbors()``.

Conclusions
...........

The simplest and safest nbor list handling approach is based on these
principles:

1) Threads need to lock their nbor lists while updating them.  This has the 
   extra benefit that the updating does not require copying the nbor list in
   the mean time, except for the benefit of making the lock time as short as
   possible.  However, considering how relatively seldom nbor lists need to
   be updated, locking during update is the simplest and safest.

2) While accessing the nbor lists of other tasks (e.g. during ``check_nbors()``),
   threads must acquire a lock on the nbor list while using it, to ensure it
   isn't being changed in the meantime.  However, it must not modify the nbor
   lists of other tasks

3) The nbor lists presented to other tasks need not be completely up-to-date,
   but they must be *valid*, in the sense that the nbors exist and have
   the expected time advance.

4) If at some point in (wall clock) time a task (A) has an nbor list that
   contains an nbor (B) that does *not* have task A in its nbor list, then to
   avoid deadlocks, a ``check_ready()`` on task A must happen as a consequence
   of it being added to the nbor list of task B.


.. toctree::
   :maxdepth: 4
