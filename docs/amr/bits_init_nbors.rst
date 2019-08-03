Handling ``bits%init_nbors``
----------------------------

When a thread updating a boundary task discovers a ``bits%init_nbors``,
it should not clear that bit until after the ``send_to_vnbors()`` call.

When should the thread create the new nbor list?  The thread is updating the
task because it was deemed ready to update based on the existing nbor list,
which has then been used (before the call to ``task%update()``) to load data
into the guard zones.  The existing nbor list might also be used for some
unspecified reason as part of the update, so it should not be touched until
after the task has updated.

Hence the ``init_nbors()`` call that should be the result of detecting that
bit should *not* be made by the unpack procedure; that would violate the
principle of identical action.  Instead, it should be done as part of the
task handling, presumable by ``task_list_t%update()``, and *after* its call
to ``task%update()``.

It *must* be done before the call to ``load_balance()``, because otherwise the
task might be passed on to another rank, where its nbor list will be generated
from scratch during unpack.

The state of that task after the ``task%update()`` and ``init_nbors()`` is that
of a "dormant" task, ready to be checked as a candidate for the ready queue by 
any thread updating one of its (possibly new) nbors.  That is at it should be.


.. toctree::
   :maxdepth: 4
