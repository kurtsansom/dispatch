Send procedures
----------------

Send requests are started with ``MPI_ISEND()`` calls, and the
requests are held on lists until each set of sends (each task is in general
sent to multiple ranks) is completed (as tested via ``MPI_TEST_ALL()``), at
which time the send buffer is deallocated and the send request is removed
from the ``sent_list`` and deleted.

It is an advantage to arrange this so that each thread maintain a thread-
private ``sent_list`` -- this avoids the need for OMP critical regions or
locks.

.. toctree::
   :maxdepth: 4

