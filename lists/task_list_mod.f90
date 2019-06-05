!===============================================================================
!> Task list data type, with methods for startup and updates.  Message handling
!> is inherited from the task_mesg_t data type.  At this level, task handling
!> is aware of experiment_t data types, and all data types it has inheritet from.
!===============================================================================
MODULE task_list_mod
  USE io_mod
  USE trace_mod
  USE omp_mod
  USE omp_timer_mod
  USE timer_mod
  USE mpi_mod
  USE mpi_mesg_mod
  USE list_mod
  USE task_mod
  USE patch_mod
  USE experiment_mod
  USE bits_mod
  USE refine_mod
  USE load_balance_mod
  USE task_mesg_mod
  USE timer_mod
  USE validate_mod
  USE download_mod
  implicit none
  private
  type, public, extends(task_mesg_t):: task_list_t
    logical:: syncing=.false.
    logical:: dispatcher=.false.
    real(8):: sync_time=-1.0, sync_next=-1.0
  contains
    procedure:: init
    procedure:: initialize
    procedure:: init_levels
    procedure:: startup
    procedure:: execute
    procedure:: update
    procedure:: average
    procedure:: append_task_list
    procedure:: print => print1
    procedure:: info
    procedure:: init_levelstats
  end type
  public:: task2patch
  type(task_list_t), target:: task_list
CONTAINS

!===============================================================================
!> Initialize the task list, by first initializing the list tasks, then making
!> neighbor lists, and finally checking if they are ready to execute
!===============================================================================
SUBROUTINE init (self, name)
  class(task_list_t):: self
  character(len=*), optional:: name
  class(link_t), pointer:: link
  integer:: iostat
  real(8), save:: job_seconds=1d30
  real(8), save:: sync_time=0.0
  integer, save:: verbose=0
  namelist /task_list_params/ sync_time, job_seconds, verbose
  !-----------------------------------------------------------------------------
  ! An optional namelist can be used to turn debugging on
  !-----------------------------------------------------------------------------
  call trace%begin('task_list_t%init')
  call self%list_t%init (name)
  call self%lock%append
  call mpi_mesg%init
  call task_mesg%init
  rewind (io%input); read(io%input, task_list_params, iostat=iostat)
  if (io%master) write (*,task_list_params)
  if (sync_time > 0.0) &
    self%sync_next = sync_time
  self%sync_time = sync_time
  io%job_seconds = job_seconds
  self%verbose = verbose
  write (io_unit%mpi,*) 'task_list%init: next_sync =', self%sync_next, self%sync_time
  call load_balance%init
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Initialize the task list, so it is ready to execute
!===============================================================================
SUBROUTINE initialize (self)
  class(task_list_t):: self
  !-----------------------------------------------------------------------------
  call self%init_bdries         ! initialize boundaries, based on position
  call self%init_all_nbors      ! make nbor lists, based on position
  call self%reset_status        ! set status buts, based on nbor lists
  call self%count_status        ! count number of different tasks
END SUBROUTINE initialize

!===============================================================================
SUBROUTINE init_levels (self)
  class(task_list_t):: self
  type(link_t), pointer:: link
  !-----------------------------------------------------------------------------
  call trace%begin ('task_list_t%init_levels')
  call refine%init                      ! refinement parameters
  !-----------------------------------------------------------------------------
  ! FIXME: For now, active AMR cannot use recv_priv()
  !-----------------------------------------------------------------------------
  if (refine%on)  then
    !$omp atomic write
    mpi_mesg%recv_priv = .false.
  end if
  !$omp atomic write
  refine%levelmin = refine%levelmax     ! ignore input value!
  link => self%head
  do while (associated(link))
    associate (patch => link%task)
    select type (patch)
    class is (patch_t)
      !-------------------------------------------------------------------------
      ! Set up the refinement "bill board", where requests for refinement are made
      !-------------------------------------------------------------------------
      if (refine%on) then
        allocate (patch%xrefine(patch%gn(1),patch%gn(2),patch%gn(3)))
        call io%bits_mem (storage_size(patch%xrefine), &
                         product(shape(patch%xrefine)), 'xrefine')
        patch%xrefine = -1
      end if
      patch%level = minval(nint(log(patch%box/patch%ds)/log(real(refine%ratio))), &
        mask=patch%n > 1)
    end select
    refine%levelmin = min(refine%levelmin,patch%level)
    refine%levelmax = max(refine%levelmax,patch%level)
    link => link%next
    end associate
  end do
  if (refine%levelmax > refine%levelmin) &
    download%check_filled = .true.
  write (io%output,*) &
    'task_list_t%init_levels: levelmin,max =', refine%levelmin, refine%levelmax
  call self%init_levelstats
  call trace%end
END SUBROUTINE init_levels

!===============================================================================
!> Execute the task list, updating it until it is empty.  With !$omp parallel here,
!> everything local to self%update is thread private.
!===============================================================================
SUBROUTINE startup (self)
  class(task_list_t):: self
  !.............................................................................
  class(link_t), pointer:: link, rem
  class(task_t), pointer:: task
  real(8):: sec
  integer:: dims(4)
  integer, save:: itimer=0
  logical:: was_refined
  !-----------------------------------------------------------------------------
  call trace%begin ('task_list_t%startup', itimer=itimer)
  if (.not.self%dispatcher) call mpi%abort ("no dispatcher is used to execute the task list. Execution is unpredictable, therefore aborting.")
  !-----------------------------------------------------------------------------
  ! Initialize nbor lists and setting virtual and boundary bit should NOT be
  ! done here, since at this point the task list may be a splice between the
  ! task lists of several components, with nbor relations that should not be
  ! overwritten.  Each component should return a task list with appropriate
  ! nbor lists, which may then be combined into a more comprehensive setup, with
  ! special handling of nbor relations
  !-----------------------------------------------------------------------------
  !call self%init_bdries
  !call self%init_all_nbors
  !call self%reset_status
  call self%info
  call self%init_levels
  call validate%init
  !-----------------------------------------------------------------------------
  ! When restarting, the virtual tasks need to be received and unpacked before
  ! the first guard zone downloads.  We therefore trigger a send of all bndry
  ! patches, initialize the recv mechanism, and check for corresponding incoming
  ! messages.
  !-----------------------------------------------------------------------------
  if (mpi_mesg%recv_active .and. mpi_mesg%recv_priv) &
    call self%init_virtual ()
  if (io%restart >= 0) then
    call self%resend_bdry ()
  end if
  !-----------------------------------------------------------------------------
  ! Check the ready status of all tasks.  At this point, virtual patches that
  ! have not been received yet should have negative initial times, so should
  ! prevent boundary tasks from trying to update.
  !-----------------------------------------------------------------------------
  link => self%head
  do while (associated(link))
    call link%task%clear(bits%ready)
    if (link%task%is_set(bits%virtual)) &
      link%task%wc_last = wallclock()
    call self%check_ready (link)
    if (self%sync_time > 0.0) &
      link%task%sync_time = link%task%time + self%sync_time
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! Initialize the data I/O, and start timer
  !-----------------------------------------------------------------------------
  if (mpi%master) &
    print '(a,f8.3,a)', ' Memory per process:', io%gb, ' GB'
  call mpi%barrier ('self%execute')
  call tic (time=sec)
  call timer%print()
  self%n_tasks = self%na
  call trace%end (itimer)
END SUBROUTINE startup

!===============================================================================
!> Execute the task list, updating it until it is empty.  With !$omp parallel here,
!> everything local to self%update is thread private.
!===============================================================================
SUBROUTINE execute (self)
  USE omp_timer_mod
  class(task_list_t):: self
  real(8):: sec
  integer:: dims(4)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('task_list_t%execute', itimer=itimer)
  call self%startup
  call tic (time=sec)
  call timer%print()
  !-----------------------------------------------------------------------------
  ! Loop over tasks
  !-----------------------------------------------------------------------------
  !$omp parallel
  do while (self%na > 0 .and. .not.io%do_stop)
    call self%update (self%head, .false.)
  end do
  call timer%print()
  call mpi_mesg%diagnostics(1)
  write (io_unit%log,*) 'thread',omp%thread,' arrived'
  flush (io_unit%log)
  !$omp end parallel
  call toc ('wall time', timer%n_update, time=sec)
  call mpi%barrier ('end')
  !if (io%master) &
  !  write (io_unit%log,*) 'download_differ % =', 100.*download%n_differ/ &
  !    max(1,(download%n_differ+download%n_same))
  call trace%end (itimer)
END SUBROUTINE execute

!===============================================================================
!> Common part of task list update, used by dispatchers
!===============================================================================
SUBROUTINE update (self, head, test, was_refined, was_derefined)
  class(task_list_t):: self
  class(link_t), pointer:: head
  logical:: test
  logical, optional:: was_refined, was_derefined
  !.............................................................................
  class(task_t), pointer:: task
  real(8):: wc
  logical:: refined, derefined
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  ! Check for flag files.  If io%out_time is set, set a new out_next.
  ! If task%time is slightly smaller than a multiple of out_time, then
  ! out_next will becoe that multiple.
  !----------------------------------------------------------------------------
  call trace%begin('task_list_t%update', itimer=itimer)
  task => head%task
  call io%check_flags
  if (task%id==io%id_track) then
    task%track =  .not.task%track
    io%id_track = 0
  end if
  !---------------------------------------------------------------------------
  ! If the task is frozen, make sure it is marked "not ready", and then return
  !---------------------------------------------------------------------------
  if (task%is_set (bits%frozen)) then
    call task%clear (bits%ready)
    call trace%end (itimer)
    return
  end if
  !-----------------------------------------------------------------------------
  ! Check if refinement is needed on the task; if so this will push new tasks
  ! onto the queue, with the same task time; i.e., to the head of the queue.
  ! If the task is virtual, it will not be checked by this rank.
  !-----------------------------------------------------------------------------
  call refine%check_current(self, head, refined, derefined)
  !-----------------------------------------------------------------------------
  if (present(was_refined)) then
    was_refined = refined
  end if
  if (present(was_derefined)) then
    was_derefined = derefined
  end if
  !-----------------------------------------------------------------------------
  ! If the task was derefined it does not exist, so bail out
  !-----------------------------------------------------------------------------
  if (derefined) then
    call trace%end (itimer)
    return
  end if
  !-----------------------------------------------------------------------------
  ! Download nbor info -- this may be an experiment_t procedure, but if not,
  !-- the call should be answered by a solver-specific procedure (which may or
  ! may not choose to call the generic patch guard zone handler in tasks/)
  !-----------------------------------------------------------------------------
  if (.not.test) &
    call task%dnload                                  ! download nbor info
  !-----------------------------------------------------------------------------
  ! Update the task, whatever that means (may include call to task%output)
  !-----------------------------------------------------------------------------
  if (mpi_mesg%delay_ms==0.0) then
    wc = wallclock()
  end if
  task%rotated = .false.
  if (test) then
    call task%test_update
  else
    select type (task)
    class is (patch_t)
    call validate%check (head, task%mem(:,:,:,task%idx%d,task%it,1), 'before update')
    end select
    !---------------------------------------------------------------------------
    ! IMPORTANT:  The globally available io%ntask must reflect the total number
    ! of active tasks on the rank (io%nwrite is set in data_io_mod)
    !---------------------------------------------------------------------------
    !$omp atomic write
    io%ntask  = self%na                                        ! Authoritative !
    !---------------------------------------------------------------------------
    call task%update                                  ! update the task
    if (self%verbose > 0) then
      associate (unit => merge (io_unit%log, io_unit%mpi, self%verbose > 1))
      write (unit,'(f12.6,2x,a,i6,1p,g14.6)') &
        wallclock(), 'update:', task%id, task%time+task%dtime
      flush (unit)
      end associate
    end if
  end if
  !----------------------------------------------------------------------------
  ! Default delay time = first update time
  !----------------------------------------------------------------------------
  if (mpi_mesg%delay_ms==0.0) then
    !$omp atomic write
    mpi_mesg%delay_ms = 1d3*(wallclock()-wc)
  end if
  !$omp atomic
  mpi_mesg%n_update = mpi_mesg%n_update+1
  if (io%verbose>1) &
    write (io_unit%log,'(a,i4,2x,a,i7,2x,a,2x,a,i7,2x,a,1p,2g14.6,2x,a,2i5,l3)') &
    'thread', omp_mythread, 'task', task%id, trim(task%type), &
    'step', task%istep, 'dt, time:', task%dtime, task%time, &
    'n, nq', self%n, self%nq, associated(head%nbor)
  !-----------------------------------------------------------------------------
  ! The rotate procedure in the generic patch_t module is responsible for
  ! updating the time, and rotating the memory slots where information about
  ! previous time steps are saved.  Task may signal that they have already done
  ! the rotate internally, or else are not yet ready to do so, by setting the
  ! task%rotated flag
  !-----------------------------------------------------------------------------
  if (.not.task%rotated) then
    call task%rotate
  end if
  select type (task)
  class is (patch_t)
  call validate%check (head, task%mem(:,:,:,task%idx%d,task%it,1), ' after update')
  end select
  call task%info (self%nq, self%na)                   ! print info on stdout
  if (io%log_sent > 0) then
    !$omp critical (log_sent_cr)
    call mpi_mesg%log_files ()
    write (io_unit%sent,'(f16.6,i4,2x,a,i6,f16.6,l5)') &
      wallclock(), omp%thread, 'upd', task%id, task%time, task%is_set(bits%boundary)
    flush (io_unit%sent)
    !$omp end critical (log_sent_cr)
  end if
  !-----------------------------------------------------------------------------
  ! If the task is a boundary patch, first check if it should be given to
  ! another rank.  A patch that has been given to another rank becomes
  ! a virtual patch, and the load balance procedure has already sent those
  ! of its nbors that became new boundary patches over to relevant ranks,
  ! while the patch itself is sent here, with bits%virtual+swap_reqest set.
  !-----------------------------------------------------------------------------
  if (task%is_set(bits%boundary)) then                ! boundary patch?
    head%task%nq = self%nq                            ! make sure to pass on
    if (load_balance%check_load (head)) then          ! sell?
       call self%count_status                         ! cheap & robust overkill
    end if
    if (task%id == io%id_debug) &
      write(io_unit%mpi,*) &
        'DBG task_list_t%update: calling send_to_vnbors', task%id
    call self%send_to_vnbors (head)                   ! send to virtual nbors
  end if
  !-----------------------------------------------------------------------------
  ! Tasks that have finished are subtracted from the task list count but are
  ! not removed, since their data may be needed by other tasks (including on
  ! other ranks). In principle it should not be necessary to put a finished
  ! task back into the ready queue, since it has a time > end_time, and thus
  ! should not be able to stop another task from being considered ready.
  !-----------------------------------------------------------------------------
  if (task%has_finished()) then                       ! just finished:
    if (timer%dead_mans_hand == 0.0) then
      timer%dead_mans_hand = wallclock()              ! cf. dispatcher
    end if
    call load_balance%active (.false.)                ! turn off load balancing
    !---------------------------------------------------------------------------
    ! Use bits%frozen to make sure the task only decrements from self%na once
    !---------------------------------------------------------------------------
    if (task%is_clear (bits%frozen)) then
      call task%set (bits%frozen)
      call self%count_status ('has_finished')         ! TEST
      if (io%verbose > 0) then
        write (io_unit%log,*) task%id, 'has finished, na =', self%na
        flush (io_unit%log)
      end if
    end if
  end if
  if (task%time > io%end_time .and. task%is_clear(bits%frozen)) then
    print *, mpi%rank, omp%thread, 'ERROR: task should be frozen', &
      task%time, io%end_time, task%is_set(bits%frozen), task%has_finished()
    call task%set (bits%frozen)
  end if
  !-----------------------------------------------------------------------------
  ! Periodic task syncronization:  The first time a task%time is exactly equal
  ! to load_balance%sync_next, we wait for all ranks to arrive at the same
  ! time. All tasks within a rank are then also at this time, since as long as
  ! one has not arrived there yet, the head%ask%time is smaller.
  !-----------------------------------------------------------------------------
  if (self%syncing) then
    write (io_unit%mpi,*) task%id,omp%thread, &
      ' is waiting on a sync at t =', self%sync_next
    !call trace%end; call trace%begin ('mpi%barrier', itimer=itimer)
    call mpi%barrier ('sync')
    !call trace%end (itimer); call trace%begin('dispatcher0_t%update')
    write (io_unit%mpi,*) task%id,omp%thread, &
      ' finished wating on a sync at t =', self%sync_next
    self%sync_next = self%sync_next + self%sync_time
    self%syncing = .false.
  end if
  task%sync_time = self%sync_next
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Average over a variable with index idx
!===============================================================================
FUNCTION average (self, idx, time)
  class (task_list_t):: self
  integer:: idx
  real(8):: average, time
  !.............................................................................
  class(link_t) , pointer:: link
  class(task_t) , pointer:: task
  real(8), dimension(:,:,:), allocatable:: buffer
  integer:: n, m(3), jt(2), l(3), u(3)
  real:: pt(2)
  !-----------------------------------------------------------------------------
  ! Compute smallest time on the list
  !-----------------------------------------------------------------------------
  call trace%begin ('task_list_t%average')
  link => self%head
  time = link%task%time
  do while (associated(link))
    time = min(time, link%task%time)
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! Interpolate to that time and sum up averages
  !-----------------------------------------------------------------------------
  link => self%head
  n = 0
  average = 0d0
  do while (associated(link))
    task => link%task
    select type (task)
    class is (patch_t)
      m = task%ncell
      l = task%li
      u = l + m - 1
      call task%time_interval (time, jt, pt)
      if (.not.allocated(buffer)) &
        allocate (buffer(m(1),m(2),m(3)))
      buffer = task%mem(l(1):u(1),l(2):u(2),l(3):u(3),idx,jt(1),1)*pt(1) &
             + task%mem(l(1):u(1),l(2):u(2),l(3):u(3),idx,jt(2),1)*pt(2)
      average = average + sum(buffer)
      n = n+1
      !print '(5i5,2f8.5,3f12.6,1p,e12.4)', idx, n, task%id, jt, pt, task%t(jt(1)), time, task%t(jt(2)), average
    end select
    link => link%next
  end do
  average = average/n/product(shape(buffer))
  call trace%end()
END FUNCTION average

!===============================================================================
!> Append a task list to self
!===============================================================================
SUBROUTINE append_task_list (self, task_list)
  class (task_list_t):: self
  class (task_list_t), pointer:: task_list
  class (list_t), pointer:: list
  !-----------------------------------------------------------------------------
  call trace%begin ('task_list_t%append_task_list')
  list => task_list
  call self%append_list (list)
  call trace%end()
END SUBROUTINE append_task_list

!===============================================================================
!> Print a table with task id and all neighbors
!===============================================================================
SUBROUTINE print1 (self, label)
  class (task_list_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  call trace%begin ('task_list_t%print_tasklist')
  if (io%master.and.present(label)) write (io_unit%output,*) label
  link => self%head
  do while (associated(link))
    call link%print_nbors ('')
    link => link%next
  end do
  call self%print_queue
  call trace%end()
END SUBROUTINE print1

!===============================================================================
!> Dump nbor relations, flags and status bits around io%id_debug
!===============================================================================
SUBROUTINE info (self)
  class(task_list_t):: self
  class(link_t), pointer:: link
  class(patch_t), pointer:: task
  !-----------------------------------------------------------------------------
  link => self%head
  do while (associated(link))
    task => task2patch(link%task)
    call task%nbors_info()
    link => link%next
  end do
END SUBROUTINE info

!===============================================================================
!> Upgrade task pointer to experiment level
!===============================================================================
FUNCTION task2experiment (task) RESULT (exper)
  class(task_t), pointer:: task
  class(experiment_t), pointer:: exper
  !...........................................................................
  select type (task)
  class is (experiment_t)
  exper => task
  end select
END FUNCTION task2experiment

!===============================================================================
!> Initialize level statistics, making sure to unly do it once.  Note that, by
!> placing the allocation of levelcost last, threads are encouraged to get
!> caught on this critical region until everything is ready.
!===============================================================================
SUBROUTINE init_levelstats (self)
  class(task_list_t):: self
  class(link_t), pointer:: link
  !-----------------------------------------------------------------------------
  call trace%begin ('task_list_t%init_levelstats')
  !$omp critical (levelstat_cr)
  if (.not. allocated(timer%levelcost)) then
    allocate (timer%levelpop  (refine%levelmin:refine%levelmax))
    timer%levelmin = refine%levelmin
    timer%levelmax = refine%levelmax
    timer%levelpop = 0
    link => self%head
    do while (associated(link))
      !$omp atomic update
      timer%levelpop(link%task%level) = &
      timer%levelpop(link%task%level) + 1
      link => link%next
    end do
    allocate (timer%levelcost(refine%levelmin:refine%levelmax))
    timer%levelcost = 0d0
    if (io%verbose > 1) &
      write (io_unit%output,*) &
        'init_levelstats: levelpop =', &
        refine%levelmin, refine%levelmax, timer%levelpop
  end if
  !$omp end critical (levelstat_cr)
  call trace%end
END SUBROUTINE init_levelstats

END MODULE task_list_mod
