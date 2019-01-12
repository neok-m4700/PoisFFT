type(SOLVER), pointer :: f_self
type(c_ptr), intent(out) :: self
integer(c_int) :: nxyz(DIM), bcs(2 * DIM)
real(RPC) :: lxyz(DIM)
integer(c_int), optional :: approximation, gnxyz(DIM), offs(DIM)
type(c_ptr), value :: comm_ptr
integer(c_int), value :: nthreads
integer :: f_comm, appr
integer, parameter :: idxs(6) = [2, 1, 4, 3, 6, 5]

#ifdef MPI
f_comm = mpi_comm_c2f(comm_ptr)
#else
 !to avoid unused variable warning
if (c_associated(comm_ptr)) f_comm = 0
#endif

allocate(f_self)

if (nthreads < 1) nthreads = 1

if (present(approximation)) then
   appr = int(approximation)
else
   appr = 0
end if

if (present(gnxyz) .and. present(offs)) then
   f_self = SOLVER(&
      int(nxyz(DIM:1:-1)), lxyz(DIM:1:-1), int(bcs(idxs(2 * DIM:1:-1))), &
      appr, int(gnxyz(DIM:1:-1)), int(offs(DIM:1:-1)), f_comm, int(nthreads))
else
   f_self = SOLVER(int(nxyz(DIM:1:-1)), lxyz(DIM:1:-1), int(bcs(idxs(2 * DIM:1:-1))), appr, nthreads=int(nthreads))
end if

self = c_loc(f_self)
