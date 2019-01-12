type(solver), pointer :: f_D
type(c_ptr), intent(out) :: D
integer(c_int) :: nxyz(dims), BCs(2 * dims)
real(rp) :: Lxyz(dims)
integer(c_int), optional :: approximation, gnxyz(dims), offs(dims)
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

allocate(f_d)

if (nthreads < 1) nthreads = 1

if (present(approximation)) then
   appr = int(approximation)
else
   appr = 0
end if

if (present(gnxyz) .and. present(offs)) then
   f_d = solver(int(nxyz(dims:1:-1)), &
      lxyz(dims:1:-1), &
      int(bcs(idxs(2 * dims:1:-1))), &
      appr, &
      int(gnxyz(dims:1:-1)), &
      int(offs(dims:1:-1)), &
      f_comm, &
      int(nthreads))
else
   f_d = solver(int(nxyz(dims:1:-1)), &
      lxyz(dims:1:-1), &
      int(bcs(idxs(2 * dims:1:-1))), &
      appr, &
      nthreads = int(nthreads))
end if

d = c_loc(f_d)
