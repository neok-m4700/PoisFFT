module poisfft
#ifdef MPI
   use pfft
#endif
   use poisfft_constants
   use poisfft_sp, &
      poisfft_solver1d_sp => poisfft_solver1d, &
      poisfft_solver2d_sp => poisfft_solver2d, &
      poisfft_solver3d_sp => poisfft_solver3d
   use poisfft_dp, &
      poisfft_solver1d_dp => poisfft_solver1d, &
      poisfft_solver2d_dp => poisfft_solver2d, &
      poisfft_solver3d_dp => poisfft_solver3d
   implicit none

contains
#ifdef MPI
   subroutine poisfft_initmpigrid(mpi_comm, np, poisfft_comm, ierr)
      use mpi, only : mpi_comm_size
      integer, intent(in) :: mpi_comm
      integer(c_int), intent(in) :: np(:)
      integer, intent(out) :: poisfft_comm, ierr
      integer :: comm_size

      call mpi_comm_size(mpi_comm, comm_size, ierr)
      if (ierr /= 0 .or. comm_size /= product(np)) then
         write(ferr, *) '[poisfft_initmpigrid] ERROR: ierr=', ierr, ' comm_size=', comm_size, ' product(np)=', product(np)
         ierr = 1; return
      end if
      ! this will fail to compile if integer /= 4 bytes as pfft assumes.
      ierr = pfft_create_procmesh_2d(mpi_comm, np(1), np(2), poisfft_comm)
   end subroutine

   subroutine poisfft_init(mpi_comm, dim, np, poisfft_comm, ierr)
      use mpi, only : mpi_comm_size
      integer, intent(in) :: mpi_comm, dim
      integer(c_int), intent(in) :: np(:)
      integer, intent(out) :: poisfft_comm, ierr
      integer :: comm_size

      call mpi_comm_size(mpi_comm, comm_size, ierr)
      if (ierr /= 0 .or. comm_size /= product(np)) then
         write(ferr, *) '[poisfft_init] ERROR: ierr=', ierr, ' comm_size=', comm_size, ' product(np)=', product(np)
         ierr = 1; return
      end if
      ierr = pfft_create_procmesh(dim, mpi_comm, np, poisfft_comm)
   end subroutine

   subroutine poisfft_localgridsize(rnk, nos, poisfft_comm, local_ni, local_i_start, local_no, local_o_start)
      integer(c_int), value :: rnk
      integer(c_intptr_t), intent(in) :: nos(:)
      integer, intent(in) :: poisfft_comm
      integer(c_intptr_t), intent(out) :: local_ni(:), local_i_start(:)
      integer(c_intptr_t), intent(out) :: local_no(:), local_o_start(:)
      integer(c_intptr_t) :: alloc_local, nback(1:size(nos))

      nback = nos(rnk:1:-1) ! not necessary, but to avoid warnings about temporary array passed
      alloc_local = pfft_local_size_dft(rnk, nback, poisfft_comm, pfft_transposed_none, local_ni, local_i_start, local_no, local_o_start)
      ! alloc_local = pfft_local_size_dft_3d(nback, poisfft_comm, pfft_transposed_none, local_ni, local_i_start, local_no, local_o_start)
      local_ni = local_ni(rnk:1:-1)
      local_i_start = local_i_start(rnk:1:-1)
      local_no = local_no(rnk:1:-1)
      local_o_start = local_o_start(rnk:1:-1)
   end subroutine

#endif
end module
