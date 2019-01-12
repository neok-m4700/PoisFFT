module poisfft_c_binding
   use iso_c_binding
   use poisfft

   implicit none

#ifdef MPI
   interface
      integer function mpi_comm_c2f(c_handle) bind(c, name="f_MPI_Comm_c2f")
         use iso_c_binding
         type(c_ptr), value :: c_handle
      end function
   end interface

#endif

contains
#define rp c_double
   subroutine poisfft_solver1d_new(d, nxyz, lxyz, bcs, approximation, gnxyz, offs, comm_ptr, nthreads) bind(c, name="poisfft_solver1d_new")
#define dims 1
#define solver poisfft_solver1d_dp
#include "c_new-inc.f90"
#undef solver
#undef dims
   end subroutine

   subroutine poisfft_solver2d_new(d, nxyz, lxyz, bcs, approximation, gnxyz, offs, comm_ptr, nthreads) bind(c, name="poisfft_solver2d_new")
#define dims 2
#define solver poisfft_solver2d_dp
#include "c_new-inc.f90"
#undef solver
#undef dims
   end subroutine

   subroutine poisfft_solver3d_new(d, nxyz, lxyz, bcs, approximation, gnxyz, offs, comm_ptr, nthreads) bind(c, name="poisfft_solver3d_new")
#define dims 3
#define solver poisfft_solver3d_dp
#include "c_new-inc.f90"
#undef solver
#undef dims
   end subroutine
#undef rp

#define rp c_float
   subroutine poisfft_solver1d_f_new(d, nxyz, lxyz, bcs, approximation, gnxyz, offs, comm_ptr, nthreads) bind(c, name="poisfft_solver1d_f_new")
#define dims 1
#define solver poisfft_solver1d_sp
#include "c_new-inc.f90"
#undef solver
#undef dims
   end subroutine

   subroutine poisfft_solver2d_f_new(d, nxyz, lxyz, bcs, approximation, gnxyz, offs, comm_ptr, nthreads) bind(c, name="poisfft_solver2d_f_new")
#define dims 2
#define solver poisfft_solver2d_sp
#include "c_new-inc.f90"
#undef solver
#undef dims
   end subroutine

   subroutine poisfft_solver3d_f_new(d, nxyz, lxyz, bcs, approximation, gnxyz, offs, comm_ptr, nthreads) bind(c, name="poisfft_solver3d_f_new")
#define dims 3
#define solver poisfft_solver3d_sp
#include "c_new-inc.f90"
#undef solver
#undef dims
   end subroutine

#undef rp

#define rp c_double
   subroutine poisfft_solver1d_execute(d, phi, rhs, ngphi, ngrhs) bind(c, name="poisfft_solver1d_execute")
      type(poisfft_solver1d_dp), pointer :: f_d
#define dims 1
#include "c_execute-inc.f90"
#undef dims
   end subroutine

   subroutine poisfft_solver2d_execute(d, phi, rhs, ngphi, ngrhs) bind(c, name="poisfft_solver2d_execute")
      type(poisfft_solver2d_dp), pointer :: f_d
#define dims 2
#include "c_execute-inc.f90"
#undef dims
   end subroutine

   subroutine poisfft_solver3d_execute(d, phi, rhs, ngphi, ngrhs) bind(c, name="poisfft_solver3d_execute")
      type(poisfft_solver3d_dp), pointer :: f_d
#define dims 3
#include "c_execute-inc.f90"
#undef dims
   end subroutine

#undef dims
#undef rp

#define rp c_float
   subroutine poisfft_solver1d_f_execute(d, phi, rhs, ngphi, ngrhs) bind(c, name="poisfft_solver1d_f_execute")
      type(poisfft_solver1d_sp), pointer :: f_d
#define dims 1
#include "c_execute-inc.f90"
#undef dims
   end subroutine

   subroutine poisfft_solver2d_f_execute(d, phi, rhs, ngphi, ngrhs) bind(c, name="poisfft_solver2d_f_execute")
      type(poisfft_solver2d_sp), pointer :: f_d
#define dims 2
#include "c_execute-inc.f90"
#undef dims
   end subroutine

   subroutine poisfft_solver3d_f_execute(d, phi, rhs, ngphi, ngrhs) bind(c, name="poisfft_solver3d_f_execute")
      type(poisfft_solver3d_sp), pointer :: f_d
#define dims 3
#include "c_execute-inc.f90"
#undef dims
   end subroutine

#undef dims
#undef rp
   subroutine poisfft_solver1d_finalize(d) bind(c, name="poisfft_solver1d_finalize")
      type(poisfft_solver1d_dp), pointer :: f_d
#include "c_finalize-inc.f90"
   end subroutine

   subroutine poisfft_solver2d_finalize(d) bind(c, name="poisfft_solver2d_finalize")
      type(poisfft_solver2d_dp), pointer :: f_d
#include "c_finalize-inc.f90"
   end subroutine

   subroutine poisfft_solver3d_finalize(d) bind(c, name="poisfft_solver3d_finalize")
      type(poisfft_solver3d_dp), pointer :: f_d
#include "c_finalize-inc.f90"
   end subroutine

   subroutine poisfft_solver1d_f_finalize(d) bind(c, name="poisfft_solver1d_f_finalize")
      type(poisfft_solver1d_sp), pointer :: f_d
#include "c_finalize-inc.f90"
   end subroutine

   subroutine poisfft_solver2d_f_finalize(d) bind(c, name="poisfft_solver2d_f_finalize")
      type(poisfft_solver2d_sp), pointer :: f_d
#include "c_finalize-inc.f90"
   end subroutine

   subroutine poisfft_solver3d_f_finalize(d) bind(c, name="poisfft_solver3d_f_finalize")
      type(poisfft_solver3d_sp), pointer :: f_d
#include "c_finalize-inc.f90"
   end subroutine

end module
