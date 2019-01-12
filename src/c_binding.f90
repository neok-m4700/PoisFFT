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
#define RPC c_double
   subroutine poisfft_solver1d_new(self, nxyz, lxyz, bcs, approximation, gnxyz, offs, comm_ptr, nthreads) bind(c, name="poisfft_solver1d_new")
#define DIM 1
#define SOLVER poisfft_solver1d_dp
#include "c_new-inc.f90"
#undef SOLVER
#undef DIM
   end subroutine

   subroutine poisfft_solver2d_new(self, nxyz, lxyz, bcs, approximation, gnxyz, offs, comm_ptr, nthreads) bind(c, name="poisfft_solver2d_new")
#define DIM 2
#define SOLVER poisfft_solver2d_dp
#include "c_new-inc.f90"
#undef SOLVER
#undef DIM
   end subroutine

   subroutine poisfft_solver3d_new(self, nxyz, lxyz, bcs, approximation, gnxyz, offs, comm_ptr, nthreads) bind(c, name="poisfft_solver3d_new")
#define DIM 3
#define SOLVER poisfft_solver3d_dp
#include "c_new-inc.f90"
#undef SOLVER
#undef DIM
   end subroutine

#undef RPC
#define RPC c_float
   subroutine poisfft_solver1d_f_new(self, nxyz, lxyz, bcs, approximation, gnxyz, offs, comm_ptr, nthreads) bind(c, name="poisfft_solver1d_f_new")
#define DIM 1
#define SOLVER poisfft_solver1d_sp
#include "c_new-inc.f90"
#undef SOLVER
#undef DIM
   end subroutine

   subroutine poisfft_solver2d_f_new(self, nxyz, lxyz, bcs, approximation, gnxyz, offs, comm_ptr, nthreads) bind(c, name="poisfft_solver2d_f_new")
#define DIM 2
#define SOLVER poisfft_solver2d_sp
#include "c_new-inc.f90"
#undef SOLVER
#undef DIM
   end subroutine

   subroutine poisfft_solver3d_f_new(self, nxyz, lxyz, bcs, approximation, gnxyz, offs, comm_ptr, nthreads) bind(c, name="poisfft_solver3d_f_new")
#define DIM 3
#define SOLVER poisfft_solver3d_sp
#include "c_new-inc.f90"
#undef SOLVER
#undef DIM
   end subroutine

#undef RPC
#define RPC c_double
   subroutine poisfft_solver1d_execute(self, phi, rhs, ngphi, ngrhs) bind(c, name="poisfft_solver1d_execute")
      type(poisfft_solver1d_dp), pointer :: f_self
#define DIM 1
#include "c_execute-inc.f90"
#undef DIM
   end subroutine

   subroutine poisfft_solver2d_execute(self, phi, rhs, ngphi, ngrhs) bind(c, name="poisfft_solver2d_execute")
      type(poisfft_solver2d_dp), pointer :: f_self
#define DIM 2
#include "c_execute-inc.f90"
#undef DIM
   end subroutine

   subroutine poisfft_solver3d_execute(self, phi, rhs, ngphi, ngrhs) bind(c, name="poisfft_solver3d_execute")
      type(poisfft_solver3d_dp), pointer :: f_self
#define DIM 3
#include "c_execute-inc.f90"
#undef DIM
   end subroutine

#undef RPC
#define RPC c_float
   subroutine poisfft_solver1d_f_execute(self, phi, rhs, ngphi, ngrhs) bind(c, name="poisfft_solver1d_f_execute")
      type(poisfft_solver1d_sp), pointer :: f_self
#define DIM 1
#include "c_execute-inc.f90"
#undef DIM
   end subroutine

   subroutine poisfft_solver2d_f_execute(self, phi, rhs, ngphi, ngrhs) bind(c, name="poisfft_solver2d_f_execute")
      type(poisfft_solver2d_sp), pointer :: f_self
#define DIM 2
#include "c_execute-inc.f90"
#undef DIM
   end subroutine

   subroutine poisfft_solver3d_f_execute(self, phi, rhs, ngphi, ngrhs) bind(c, name="poisfft_solver3d_f_execute")
      type(poisfft_solver3d_sp), pointer :: f_self
#define DIM 3
#include "c_execute-inc.f90"
#undef DIM
   end subroutine

#undef RPC
   subroutine poisfft_solver1d_finalize(self) bind(c, name="poisfft_solver1d_finalize")
      type(poisfft_solver1d_dp), pointer :: f_self
#include "c_finalize-inc.f90"
   end subroutine

   subroutine poisfft_solver2d_finalize(self) bind(c, name="poisfft_solver2d_finalize")
      type(poisfft_solver2d_dp), pointer :: f_self
#include "c_finalize-inc.f90"
   end subroutine

   subroutine poisfft_solver3d_finalize(self) bind(c, name="poisfft_solver3d_finalize")
      type(poisfft_solver3d_dp), pointer :: f_self
#include "c_finalize-inc.f90"
   end subroutine

   subroutine poisfft_solver1d_f_finalize(self) bind(c, name="poisfft_solver1d_f_finalize")
      type(poisfft_solver1d_sp), pointer :: f_self
#include "c_finalize-inc.f90"
   end subroutine

   subroutine poisfft_solver2d_f_finalize(self) bind(c, name="poisfft_solver2d_f_finalize")
      type(poisfft_solver2d_sp), pointer :: f_self
#include "c_finalize-inc.f90"
   end subroutine

   subroutine poisfft_solver3d_f_finalize(self) bind(c, name="poisfft_solver3d_f_finalize")
      type(poisfft_solver3d_sp), pointer :: f_self
#include "c_finalize-inc.f90"
   end subroutine

end module
