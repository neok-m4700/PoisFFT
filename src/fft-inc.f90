#ifndef NO_CONTIGUOUS
#define CONTIG , contiguous
#else
#define CONTIG
#endif

#if (PREC==2)
#define RP DRP
#define CP DCP
#define FFTW_EXECUTE_COMPLEX fftw_execute_dft
#define FFTW_MPI_EXECUTE_COMPLEX fftw_mpi_execute_dft
#define FFTW_EXECUTE_REAL fftw_execute_r2r
#define PFFT_EXECUTE_GEN pfft_execute
#else
#define RP SRP
#define CP SCP
#define FFTW_EXECUTE_COMPLEX fftwf_execute_dft
#define FFTW_MPI_EXECUTE_COMPLEX fftwf_mpi_execute_dft
#define FFTW_EXECUTE_REAL fftwf_execute_r2r
#define PFFT_EXECUTE_GEN pfftf_execute
#endif
use iso_c_binding
use poisfft_precisions
use fftw3
use iso_fortran_env, only : &
   fin => input_unit, &
   fout => output_unit, &
   ferr => error_unit
#ifdef MPI
use pfft
#endif
implicit none

integer, parameter :: fft_complex = 0
integer, parameter :: fft_realeven00 = fftw_redft00
integer, parameter :: fft_realeven01 = fftw_redft01
integer, parameter :: fft_realeven10 = fftw_redft10
integer, parameter :: fft_realeven11 = fftw_redft11
integer, parameter :: fft_realodd00 = fftw_rodft00
integer, parameter :: fft_realodd01 = fftw_rodft01
integer, parameter :: fft_realodd10 = fftw_rodft10
integer, parameter :: fft_realodd11 = fftw_rodft11

integer, parameter :: fft_distributed_fftw = 1
integer, parameter :: fft_distributed_pfft = 2

type poisfft_plan1d
   type(c_ptr) :: planptr = c_null_ptr
   logical :: planowner = .false.
   logical :: distributed = .false.
   integer(c_int) :: dir
end type


type poisfft_plan2d
   type(c_ptr) :: planptr = c_null_ptr
   logical :: planowner = .false.
   logical :: distributed = .false.
   integer :: method = fft_distributed_pfft
   integer(c_int) :: dir
end type


type poisfft_plan3d
   type(c_ptr) :: planptr = c_null_ptr
   logical :: planowner = .false.
   logical :: distributed = .false.
   integer(c_int) :: dir
end type


type mpi_vars_1d
   integer :: comm = -1 ! MPI communicator for the exchange
   integer :: rank ! our rank in comm
   integer :: np ! number of processes in comm
   ! s for dimensions of send buffers, r for receive buffers
   ! nx, nz are the dimensions of individual parts which will be sent to or received from other processes in comm
   ! displs are the displacements (offsets) of individual pieces in the whole 1D buffer
   ! counts are the the number of elements in each piece
   ! sumrnzs(i) is the sum of rnzs in pieces 1..i-1
   integer, dimension(:), allocatable :: snxs, snzs, rnxs, rnzs, sdispls, scounts, rdispls, rcounts, sumrnzs
   !contiguous 3D buffer for locally transposed RHS
   real(RP), allocatable :: tmp1(:, :, :)
   !1D buffer for globally transposed blocks of tmp1 from different processes after MPI_Alltoallv
   real(RP), allocatable :: tmp2(:)
   !3D buffer for contiguous 1D rows on which 1D FFT can be performed locally,
   ! constructed by local reordering of tmp2
   real(RP), allocatable :: rwork(:, :, :)
end type

type poisfft_solver1d
   real(RP) :: lx
   integer(c_int) :: nx, gnx, nxyz(1)
   integer(c_int) :: offx = 0 !offset from global index
   integer(c_size_t) :: cnt, gcnt
   real(RP) :: norm_factor
   integer(c_int), dimension(2) :: bcs
   real(RP), allocatable, dimension(:) :: denomx
   integer :: approximation = 0

   type(poisfft_plan1d) :: forward, backward
   complex(CP), dimension(:), pointer CONTIG :: cwork => null()
   real(RP), dimension(:), pointer CONTIG :: rwork => null()

   logical :: mpi_transpose_needed = .false.
   integer :: nthreads = 1
   type(mpi_vars_1D) :: mpi
end type

type mpi_vars_2d
   integer :: rank
   integer :: np
   integer :: comm = -1
   integer :: comm_dim = 2
end type

type poisfft_solver2d
   real(RP) :: lx, ly
   integer(c_int) :: nx, ny, gnx, gny, nxyz(2)
   integer(c_int) :: offx = 0, offy = 0 !offsets from global indexes
   integer(c_size_t) :: cnt, gcnt
   real(RP) :: norm_factor
   integer(c_int), dimension(4) :: bcs
   real(RP), allocatable, dimension(:) :: denomx, denomy
   integer :: approximation = 0

   type(poisfft_plan2d) :: forward, backward
   complex(CP), dimension(:, :), pointer CONTIG :: cwork => null()
   real(RP), dimension(:, :), pointer CONTIG :: rwork => null()
   integer :: nthreads = 1
   type(mpi_vars_2d) :: mpi
end type

type mpi_vars_3d
   integer :: comm = -1
end type

type poisfft_solver3d
   real(RP) :: lx, ly, lz
   integer(c_int) :: nx, ny, nz, gnx, gny, gnz, nxyz(3)
   integer(c_int) :: offx = 0, offy = 0, offz = 0 !offsets from global indexes
   integer(c_size_t) :: cnt, gcnt
   real(RP) :: norm_factor
   integer(c_int), dimension(6) :: bcs
   real(RP), allocatable, dimension(:) :: denomx, denomy, denomz
   integer :: approximation = 0

   type(poisfft_plan3d) :: forward, backward
   complex(CP), dimension(:, :, :), pointer CONTIG :: cwork => null()
   real(RP), dimension(:, :, :), pointer CONTIG :: rwork => null()
   integer :: nthreads = 1
   !will be used in splitting for some boundary conditions
   type(poisfft_solver1d), dimension(:), allocatable :: solvers1d
   type(poisfft_solver2d), dimension(:), allocatable :: solvers2d
   type(mpi_vars_3d) :: mpi
end type

interface deallocate_fftw
   module procedure deallocate_fftw_1d_complex
   module procedure deallocate_fftw_1d_real
   module procedure deallocate_fftw_2d_complex
   module procedure deallocate_fftw_2d_real
   module procedure deallocate_fftw_3d_complex
   module procedure deallocate_fftw_3d_real
end interface

interface allocate_fftw_real
   module procedure allocate_fftw_1d_real
   module procedure allocate_fftw_2d_real
   module procedure allocate_fftw_3d_real
end interface

interface allocate_fftw_complex
   module procedure allocate_fftw_1d_complex
   module procedure allocate_fftw_2d_complex
   module procedure allocate_fftw_3d_complex
end interface

interface execute
   module procedure poisfft_plan1d_execute_real
   module procedure poisfft_plan1d_execute_complex
   module procedure poisfft_plan2d_execute_real
   module procedure poisfft_plan2d_execute_complex
   module procedure poisfft_plan3d_execute_real
   module procedure poisfft_plan3d_execute_complex
end interface

#ifdef MPI

interface execute_mpi
   module procedure poisfft_plan1d_execute_mpi
   module procedure poisfft_plan2d_execute_mpi
   module procedure poisfft_plan3d_execute_mpi
end interface

#endif

interface finalize
   module procedure poisfft_plan1d_finalize
   module procedure poisfft_plan2d_finalize
   module procedure poisfft_plan3d_finalize
end interface

interface poisfft_plan1d
   module procedure poisfft_plan1d_new
end interface

interface poisfft_plan2d
   module procedure poisfft_plan2d_new
end interface

interface poisfft_plan3d
   module procedure poisfft_plan3d_new
end interface

contains

function poisfft_plan1d_new(self, plantypes, distributed) result(plan)
#define dimensions 1
#include "plan_new-inc.f90"
#undef dimensions
end function

function poisfft_plan2d_new(self, plantypes, distributed) result(plan)
#define dimensions 2
#include "plan_new-inc.f90"
#undef dimensions
end function

function poisfft_plan3d_new(self, plantypes, distributed) result(plan)
#define dimensions 3
#include "plan_new-inc.f90"
#undef dimensions
end function

subroutine allocate_fftw_1D_complex(self)
#define dimensions 1
#define realcomplex 2

#include "allocate_fftw-inc.f90"

#undef dimensions
#undef realcomplex
end subroutine

subroutine allocate_fftw_1d_real(self)
#define dimensions 1
#define realcomplex 1

#include "allocate_fftw-inc.f90"

#undef dimensions
#undef realcomplex
end subroutine

subroutine allocate_fftw_2d_complex(self)
#define dimensions 2
#define realcomplex 2

#include "allocate_fftw-inc.f90"

#undef dimensions
#undef realcomplex
end subroutine

subroutine allocate_fftw_2d_real(self)
#define dimensions 2
#define realcomplex 1

#include "allocate_fftw-inc.f90"

#undef dimensions
#undef realcomplex
end subroutine

subroutine allocate_fftw_3d_complex(self)
#define dimensions 3
#define realcomplex 2

#include "allocate_fftw-inc.f90"

#undef dimensions
#undef realcomplex
end subroutine

subroutine allocate_fftw_3d_real(self)
#define dimensions 3
#define realcomplex 1

#include "allocate_fftw-inc.f90"

#undef dimensions
#undef realcomplex
end subroutine

subroutine poisfft_plan1d_execute_complex(plan, data)
   type(poisfft_plan1d), intent(in) :: plan
   complex(CP), dimension(:) CONTIG :: data
   call FFTW_EXECUTE_COMPLEX(plan % planptr, data, data)
end subroutine

subroutine poisfft_plan1d_execute_real(plan, data)
   type(poisfft_plan1d), intent(in) :: plan
   real(RP), dimension(:) CONTIG :: data
   call FFTW_EXECUTE_REAL(plan % planptr, data, data)
end subroutine

subroutine poisfft_plan2d_execute_complex(plan, data)
   type(poisfft_plan2d), intent(in) :: plan
   complex(CP), dimension(:, :) CONTIG :: data
   call FFTW_EXECUTE_COMPLEX(plan % planptr, data, data)
end subroutine

subroutine poisfft_plan2d_execute_real(plan, data)
   type(poisfft_plan2d), intent(in) :: plan
   real(RP), dimension(:, :) CONTIG :: data
   ! write(ferr, *) 'shape(data)=', shape(data)
   call FFTW_EXECUTE_REAL(plan % planptr, data, data)
end subroutine

subroutine poisfft_plan3d_execute_complex(plan, data)
   type(poisfft_plan3d), intent(in) :: plan
   complex(CP), dimension(:, :, :) CONTIG :: data
   call FFTW_EXECUTE_COMPLEX(plan % planptr, data, data)
end subroutine

subroutine poisfft_plan3d_execute_real(plan, data)
   type(poisfft_plan3d), intent(in) :: plan
   real(RP), dimension(:, :, :) CONTIG :: data
   call FFTW_EXECUTE_REAL(plan % planptr, data, data)
end subroutine

#ifdef MPI
subroutine poisfft_plan1d_execute_mpi(plan)
   type(poisfft_plan1d), intent(in) :: plan
   call PFFT_EXECUTE_GEN(plan % planptr)
end subroutine

subroutine poisfft_plan2d_execute_mpi(plan, data)
   type(poisfft_plan2d), intent(in) :: plan
   complex(CP), dimension(:, :), optional CONTIG :: data

   if (plan % method == fft_distributed_fftw) then
      if (present(data)) then
         call FFTW_MPI_EXECUTE_COMPLEX(plan % planptr, data, data)
      else
         stop 'Error, missing `data` argument in a call to Execute_MPI with FFTW.'
      end if
   else
      call PFFT_EXECUTE_GEN(plan % planptr)
   end if
end subroutine

subroutine poisfft_plan3d_execute_mpi(plan)
   type(poisfft_plan3d), intent(in) :: plan
   call PFFT_EXECUTE_GEN(plan % planptr)
end subroutine

#endif

subroutine deallocate_fftw_1d_complex(data)
   complex(CP), dimension(:), pointer CONTIG :: data
   type(c_ptr) :: p
   p = c_loc(data(1)); call fftw_free(p); nullify(data)
end subroutine

subroutine deallocate_fftw_1d_real(data)
   real(RP), dimension(:), pointer CONTIG :: data
   type(c_ptr) :: p
   p = c_loc(data(1)); call fftw_free(p); nullify(data)
end subroutine

subroutine deallocate_fftw_2d_complex(data)
   complex(CP), dimension(:, :), pointer CONTIG :: data
   type(c_ptr) :: p
   p = c_loc(data(1, 1)); call fftw_free(p); nullify(data)
end subroutine

subroutine deallocate_fftw_2d_real(data)
   real(RP), dimension(:, :), pointer CONTIG :: data
   type(c_ptr) :: p
   p = c_loc(data(1, 1)); call fftw_free(p); nullify(data)
end subroutine

subroutine deallocate_fftw_3d_complex(data)
   complex(CP), dimension(:, :, :), pointer CONTIG :: data
   type(c_ptr) :: p
   p = c_loc(data(1, 1, 1)); call fftw_free(p); nullify(data)
end subroutine

subroutine deallocate_fftw_3d_real(data)
   real(RP), dimension(:, :, :), pointer CONTIG :: data
   type(c_ptr) :: p
   p = c_loc(data(1, 1, 1)); call fftw_free(p); nullify(data)
end subroutine

subroutine poisfft_plan1d_finalize(plan)
   type(PoisFFT_Plan1D) :: plan

   if (c_associated(plan % planptr) .and. plan % planowner) call fftw_destroy_plan(plan % planptr)
   plan % planptr = c_null_ptr
end subroutine

subroutine poisfft_plan2d_finalize(plan)
   type(PoisFFT_Plan2D) :: plan

   if (c_associated(plan % planptr) .and. plan % planowner) then
#ifdef MPI
      if (plan % distributed .and. plan % method == FFT_DISTRIBUTED_PFFT) then
         call pfft_destroy_plan(plan % planptr)
      else
         call fftw_destroy_plan(plan % planptr)
      end if
#else
      call fftw_destroy_plan(plan % planptr)
#endif
   end if
   plan % planptr = c_null_ptr
end subroutine

subroutine poisfft_plan3d_finalize(plan)
   type(PoisFFT_Plan3D) :: plan

   if (c_associated(plan % planptr) .and. plan % planowner) then
#ifdef MPI
      if (plan % distributed) then
         call pfft_destroy_plan(plan % planptr)
      else
         call fftw_destroy_plan(plan % planptr)
      end if
#else
      call fftw_destroy_plan(plan % planptr)
#endif
   end if
   plan % planptr = c_null_ptr
end subroutine

subroutine poisfft_initthreads(nthreads) !instructs fftw to plan to use nthreads threads
   integer, intent(in) :: nthreads
#if defined(_OPENMP) || defined(ENABLE_PTHREADS)
   integer(c_int) :: error
   error = fftw_init_threads()

   if (error == 0) then
      write(*, *) 'Error when initializing FFTW for threads.'
   else
      call fftw_plan_with_nthreads(int(nthreads, c_int))
   end if
#endif
end subroutine

#ifdef MPI
subroutine poisfft_pfft_init
   call pfft_init
end subroutine

subroutine poisfft_pfft_initthreads(nthreads) !instructs PFFT and fftw to plan to use nthreads threads
   integer, intent(in) :: nthreads
#if defined(_OPENMP) || defined(ENABLE_PTHREADS)
   call pfft_plan_with_nthreads(int(nthreads, c_int))
#endif
end subroutine

#endif

#undef RP
#undef CP
#undef FFTW_EXECUTE_COMPLEX
#undef FFTW_EXECUTE_REAL
#undef PFFT_EXECUTE_GEN
#undef FFTW_MPI_EXECUTE_COMPLEX
