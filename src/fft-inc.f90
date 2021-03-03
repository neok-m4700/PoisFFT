#ifndef NO_CONTIGUOUS
#define CONTIG ,contiguous
#else
#define CONTIG
#endif

#if (PREC==2)
#define RP drp
#define CP dcp
#define FFTW_EXECUTE_COMPLEX fftw_execute_dft
#define FFTW_MPI_EXECUTE_COMPLEX fftw_mpi_execute_dft
#define FFTW_EXECUTE_REAL fftw_execute_r2r
#define PFFT_EXECUTE_GEN pfft_execute
#else
#define RP srp
#define CP scp
#define FFTW_EXECUTE_COMPLEX fftwf_execute_dft
#define FFTW_MPI_EXECUTE_COMPLEX fftwf_mpi_execute_dft
#define FFTW_EXECUTE_REAL fftwf_execute_r2r
#define PFFT_EXECUTE_GEN pfftf_execute
#endif
use iso_c_binding
use poisfft_constants
use fftw3
#ifdef MPI
  use pfft
#endif
  implicit none



  integer, parameter :: FFT_Complex = 0
  integer, parameter :: FFT_RealEven00  = FFTW_REDFT00
  integer, parameter :: FFT_RealEven01  = FFTW_REDFT01
  integer, parameter :: FFT_RealEven10  = FFTW_REDFT10
  integer, parameter :: FFT_RealEven11  = FFTW_REDFT11
  integer, parameter :: FFT_RealOdd00   = FFTW_RODFT00
  integer, parameter :: FFT_RealOdd01   = FFTW_RODFT01
  integer, parameter :: FFT_RealOdd10   = FFTW_RODFT10
  integer, parameter :: FFT_RealOdd11   = FFTW_RODFT11

  integer, parameter :: FFT_DISTRIBUTED_FFTW = 1
  integer, parameter :: FFT_DISTRIBUTED_PFFT = 2

  type PoisFFT_Plan1D
    type(c_ptr)         :: planptr = c_null_ptr
    logical             :: planowner = .false.
    logical             :: distributed = .false.
    integer(c_int)      :: dir
  end type PoisFFT_Plan1D


  type PoisFFT_Plan2D
    type(c_ptr)         :: planptr = c_null_ptr
    logical             :: planowner = .false.
    logical             :: distributed = .false.
    integer             :: method = FFT_DISTRIBUTED_PFFT
    integer(c_int)      :: dir
  end type PoisFFT_Plan2D


  type PoisFFT_Plan3D
    type(c_ptr)         :: planptr = c_null_ptr
    logical             :: planowner = .false.
    logical             :: distributed = .false.
    integer(c_int)      :: dir
  end type PoisFFT_Plan3D


  type mpi_vars_1d
    !MPI communicator for the exchange
    integer :: comm = -1
    !our rank in comm
    integer :: rank
    !number of processes in comm
    integer :: np
    ! s for dimensions of send buffers, r for receive buffers
    ! nx, nz are the dimensions of individual parts which will be sent to or received from other processes in comm
    ! displs are the displacements (offsets) of individual pieces in the whole 1D buffer
    ! counts are the the number of elements in each piece
    ! sumrnzs(i) is the sum of rnzs in pieces 1..i-1
    integer,dimension(:),allocatable :: snxs, snzs, &
                                        rnxs, rnzs, &
                                        sdispls, scounts, &
                                        rdispls, rcounts, &
                                        sumrnzs
    !contiguous 3D buffer for locally transposed RHS
    real(RP), allocatable :: tmp1(:,:,:)
    !1D buffer for globally transposed blocks of tmp1 from different processes after MPI_Alltoallv
    real(RP), allocatable :: tmp2(:)
    !3D buffer for contiguous 1D rows on which 1D FFT can be performed locally,
    ! constructed by local reordering of tmp2
    real(RP), allocatable :: rwork(:,:,:)
    ! global indexes corresponding to the given y or z line in the rwork array.
    integer, allocatable :: glob_i(:,:) ! value of global i, (:,:) local indexes of rwork
    
  end type

  type PoisFFT_Solver1D
    real(RP) :: Lx
    integer(c_int) :: nxyz(1)
    integer(c_int) :: nx
    integer(c_int) :: gnx
    integer(c_int) :: offx = 0 !offset from global index
    integer(c_size_t) :: cnt
    integer(c_size_t) :: gcnt
    real(RP) :: norm_factor
    integer(c_int), dimension(2) :: BCs
    real(RP), allocatable, dimension(:) :: denomx
    integer :: approximation = 0
    
    type(PoisFFT_Plan1D) :: forward, backward
    complex(CP), dimension(:), &
#ifndef NO_CONTIGUOUS
    contiguous, &
#endif
                       pointer :: cwork => null()
    real(RP), dimension(:), &
#ifndef NO_CONTIGUOUS
    contiguous, &
#endif
                       pointer :: rwork => null()
    logical :: mpi_transpose_needed = .false.
    integer :: nthreads = 1
    type(mpi_vars_1D) :: mpi
  end type PoisFFT_Solver1D

  type mpi_vars_2d
    integer :: rank
    integer :: np
    integer :: comm = -1
    integer :: comm_dim = 2
  end type

  type PoisFFT_Solver2D
    real(RP) :: Lx, Ly
    integer(c_int) :: nxyz(2)
    integer(c_int) :: nx, ny
    integer(c_int) :: gnx, gny
    integer(c_int) :: offx = 0, offy = 0 !offsets from global indexes
    integer(c_size_t) :: cnt
    integer(c_size_t) :: gcnt
    real(RP) :: norm_factor
    integer(c_int), dimension(4) :: BCs
    real(RP), allocatable, dimension(:) :: denomx, denomy
    integer :: approximation = 0
    
    type(PoisFFT_Plan2D) :: forward, backward
    complex(CP), dimension(:,:), &
#ifndef NO_CONTIGUOUS
    contiguous, &
#endif
                       pointer :: cwork => null()
    real(RP), dimension(:,:), &
#ifndef NO_CONTIGUOUS
    contiguous, &
#endif
                       pointer :: rwork => null()
    integer :: nthreads = 1
    type(mpi_vars_2D) :: mpi
  end type PoisFFT_Solver2D
  
  type mpi_vars_3D
    integer :: comm = -1
  end type

  type PoisFFT_Solver3D
    real(RP) :: Lx, Ly, Lz
    integer(c_int) :: nxyz(3)
    integer(c_int) :: nx, ny, nz
    integer(c_int) :: gnx, gny, gnz
    integer(c_int) :: offx = 0, offy = 0, offz = 0 !offsets from global indexes
    integer(c_size_t) :: cnt
    integer(c_size_t) :: gcnt
    real(RP) :: norm_factor
    integer(c_int), dimension(6) :: BCs
    real(RP), allocatable, dimension(:) :: denomx, denomy, denomz
    integer :: approximation = 0
    
    type(PoisFFT_Plan3D) :: forward, backward
    complex(CP), dimension(:,:,:), &
#ifndef NO_CONTIGUOUS
    contiguous, &
#endif
                       pointer :: cwork => null()
    real(RP), dimension(:,:,:), &
#ifndef NO_CONTIGUOUS
    contiguous, &
#endif
                       pointer :: rwork => null()
    integer :: nthreads = 1
    !will be used in splitting for some boundary conditions
    type(PoisFFT_Solver1D),dimension(:),allocatable :: Solvers1D
    type(PoisFFT_Solver2D),dimension(:),allocatable :: Solvers2D
    type(mpi_vars_3D) :: mpi
  end type PoisFFT_Solver3D
 
  type, extends(PoisFFT_Solver3D) :: PoisFFT_Solver3D_nonuniform_Z
    real(RP) :: z_start, z_end
    real(RP), allocatable :: z(:), z_u(:)
    real(RP), allocatable :: mat_a(:), mat_b(:), mat_c(:) !tridiagonal matrix elements, variable in z
  end type
 
 

  interface deallocate_fftw
    module procedure deallocate_fftw_1D_complex
    module procedure deallocate_fftw_1D_real
    module procedure deallocate_fftw_2D_complex
    module procedure deallocate_fftw_2D_real
    module procedure deallocate_fftw_3D_complex
    module procedure deallocate_fftw_3D_real
  end interface

  interface allocate_fftw_real
    module procedure allocate_fftw_1D_real
    module procedure allocate_fftw_2D_real
    module procedure allocate_fftw_3D_real
  end interface

  interface allocate_fftw_complex
    module procedure allocate_fftw_1D_complex
    module procedure allocate_fftw_2D_complex
    module procedure allocate_fftw_3D_complex
  end interface

  interface Execute
    module procedure PoisFFT_Plan1D_Execute_Real
    module procedure PoisFFT_Plan1D_Execute_Complex
    module procedure PoisFFT_Plan2D_Execute_Real
    module procedure PoisFFT_Plan2D_Execute_Complex
    module procedure PoisFFT_Plan3D_Execute_Real
    module procedure PoisFFT_Plan3D_Execute_Complex
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
#define DIM 1
#include "plan_new-inc.f90"
#undef DIM
end function

function poisfft_plan2d_new(self, plantypes, distributed) result(plan)
#define DIM 2
#include "plan_new-inc.f90"
#undef DIM
end function

function poisfft_plan3d_new(self, plantypes, distributed) result(plan)
#define DIM 3
#include "plan_new-inc.f90"
#undef DIM
end function

subroutine allocate_fftw_1D_complex(self)
#define DIM 1
#define ISREAL 0
#include "allocate_fftw-inc.f90"
#undef DIM
#undef ISREAL
end subroutine

subroutine allocate_fftw_1d_real(self)
#define DIM 1
#define ISREAL 1
#include "allocate_fftw-inc.f90"
#undef DIM
#undef ISREAL
end subroutine

subroutine allocate_fftw_2d_complex(self)
#define DIM 2
#define ISREAL 0
#include "allocate_fftw-inc.f90"
#undef DIM
#undef ISREAL
end subroutine

subroutine allocate_fftw_2d_real(self)
#define DIM 2
#define ISREAL 1
#include "allocate_fftw-inc.f90"
#undef DIM
#undef ISREAL
end subroutine

subroutine allocate_fftw_3d_complex(self)
#define DIM 3
#define ISREAL 0
#include "allocate_fftw-inc.f90"
#undef DIM
#undef ISREAL
end subroutine

subroutine allocate_fftw_3d_real(self)
#define DIM 3
#define ISREAL 1
#include "allocate_fftw-inc.f90"
#undef DIM
#undef ISREAL
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
   ! write(ferr, *) 'shape(data)=', shape(data), ' c_associated(planptr)=', c_associated(plan % planptr)
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
         stop 'error, missing `data` argument in a call to execute_mpi with fftw.'
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
   type(poisfft_plan1d) :: plan

   if (c_associated(plan % planptr) .and. plan % planowner) call fftw_destroy_plan(plan % planptr)
   plan % planptr = c_null_ptr
end subroutine

subroutine poisfft_plan2d_finalize(plan)
   type(poisfft_plan2d) :: plan

   if (c_associated(plan % planptr) .and. plan % planowner) then
#ifdef MPI
      if (plan % distributed .and. plan % method == fft_distributed_pfft) then
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
   type(poisfft_plan3d) :: plan

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

subroutine poisfft_initthreads(nthreads) ! instructs fftw to plan to use nthreads threads
   integer, intent(in) :: nthreads
#if defined(_OPENMP) || defined(ENABLE_PTHREADS)
   integer(c_int) :: error
   error = fftw_init_threads()

   if (error == 0) then
      write(ferr, *) 'error when initializing fftw for threads.'
   else
      call fftw_plan_with_nthreads(int(nthreads, c_int))
   end if
#endif
end subroutine

#ifdef MPI
subroutine poisfft_pfft_init
   call pfft_init
end subroutine

subroutine poisfft_pfft_initthreads(nthreads) ! instructs PFFT and fftw to plan to use nthreads threads
   integer, intent(in) :: nthreads
#if defined(_OPENMP) || defined(ENABLE_PTHREADS)
   call pfft_plan_with_nthreads(int(nthreads, c_int))
#endif
end subroutine

#endif

#undef RP
#undef CP
#undef fftw_execute_complex
#undef fftw_execute_real
#undef pfft_execute_gen
#undef fftw_mpi_execute_complex
