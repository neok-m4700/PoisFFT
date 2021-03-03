#if (PREC==2)

#define RP drp
#define CP dcp
#define MPI_RP mpi_double_precision

#else

#define RP srp
#define CP scp
#define MPI_RP mpi_real

#endif

use iso_c_binding
!$ use omp_lib
use poisfft_constants
#ifdef MPI
  use mpi
#endif
  implicit none

  private
  public :: PoisFFT_Solver1D, PoisFFT_Solver2D, PoisFFT_Solver3D, &
            PoisFFT_Solver3D_nonuniform_z, &
            Finalize, Execute

  real(RP), parameter, private :: pi = 3.141592653589793238462_RP


  interface PoisFFT_Solver1D
    module procedure PoisFFT_Solver1D__New
  end interface

  interface PoisFFT_Solver2D
    module procedure PoisFFT_Solver2D__New
  end interface

  interface PoisFFT_Solver3D
    module procedure PoisFFT_Solver3D__New
  end interface

  interface PoisFFT_Solver3D_nonuniform_z
    module procedure PoisFFT_Solver3D_nonuniform_z__New
  end interface

  interface Finalize
    module procedure PoisFFT_Solver1D__Finalize
    module procedure PoisFFT_Solver2D__Finalize
    module procedure PoisFFT_Solver3D__Finalize
    module procedure PoisFFT_Solver3D_nonuniform_z__Finalize
  end interface Finalize


  interface Init
    module procedure PoisFFT_Solver1D_Init
    module procedure PoisFFT_Solver2D_Init
    module procedure PoisFFT_Solver3D_Init
    module procedure PoisFFT_Solver3D_nonuniform_z_Init
  end interface Init

  interface Execute
    module procedure PoisFFT_Solver1D__Execute
    module procedure PoisFFT_Solver2D__Execute
    module procedure PoisFFT_Solver3D__Execute
    module procedure PoisFFT_Solver3D_nonuniform_z__Execute
  end interface Execute


  contains


    function PoisFFT_Solver3D__New(nxyz,Lxyz,BCs,approximation, &
                                  gnxyz,offs,mpi_comm,nthreads) result(D)
      type(PoisFFT_Solver3D) :: D

      integer, intent(in)   :: nxyz(3)
      real(RP), intent(in)  :: Lxyz(3)
      integer, intent(in)   :: bcs(6)
      integer, intent(in), optional :: approximation
      integer, intent(in), optional :: gnxyz(3)
      integer, intent(in), optional :: offs(3)
      integer, intent(in), optional :: mpi_comm
      integer, intent(in), optional :: nthreads

      D%nxyz = nxyz

      D%Lx = Lxyz(1)
      D%Ly = Lxyz(2)
      D%Lz = Lxyz(3)

      D%nx = nxyz(1)
      D%ny = nxyz(2)
      D%nz = nxyz(3)
      
      if (present(gnxyz)) then
        D%gnx = gnxyz(1)
        D%gny = gnxyz(2)
        D%gnz = gnxyz(3)
      else
        D%gnx = D%nx
        D%gny = D%ny
        D%gnz = D%nz
      end if

      if (present(offs)) then
        D%offx = offs(1)
        D%offy = offs(2)
        D%offz = offs(3)
      end if

      D%cnt = product(D%nxyz)
      
      D%gcnt = product(int([D%gnx, D%gny, D%gnz], kind(D%gcnt)))

      D%BCs = BCs
      
      if (present(approximation)) D%approximation = approximation
      
      if (present(mpi_comm)) then
        D%mpi%comm = mpi_comm
#ifdef MPI      
      else
        stop "No PFFT comm present in PoisFFT_Solver3D__New."
#endif
      end if

      if (present(nthreads)) then
        D%nthreads = nthreads
      else
        D%nthreads = 1
      endif

      !create fftw plans and allocate working arrays
      call Init(D)
    end function PoisFFT_Solver3D__New
    
    


    subroutine PoisFFT_Solver3D_Init(D)
      type(PoisFFT_Solver3D), intent(inout) :: D
      integer :: real_forw, real_back
      integer :: real_forw_xy, real_back_xy
      integer :: real_forw_z, real_back_z
      integer :: i
      !$omp parallel
      !$omp single
      !$ D%nthreads = omp_get_num_threads()
      !$omp end single
      !$omp end parallel

      D%norm_factor = norm_factor(D%gnx,D%BCs(1:2)) * &
                      norm_factor(D%gny,D%BCs(3:4)) * &
                      norm_factor(D%gnz,D%BCs(5:6))
      
      allocate(D%denomx(D%nx))
      allocate(D%denomy(D%ny))
      allocate(D%denomz(D%nz))
      
      
#ifdef MPI
   else
      stop 'no pfft comm present in poisfft_solver3d__new.'
#endif
   end if

   self % nthreads = merge(nthreads, 1, present(nthreads))
   call init(self) ! create fftw plans and allocate working array
end function

subroutine poisfft_solver3d_init(self)
   type(poisfft_solver3d), intent(inout) :: self
   integer :: real_forw, real_back, i
   !$omp parallel
   !$omp single
   !$ self%nthreads = omp_get_num_threads()
   !$omp end single
   !$omp end parallel

   self % norm_factor = &
      norm_factor(self % gnx, self % BCs(1:2)) * &
      norm_factor(self % gny, self % BCs(3:4)) * &
      norm_factor(self % gnz, self % BCs(5:6))

   allocate(self % denomx(self % nx), self % denomy(self % ny), self % denomz(self % nz))

#ifdef MPI
   call poisfft_pfft_init
#else
   if (self % nthreads > 1) call poisfft_initthreads(1)
#endif

        call allocate_fftw_real(D)
        
        real_forw = real_transform_type_forward(D%BCs(1:2))
        real_back = real_transform_type_backward(D%BCs(1:2))

        D%forward = PoisFFT_Plan3D(D, [(real_forw, i=1,3)])
        D%backward = PoisFFT_Plan3D(D, [(real_back, i=1,3)])

      else if (all(D%BCs(1:4)==PoisFFT_Periodic) .and. &
               (all(D%BCs(5:6)==PoisFFT_Neumann) .or. &
                all(D%BCs(5:6)==PoisFFT_NeumannStag) )) then

#ifdef MPI
        allocate(D%Solvers1D(3:2+D%nthreads))

        D%Solvers1D(3) = PoisFFT_Solver1D_From3D(D,3)

        call allocate_fftw_real(D%Solvers1D(3))

        real_forw = real_transform_type_forward(D%BCs(5:6))
        real_back = real_transform_type_backward(D%BCs(5:6))
        D%Solvers1D(3)%forward = PoisFFT_Plan1D(D%Solvers1D(3), [real_forw, real_forw])
        D%Solvers1D(3)%backward = PoisFFT_Plan1D(D%Solvers1D(3), [real_back, real_back])

        do i = 4, 2+D%nthreads
          D%Solvers1D(i) = D%Solvers1D(3)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_real(D%Solvers1D(i))
        end do

        if (D%ny<D%gny) then
          if (D%nthreads>1) call PoisFFT_PFFT_InitThreads(D%nthreads)

          allocate(D%Solvers2D(1))

          D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

          call allocate_fftw_complex(D%Solvers2D(1))

          D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_FORWARD])
          D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_BACKWARD])

          
        else
          allocate(D%Solvers2D(D%nthreads))

          D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

          call allocate_fftw_complex(D%Solvers2D(1))

          D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), &
                                                  [FFT_Complex, FFTW_FORWARD], &
                                                  distributed=.false.)
          D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), &
                                                   [FFT_Complex, FFTW_BACKWARD], &
                                                   distributed=.false.)

          do i = 2, D%nthreads
            D%Solvers2D(i) = D%Solvers2D(1)
            D%Solvers2D(i)%forward%planowner = .false.
            D%Solvers2D(i)%backward%planowner = .false.

            call allocate_fftw_complex(D%Solvers2D(i))
          end do
        end if

        call Init_MPI_Buffers(D, 3)

#else

        allocate(D%Solvers1D(D%nthreads))

        D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,3)

        call allocate_fftw_real(D%Solvers1D(1))

        real_forw = real_transform_type_forward(D%BCs(5:6))
        real_back = real_transform_type_backward(D%BCs(5:6))
        D%Solvers1D(1)%forward = PoisFFT_Plan1D(D%Solvers1D(1), [real_forw, real_forw])
        D%Solvers1D(1)%backward = PoisFFT_Plan1D(D%Solvers1D(1), [real_back, real_back])

        do i = 2, D%nthreads
          D%Solvers1D(i) = D%Solvers1D(1)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_real(D%Solvers1D(i))
        end do

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

        call allocate_fftw_complex(D%Solvers2D(1))

        D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_FORWARD])
        D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_BACKWARD])

        do i = 2, D%nthreads
          D%Solvers2D(i) = D%Solvers2D(1)
          D%Solvers2D(i)%forward%planowner = .false.
          D%Solvers2D(i)%backward%planowner = .false.

          call allocate_fftw_complex(D%Solvers2D(i))
        end do
        
!MPI
#endif

      else if (all(D%BCs(3:6)==PoisFFT_Periodic) .and. &
               (all(D%BCs(1:2)==PoisFFT_Neumann) .or. &
                all(D%BCs(1:2)==PoisFFT_NeumannStag))) then
                
        allocate(D%Solvers1D(D%nthreads))

        D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,1)

        call allocate_fftw_real(D%Solvers1D(1))

        real_forw = real_transform_type_forward(D%BCs(1:2))
        real_back = real_transform_type_backward(D%BCs(1:2))
        D%Solvers1D(1)%forward = PoisFFT_Plan1D(D%Solvers1D(1), [real_forw, real_forw])
        D%Solvers1D(1)%backward = PoisFFT_Plan1D(D%Solvers1D(1), [real_back, real_back])

        do i = 2, D%nthreads
          D%Solvers1D(i) = D%Solvers1D(1)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_real(D%Solvers1D(i))
        end do

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,1)

        call allocate_fftw_complex(D%Solvers2D(1))

        D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_FORWARD])
        D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_BACKWARD])

        do i = 2, D%nthreads
          D%Solvers2D(i) = D%Solvers2D(1)
          D%Solvers2D(i)%forward%planowner = .false.
          D%Solvers2D(i)%backward%planowner = .false.

          call allocate_fftw_complex(D%Solvers2D(i))
        end do                
                
      else if (all(D%BCs(1:2)==PoisFFT_Periodic) .and. &
               all(D%BCs(5:6)==PoisFFT_Periodic) .and. &
               (all(D%BCs(3:4)==PoisFFT_Neumann) .or. &
                all(D%BCs(3:4)==PoisFFT_NeumannStag) )) then

        allocate(D%Solvers1D(D%nthreads))

        D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,2)

      self % forward = poisfft_plan3d(self, [fft_complex, fftw_forward])
      self % backward = poisfft_plan3d(self, [fft_complex, fftw_backward])

        real_forw = real_transform_type_forward(D%BCs(3:4))
        real_back = real_transform_type_backward(D%BCs(3:4))
        D%Solvers1D(1)%forward = PoisFFT_Plan1D(D%Solvers1D(1), [real_forw, real_forw])
        D%Solvers1D(1)%backward = PoisFFT_Plan1D(D%Solvers1D(1), [real_back, real_back])

        do i = 2, D%nthreads
          D%Solvers1D(i) = D%Solvers1D(1)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_real(D%Solvers1D(i))
        end do

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,2)

        call allocate_fftw_complex(D%Solvers2D(1))

        D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_FORWARD])
        D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_BACKWARD])

        do i = 2, D%nthreads
          D%Solvers2D(i) = D%Solvers2D(1)
          D%Solvers2D(i)%forward%planowner = .false.
          D%Solvers2D(i)%backward%planowner = .false.

          call allocate_fftw_complex(D%Solvers2D(i))
        end do
        
      else if (all(D%BCs(1:2)==PoisFFT_Periodic) .and. &
               (all(D%BCs(3:6)==PoisFFT_Neumann) .or. &
                all(D%BCs(3:6)==PoisFFT_NeumannStag) )) then

        real_forw = real_transform_type_forward(D%BCs(3:4))
        real_back = real_transform_type_backward(D%BCs(3:4))
        
#ifdef MPI
      if (self % nthreads > 1) call poisfft_pfft_initthreads(self % nthreads)
#else
      if (self % nthreads > 1) call poisfft_initthreads(self % nthreads)
#endif




      else if ((all(D%BCs(1:4)==PoisFFT_Dirichlet) .or. &
                all(D%BCs(1:4)==PoisFFT_DirichletStag)) .and. &
               all(D%BCs(5:6)==PoisFFT_Periodic)) then

        real_forw = real_transform_type_forward(D%BCs(1:2))
        real_back = real_transform_type_backward(D%BCs(1:2))
        
! #ifdef MPI
!         !1..x, 2..y, 3..z, not different threads, but directions
!         allocate(D%Solvers1D(3*D%nthreads))
! 
!         D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,1)
!         D%Solvers1D(2) = PoisFFT_Solver1D_From3D(D,2)
!         D%Solvers1D(3) = PoisFFT_Solver1D_From3D(D,3)
!         
!         call allocate_fftw_complex(D%Solvers1D(1))
!         call allocate_fftw_real(D%Solvers1D(2))
!         call allocate_fftw_real(D%Solvers1D(3))
!         
!         D%Solvers1D(1)%forward = PoisFFT_Plan1D(D%Solvers1D(1), [FFT_Complex, FFTW_FORWARD])
!         D%Solvers1D(1)%backward = PoisFFT_Plan1D(D%Solvers1D(1), [FFT_Complex, FFTW_BACKWARD])
!         D%Solvers1D(2)%forward = PoisFFT_Plan1D(D%Solvers1D(2), [real_forw])
!         D%Solvers1D(2)%backward = PoisFFT_Plan1D(D%Solvers1D(2), [real_back])
!         D%Solvers1D(3)%forward = PoisFFT_Plan1D(D%Solvers1D(3), [real_forw])
!         D%Solvers1D(3)%backward = PoisFFT_Plan1D(D%Solvers1D(3), [real_back])
!         
!         do i = 4, 1 + 3*(D%nthreads-1), 3
!           D%Solvers1D(i) = D%Solvers1D(1)
!           D%Solvers1D(i)%forward%planowner = .false.
!           D%Solvers1D(i)%backward%planowner = .false.
!           call allocate_fftw_complex(D%Solvers1D(i))
!         end do
! 
!         do i = 5, 2 + 3*(D%nthreads-1), 3
!           D%Solvers1D(i) = D%Solvers1D(2)
!           D%Solvers1D(i)%forward%planowner = .false.
!           D%Solvers1D(i)%backward%planowner = .false.
!           call allocate_fftw_real(D%Solvers1D(i))
!         end do
!         do i = 6, 3 + 3*(D%nthreads-1), 3
!           D%Solvers1D(i) = D%Solvers1D(3)
!           D%Solvers1D(i)%forward%planowner = .false.
!           D%Solvers1D(i)%backward%planowner = .false.
!           call allocate_fftw_real(D%Solvers1D(i))
!         end do
! 
!         call Init_MPI_Buffers(D, 2)
!         call Init_MPI_Buffers(D, 3)
! 
! #else
        allocate(D%Solvers1D(D%nthreads))

        D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,3)

        call allocate_fftw_complex(D%Solvers1D(1))

        D%Solvers1D(1)%forward = PoisFFT_Plan1D(D%Solvers1D(1), [FFT_Complex, FFTW_FORWARD])
        D%Solvers1D(1)%backward = PoisFFT_Plan1D(D%Solvers1D(1), [FFT_Complex, FFTW_BACKWARD])

        do i=2,D%nthreads
          D%Solvers1D(i) = D%Solvers1D(1)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_complex(D%Solvers1D(i))
        end do

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

        call allocate_fftw_real(D%Solvers2D(1))

        D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [real_forw, real_forw])
        D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [real_back, real_back])

        do i=2,D%nthreads
          D%Solvers2D(i) = D%Solvers2D(1)
          D%Solvers2D(i)%forward%planowner = .false.
          D%Solvers2D(i)%backward%planowner = .false.

          call allocate_fftw_real(D%Solvers2D(i))

        end do
! #endif



      else if ( (all(D%BCs==PoisFFT_Neumann .or. &
                     D%BCs==PoisFFT_NeumannStag .or. &
                     D%BCs==PoisFFT_Dirichlet .or. &
                     D%BCs==PoisFFT_DirichletStag)) &
               .and. & 
               (all(D%BCs(1:2)==D%BCs(3:4)) .and. &
                any(D%BCs(1:2)/=D%BCs(5:6))) &
               .and. &
               (any(D%BCs==PoisFFT_Dirichlet .or. &
                     D%BCs==PoisFFT_DirichletStag)) &
              ) then

        real_forw_xy = real_transform_type_forward(D%BCs(1:2))
        real_back_xy = real_transform_type_backward(D%BCs(1:2))
        real_forw_z = real_transform_type_forward(D%BCs(5:6))
        real_back_z = real_transform_type_backward(D%BCs(5:6))
        
! #ifdef MPI
!         !1..x, 2..y, 3..z, not different threads, but directions
!         allocate(D%Solvers1D(3*D%nthreads))
! 
!         D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,1)
!         D%Solvers1D(2) = PoisFFT_Solver1D_From3D(D,2)
!         D%Solvers1D(3) = PoisFFT_Solver1D_From3D(D,3)
!         
!         call allocate_fftw_complex(D%Solvers1D(1))
!         call allocate_fftw_real(D%Solvers1D(2))
!         call allocate_fftw_real(D%Solvers1D(3))
!         
!         D%Solvers1D(1)%forward = PoisFFT_Plan1D(D%Solvers1D(1), [FFT_Complex, FFTW_FORWARD])
!         D%Solvers1D(1)%backward = PoisFFT_Plan1D(D%Solvers1D(1), [FFT_Complex, FFTW_BACKWARD])
!         D%Solvers1D(2)%forward = PoisFFT_Plan1D(D%Solvers1D(2), [real_forw])
!         D%Solvers1D(2)%backward = PoisFFT_Plan1D(D%Solvers1D(2), [real_back])
!         D%Solvers1D(3)%forward = PoisFFT_Plan1D(D%Solvers1D(3), [real_forw])
!         D%Solvers1D(3)%backward = PoisFFT_Plan1D(D%Solvers1D(3), [real_back])
!         
!         do i = 4, 1 + 3*(D%nthreads-1), 3
!           D%Solvers1D(i) = D%Solvers1D(1)
!           D%Solvers1D(i)%forward%planowner = .false.
!           D%Solvers1D(i)%backward%planowner = .false.
!           call allocate_fftw_complex(D%Solvers1D(i))
!         end do
! 
!         do i = 5, 2 + 3*(D%nthreads-1), 3
!           D%Solvers1D(i) = D%Solvers1D(2)
!           D%Solvers1D(i)%forward%planowner = .false.
!           D%Solvers1D(i)%backward%planowner = .false.
!           call allocate_fftw_real(D%Solvers1D(i))
!         end do
!         do i = 6, 3 + 3*(D%nthreads-1), 3
!           D%Solvers1D(i) = D%Solvers1D(3)
!           D%Solvers1D(i)%forward%planowner = .false.
!           D%Solvers1D(i)%backward%planowner = .false.
!           call allocate_fftw_real(D%Solvers1D(i))
!         end do
! 
!         call Init_MPI_Buffers(D, 2)
!         call Init_MPI_Buffers(D, 3)
! 
! #else
        allocate(D%Solvers1D(D%nthreads))

        D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,3)

        call allocate_fftw_real(D%Solvers1D(1))

        D%Solvers1D(1)%forward = PoisFFT_Plan1D(D%Solvers1D(1), [real_forw_z, real_forw_z])
        D%Solvers1D(1)%backward = PoisFFT_Plan1D(D%Solvers1D(1), [real_back_z, real_back_z])

        do i=2,D%nthreads
          D%Solvers1D(i) = D%Solvers1D(1)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_real(D%Solvers1D(i))
        end do

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

        call allocate_fftw_real(D%Solvers2D(1))

        D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [real_forw_xy, real_forw_xy])
        D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [real_back_xy, real_back_xy])

        do i=2,D%nthreads
          D%Solvers2D(i) = D%Solvers2D(1)
          D%Solvers2D(i)%forward%planowner = .false.
          D%Solvers2D(i)%backward%planowner = .false.

          call allocate_fftw_real(D%Solvers2D(i))

        end do
! #endif



      else
        stop "Unknown combination of boundary conditions."
      endif


      if (D%approximation==PoisFFT_FiniteDifference2) then
        D%denomx = eigenvalues(eig_fn_FD2, D%BCs(1:2), D%Lx, D%nx, D%gnx, D%offx)
        D%denomy = eigenvalues(eig_fn_FD2, D%BCs(3:4), D%Ly, D%ny, D%gny, D%offy)
        D%denomz = eigenvalues(eig_fn_FD2, D%BCs(5:6), D%Lz, D%nz, D%gnz, D%offz)
      else if (D%approximation==PoisFFT_FiniteDifference4) then
        D%denomx = eigenvalues(eig_fn_FD4, D%BCs(1:2), D%Lx, D%gnx, D%gnx, D%offx)
        D%denomy = eigenvalues(eig_fn_FD4, D%BCs(3:4), D%Ly, D%gny, D%gny, D%offy)
        D%denomz = eigenvalues(eig_fn_FD4, D%BCs(5:6), D%Lz, D%gnz, D%gnz, D%offz)
      else
        D%denomx = eigenvalues(eig_fn_spectral, D%BCs(1:2), D%Lx, D%nx, D%gnx, D%offx)
        D%denomy = eigenvalues(eig_fn_spectral, D%BCs(3:4), D%Ly, D%ny, D%gny, D%offy)
        D%denomz = eigenvalues(eig_fn_spectral, D%BCs(5:6), D%Lz, D%nz, D%gnz, D%offz)
      end if
      
    end subroutine PoisFFT_Solver3D_Init
    
    


    subroutine PoisFFT_Solver3D__Execute(D,Phi,RHS)
      type(PoisFFT_Solver3D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)

      integer   :: ngPhi(3), ngRHS(3)

      ngPhi = (ubound(Phi)-[D%nx,D%ny,D%nz])/2
      ngRHS = (ubound(RHS)-[D%nx,D%ny,D%nz])/2

      if (all(D%BCs==PoisFFT_Periodic)) then

        call PoisFFT_Solver3D_FullPeriodic(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))

      else if (all(D%BCs==PoisFFT_Dirichlet) .or. &
               all(D%BCs==PoisFFT_DirichletStag)) then

        call PoisFFT_Solver3D_FullDirichlet(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))

      else if (all(D%BCs==PoisFFT_Neumann) .or. &
               all(D%BCs==PoisFFT_NeumannStag)) then

        call PoisFFT_Solver3D_FullNeumann(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))

      else if (all(D%BCs(1:4)==PoisFFT_Periodic) .and. &
               (all(D%BCs(5:6)==PoisFFT_Neumann) .or. &
                all(D%BCs(5:6)==PoisFFT_NeumannStag))) then
        call PoisFFT_Solver3D_PPNs(D,&
                Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                    ngPhi(2)+1:ngPhi(2)+D%ny,&
                    ngPhi(3)+1:ngPhi(3)+D%nz),&
                RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                    ngRHS(2)+1:ngRHS(2)+D%ny,&
                    ngRHS(3)+1:ngRHS(3)+D%nz))

      else if (all(D%BCs(3:6)==PoisFFT_Periodic) .and. &
               (all(D%BCs(1:2)==PoisFFT_Neumann) .or. &
                all(D%BCs(1:2)==PoisFFT_NeumannStag))) then
        call PoisFFT_Solver3D_NsPP(D,&
                Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                    ngPhi(2)+1:ngPhi(2)+D%ny,&
                    ngPhi(3)+1:ngPhi(3)+D%nz),&
                RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                    ngRHS(2)+1:ngRHS(2)+D%ny,&
                    ngRHS(3)+1:ngRHS(3)+D%nz))

      else if (all(D%BCs(1:2)==PoisFFT_Periodic) .and. &
               (all(D%BCs(3:6)==PoisFFT_Neumann) .or. &
                all(D%BCs(3:6)==PoisFFT_NeumannStag))) then

        call PoisFFT_Solver3D_PNsNs(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))

      else if (all(D%BCs(1:2)==PoisFFT_Periodic) .and. &
               all(D%BCs(5:6)==PoisFFT_Periodic) .and. &
               (all(D%BCs(3:4)==PoisFFT_Neumann) .or. &
                all(D%BCs(3:4)==PoisFFT_NeumannStag) )) then
                
        call PoisFFT_Solver3D_PNsP(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))

      else if ((all(D%BCs(1:4)==PoisFFT_Dirichlet) .or. &
                all(D%BCs(1:4)==PoisFFT_DirichletStag)) .and. &
               all(D%BCs(5:6)==PoisFFT_Periodic)) then

        call PoisFFT_Solver3D_DDP(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))
                     
      else if ( (all(D%BCs==PoisFFT_Neumann .or. &
                     D%BCs==PoisFFT_NeumannStag .or. &
                     D%BCs==PoisFFT_Dirichlet .or. &
                     D%BCs==PoisFFT_DirichletStag)) &
               .and. & 
               (all(D%BCs(1:2)==D%BCs(3:4)) .and. &
                any(D%BCs(1:2)/=D%BCs(5:6))) &
               .and. &
               (any(D%BCs==PoisFFT_Dirichlet .or. &
                     D%BCs==PoisFFT_DirichletStag)) &
              ) then

        call PoisFFT_Solver3D_2real1real(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))
      endif

    end subroutine PoisFFT_Solver3D__Execute
    
    
    
    
    subroutine PoisFFT_Solver3D__Finalize(D)
      type(PoisFFT_Solver3D), intent(inout) :: D
      integer :: i

      call Finalize(D%forward)
      call Finalize(D%backward)

      if (associated(D%rwork)) call deallocate_fftw(D%rwork)
      if (associated(D%cwork)) call deallocate_fftw(D%cwork)

      if (allocated(D%Solvers1D)) then
        do i = lbound(D%Solvers1D,1), ubound(D%Solvers1D,1)
          call Finalize(D%Solvers1D(i))
        end do

        deallocate(D%Solvers1D)
      endif

      if (allocated(D%Solvers2D)) then
        do i = lbound(D%Solvers2D,1), ubound(D%Solvers2D,1)
          call Finalize(D%Solvers2D(i))
        end do

        deallocate(D%Solvers2D)
      endif

    endsubroutine PoisFFT_Solver3D__Finalize


      call allocate_fftw_real(self)

      real_forw = real_transform_type_forward(self % bcs(1))
      real_back = real_transform_type_backward(self % bcs(1))

      self % forward = poisfft_plan3d(self, [(real_forw, i=1, 3)])
      self % backward = poisfft_plan3d(self, [(real_back, i=1, 3)])

   else if ( &
      all(self % bcs(1:4) == poisfft_periodic) .and. &
      (all(self % bcs(5:6) == poisfft_neumann) .or. all(self % bcs(5:6) == poisfft_neumannstag))) then


    function PoisFFT_Solver1D_From3D(D3D,direction) result(D)
      type(PoisFFT_Solver1D)             :: D
      class(PoisFFT_Solver3D), intent(in) :: D3D
      integer, intent(in)                :: direction
#ifdef MPI
      allocate(self % solvers1d(3:2 + self % nthreads))
      self % solvers1d(3) = poisfft_solver1d_from3d(self, 3)
      call allocate_fftw_real(self % solvers1d(3))
      self % solvers1d(3) % forward = poisfft_plan1d(self % solvers1d(3), [fft_realeven10])
      self % solvers1d(3) % backward = poisfft_plan1d(self % solvers1d(3), [fft_realeven01])

      do i = 4, 2 + self % nthreads
         self % solvers1d(i) = self % solvers1d(3)
         self % solvers1d(i) % forward % planowner = .false.
         self % solvers1d(i) % backward % planowner = .false.
         call allocate_fftw_real(self % solvers1d(i))
      end do

      if (self % ny < self % gny) then
         if (self % nthreads > 1) call poisfft_pfft_initthreads(self % nthreads)
         allocate(self % solvers2d(1))
         self % solvers2d(1) = poisfft_solver2d_from3d(self, 3)
         call allocate_fftw_complex(self % solvers2d(1))
         self % solvers2d(1) % forward = poisfft_plan2d(self % solvers2d(1), [fft_complex, fftw_forward])
         self % solvers2d(1) % backward = poisfft_plan2d(self % solvers2d(1), [fft_complex, fftw_backward])
      else
         allocate(self % solvers2d(self % nthreads))
         self % solvers2d(1) = poisfft_solver2d_from3d(self, 3)
         call allocate_fftw_complex(self % solvers2d(1))
         self % solvers2d(1) % forward = poisfft_plan2d(self % solvers2d(1), [fft_complex, fftw_forward], distributed=.false.)
         self % solvers2d(1) % backward = poisfft_plan2d(self % solvers2d(1), [fft_complex, fftw_backward], distributed=.false.)

         do i = 2, self % nthreads
            self % solvers2d(i) = self % solvers2d(1)
            self % solvers2d(i) % forward % planowner = .false.
            self % solvers2d(i) % backward % planowner = .false.
            call allocate_fftw_complex(self % solvers2d(i))
         end do
      end if
      call init_mpi_buffers(self, 3)
#else
      allocate(self % solvers1d(self % nthreads))
      self % solvers1d(1) = poisfft_solver1d_from3d(self, 3)
      call allocate_fftw_real(self % solvers1d(1))
      self % solvers1d(1) % forward = poisfft_plan1d(self % solvers1d(1), [fft_realeven10])
      self % solvers1d(1) % backward = poisfft_plan1d(self % solvers1d(1), [fft_realeven01])

      do i = 2, self % nthreads
         self % solvers1d(i) = self % solvers1d(1)
         self % solvers1d(i) % forward % planowner = .false.
         self % solvers1d(i) % backward % planowner = .false.
         call allocate_fftw_real(self % solvers1d(i))
      end do

      allocate(self % solvers2d(self % nthreads))
      self % solvers2d(1) = poisfft_solver2d_from3d(self, 3)
      call allocate_fftw_complex(self % solvers2d(1))

      self % solvers2d(1) % forward = poisfft_plan2d(self % solvers2d(1), [fft_complex, fftw_forward])
      self % solvers2d(1) % backward = poisfft_plan2d(self % solvers2d(1), [fft_complex, fftw_backward])

      do i = 2, self % nthreads
         self % solvers2d(i) = self % solvers2d(1)
         self % solvers2d(i) % forward % planowner = .false.
         self % solvers2d(i) % backward % planowner = .false.
         call allocate_fftw_complex(self % solvers2d(i))
      end do
#endif
   else if ( &
      all(self % bcs(1:2) == poisfft_periodic) .and. &
      all(self % bcs(5:6) == poisfft_periodic) .and. &
      (all(self % bcs(3:4) == poisfft_neumann) .or. all(self % bcs(3:4) == poisfft_neumannstag))) then

      allocate(self % solvers1d(self % nthreads))
      self % solvers1d(1) = poisfft_solver1d_from3d(self, 2)
      call allocate_fftw_real(self % solvers1d(1))

      self % solvers1d(1) % forward = poisfft_plan1d(self % solvers1d(1), [fft_realeven10])
      self % solvers1d(1) % backward = poisfft_plan1d(self % solvers1d(1), [fft_realeven01])

      do i = 2, self % nthreads
         self % solvers1d(i) = self % solvers1d(1)
         self % solvers1d(i) % forward % planowner = .false.
         self % solvers1d(i) % backward % planowner = .false.
         call allocate_fftw_real(self % solvers1d(i))
      end do

      allocate(self % solvers2d(self % nthreads))
      self % solvers2d(1) = poisfft_solver2d_from3d(self, 2)
      call allocate_fftw_complex(self % solvers2d(1))

      self % solvers2d(1) % forward = poisfft_plan2d(self % solvers2d(1), [fft_complex, fftw_forward])
      self % solvers2d(1) % backward = poisfft_plan2d(self % solvers2d(1), [fft_complex, fftw_backward])

      do i = 2, self % nthreads
         self % solvers2d(i) = self % solvers2d(1)
         self % solvers2d(i) % forward % planowner = .false.
         self % solvers2d(i) % backward % planowner = .false.
         call allocate_fftw_complex(self % solvers2d(i))
      end do

   else if ( &
      all(self % bcs(1:2) == poisfft_periodic) .and. &
      (all(self % bcs(3:6) == poisfft_neumann) .or. all(self % bcs(3:6) == poisfft_neumannstag))) then

      real_forw = real_transform_type_forward(self % bcs(3))
      real_back = real_transform_type_backward(self % bcs(3))

      if (direction==1) then
        D%nx = D3D%nx
        D%gnx = D3D%gnx
        D%offx = D3D%offx
        D%BCs = D3D%BCs(1:2)
      else if (direction==2) then
        D%nx = D3D%ny
        D%gnx = D3D%gny
        D%offx = D3D%offy
        D%BCs = D3D%BCs(3:4)
      else
        D%nx = D3D%nz
        D%gnx = D3D%gnz
        D%offx = D3D%offz
        D%BCs = D3D%BCs(5:6)
      endif
      
      D%nxyz = [D%nx]

      D%cnt = D%nx
      
      D%gcnt = int(D%gnx, kind(D%gcnt))
    end function PoisFFT_Solver1D_From3D



    function PoisFFT_Solver2D_From3D(D3D,direction) result(D)
      type(PoisFFT_Solver2D)             :: D
      class(PoisFFT_Solver3D), intent(in) :: D3D
      integer, intent(in)                :: direction
#ifdef MPI
      ! 1..x, 2..y, 3..z, not different threads, but directions
      allocate(self % solvers1d(3 * self % nthreads))

      self % solvers1d(1) = poisfft_solver1d_from3d(self, 1)
      self % solvers1d(2) = poisfft_solver1d_from3d(self, 2)
      self % solvers1d(3) = poisfft_solver1d_from3d(self, 3)

      call allocate_fftw_complex(self % solvers1d(1))
      call allocate_fftw_real(self % solvers1d(2))
      call allocate_fftw_real(self % solvers1d(3))

      self % solvers1d(1) % forward = poisfft_plan1d(self % solvers1d(1), [fft_complex, fftw_forward])
      self % solvers1d(1) % backward = poisfft_plan1d(self % solvers1d(1), [fft_complex, fftw_backward])
      self % solvers1d(2) % forward = poisfft_plan1d(self % solvers1d(2), [real_forw])
      self % solvers1d(2) % backward = poisfft_plan1d(self % solvers1d(2), [real_back])
      self % solvers1d(3) % forward = poisfft_plan1d(self % solvers1d(3), [real_forw])
      self % solvers1d(3) % backward = poisfft_plan1d(self % solvers1d(3), [real_back])

      do i = 4, 1 + 3 * (self % nthreads - 1), 3
         self % solvers1d(i) = self % solvers1d(1)
         self % solvers1d(i) % forward % planowner = .false.
         self % solvers1d(i) % backward % planowner = .false.
         call allocate_fftw_complex(self % solvers1d(i))
      end do

      do i = 5, 2 + 3 * (self % nthreads - 1), 3
         self % solvers1d(i) = self % solvers1d(2)
         self % solvers1d(i) % forward % planowner = .false.
         self % solvers1d(i) % backward % planowner = .false.
         call allocate_fftw_real(self % solvers1d(i))
      end do
      do i = 6, 3 + 3 * (self % nthreads - 1), 3
         self % solvers1d(i) = self % solvers1d(3)
         self % solvers1d(i) % forward % planowner = .false.
         self % solvers1d(i) % backward % planowner = .false.
         call allocate_fftw_real(self % solvers1d(i))
      end do

      call init_mpi_buffers(self, 2)
      call init_mpi_buffers(self, 3)
#else
      allocate(self % solvers1d(self % nthreads))
      self % solvers1d(1) = poisfft_solver1d_from3d(self, 1)
      call allocate_fftw_complex(self % solvers1d(1))

      self % solvers1d(1) % forward = poisfft_plan1d(self % solvers1d(1), [fft_complex, fftw_forward])
      self % solvers1d(1) % backward = poisfft_plan1d(self % solvers1d(1), [fft_complex, fftw_backward])

      do i = 2, self % nthreads
         self % solvers1d(i) = self % solvers1d(1)
         self % solvers1d(i) % forward % planowner = .false.
         self % solvers1d(i) % backward % planowner = .false.
         call allocate_fftw_complex(self % solvers1d(i))
      end do
      allocate(self % solvers2d(self % nthreads))
      self % solvers2d(1) = poisfft_solver2d_from3d(self, 1)
      call allocate_fftw_real(self % solvers2d(1))

      self % solvers2d(1) % forward = poisfft_plan2d(self % solvers2d(1), [real_forw, real_forw])
      self % solvers2d(1) % backward = poisfft_plan2d(self % solvers2d(1), [real_back, real_back])

      do i = 2, self % nthreads
         self % solvers2d(i) = self % solvers2d(1)
         self % solvers2d(i) % forward % planowner = .false.
         self % solvers2d(i) % backward % planowner = .false.
         call allocate_fftw_real(self % solvers2d(i))
      end do
#endif
   else
      stop 'unknown combination of boundary conditions.'
   endif

   if (self % approximation == poisfft_finitedifference2) then
      self % denomx = eigenvalues(eig_fn_fd2, self % bcs(1:2), self % lx, self % nx, self % gnx, self % offx)
      self % denomy = eigenvalues(eig_fn_fd2, self % bcs(3:4), self % ly, self % ny, self % gny, self % offy)
      self % denomz = eigenvalues(eig_fn_fd2, self % bcs(5:6), self % lz, self % nz, self % gnz, self % offz)
   else if (self % approximation == poisfft_finitedifference4) then
      self % denomx = eigenvalues(eig_fn_fd4, self % bcs(1:2), self % lx, self % gnx, self % gnx, self % offx)
      self % denomy = eigenvalues(eig_fn_fd4, self % bcs(3:4), self % ly, self % gny, self % gny, self % offy)
      self % denomz = eigenvalues(eig_fn_fd4, self % bcs(5:6), self % lz, self % gnz, self % gnz, self % offz)
   else
      self % denomx = eigenvalues(eig_fn_spectral, self % bcs(1:2), self % lx, self % nx, self % gnx, self % offx)
      self % denomy = eigenvalues(eig_fn_spectral, self % bcs(3:4), self % ly, self % ny, self % gny, self % offy)
      self % denomz = eigenvalues(eig_fn_spectral, self % bcs(5:6), self % lz, self % nz, self % gnz, self % offz)
   end if

contains
   function real_transform_type_forward(bc) result(res)
      integer :: res
      integer, intent(in) :: bc
      select case(bc)
       case(poisfft_dirichlet)
         res = fft_realodd00
       case(poisfft_dirichletstag)
         res = fft_realodd10
       case(poisfft_neumann)
         res = fft_realeven00
       case(poisfft_neumannstag)
         res = fft_realeven10
       case default
         res = -1
      end select
   end function

   function real_transform_type_backward(bc) result(res)
      integer :: res
      integer, intent(in) :: bc
      select case(bc)
       case(poisfft_dirichlet)
         res = fft_realodd00
       case(poisfft_dirichletstag)
         res = fft_realodd01
       case(poisfft_neumann)
         res = fft_realeven00
       case(poisfft_neumannstag)
         res = fft_realeven01
       case default
         res = -1
      end select
   end function
end subroutine

subroutine poisfft_solver3d__execute(self, phi, rhs)
   type(poisfft_solver3d), intent(inout) :: self
   real(RP), intent(out) :: phi(:, :, :)
   real(RP), intent(in) :: rhs(:, :, :)

   integer :: ngphi(3), ngrhs(3)

   ngphi = (ubound(phi) - [self % nx, self % ny, self % nz]) / 2
   ngrhs = (ubound(rhs) - [self % nx, self % ny, self % nz]) / 2

   if (all(self % bcs == poisfft_periodic)) then
      call poisfft_solver3d_fullperiodic(self, &
         phi(ngphi(1) + 1:ngphi(1) + self % nx, ngphi(2) + 1:ngphi(2) + self % ny, ngphi(3) + 1:ngphi(3) + self % nz), &
         rhs(ngrhs(1) + 1:ngrhs(1) + self % nx, ngrhs(2) + 1:ngrhs(2) + self % ny, ngrhs(3) + 1:ngrhs(3) + self % nz))

   else if (all(self % bcs == poisfft_dirichlet) .or. all(self % bcs == poisfft_dirichletstag)) then
      call poisfft_solver3d_fulldirichlet(self, &
         phi(ngphi(1) + 1:ngphi(1) + self % nx, ngphi(2) + 1:ngphi(2) + self % ny, ngphi(3) + 1:ngphi(3) + self % nz), &
         rhs(ngrhs(1) + 1:ngrhs(1) + self % nx, ngrhs(2) + 1:ngrhs(2) + self % ny, ngrhs(3) + 1:ngrhs(3) + self % nz))

   else if (all(self % bcs == poisfft_neumann) .or. all(self % bcs == poisfft_neumannstag)) then
      call poisfft_solver3d_fullneumann(self, &
         phi(ngphi(1) + 1:ngphi(1) + self % nx, ngphi(2) + 1:ngphi(2) + self % ny, ngphi(3) + 1:ngphi(3) + self % nz), &
         rhs(ngrhs(1) + 1:ngrhs(1) + self % nx, ngrhs(2) + 1:ngrhs(2) + self % ny, ngrhs(3) + 1:ngrhs(3) + self % nz))

   else if (all(self % bcs(1:4) == poisfft_periodic) .and. (all(self % bcs(5:6) == poisfft_neumann) .or. all(self % bcs(5:6) == poisfft_neumannstag))) then
      call poisfft_solver3d_ppns(self, &
         phi(ngphi(1) + 1:ngphi(1) + self % nx, ngphi(2) + 1:ngphi(2) + self % ny, ngphi(3) + 1:ngphi(3) + self % nz), &
         rhs(ngrhs(1) + 1:ngrhs(1) + self % nx, ngrhs(2) + 1:ngrhs(2) + self % ny, ngrhs(3) + 1:ngrhs(3) + self % nz))

   else if (all(self % bcs(1:2) == poisfft_periodic) .and. (all(self % bcs(3:6) == poisfft_neumann) .or. all(self % bcs(3:6) == poisfft_neumannstag))) then
      call poisfft_solver3d_pnsns(self, &
         phi(ngphi(1) + 1:ngphi(1) + self % nx, ngphi(2) + 1:ngphi(2) + self % ny, ngphi(3) + 1:ngphi(3) + self % nz), &
         rhs(ngrhs(1) + 1:ngrhs(1) + self % nx, ngrhs(2) + 1:ngrhs(2) + self % ny, ngrhs(3) + 1:ngrhs(3) + self % nz))
   endif
end subroutine

subroutine poisfft_solver3d__finalize(self)
   type(poisfft_solver3d), intent(inout) :: self
   integer :: i

   call finalize(self % forward)
   call finalize(self % backward)

   if (associated(self % rwork)) call deallocate_fftw(self % rwork)
   if (associated(self % cwork)) call deallocate_fftw(self % cwork)

   if (allocated(self % solvers1d)) then
      do i = lbound(self % solvers1d, 1), ubound(self % solvers1d, 1)
         call finalize(self % solvers1d(i))
      end do
      deallocate(self % solvers1d)
   endif

   if (allocated(self % solvers2d)) then
      do i = lbound(self % solvers2d, 1), ubound(self % solvers2d, 1)
         call finalize(self % solvers2d(i))
      end do
      deallocate(self % solvers2d)
   endif
end subroutine

function poisfft_solver1d_from3d(d3d, direction) result(self)
   type(poisfft_solver1d) :: self
   type(poisfft_solver3d), intent(in) :: d3d
   integer, intent(in) :: direction
#ifdef MPI
   integer :: ie, dims

   call mpi_cartdim_get(d3d % mpi % comm, dims, ie)
   if (ie /= 0) stop 'error executing mpi_cartdim_get.'

   if (dims == 2) then
      ! we see the dimensions reversed in fortran!
      select case(direction)
       case(1)
         self % mpi % comm = mpi_comm_self
       case(2)
         call mpi_cart_sub(d3d % mpi % comm, [.false., .true.], self % mpi % comm, ie)
         if (ie /= 0) stop 'error executing mpi_cart_sub.'
         call mpi_comm_size(self % mpi % comm, self % mpi % np, ie)
         self % mpi_transpose_needed = self % mpi % np > 1
       case(3)
         call mpi_cart_sub(d3d % mpi % comm, [.true., .false.], self % mpi % comm, ie)
         if (ie /= 0) stop 'error executing mpi_cart_sub.'
         call mpi_comm_size(self % mpi % comm, self % mpi % np, ie)
         self % mpi_transpose_needed = self % mpi % np > 1
      end select
   else
      stop 'not implemented.'
   end if
   call mpi_comm_rank(self % mpi % comm, self % mpi % rank, ie)
#endif

   select case(direction)
    case(1)
      self % nx = d3d % nx
      self % gnx = d3d % gnx
      self % offx = d3d % offx
      self % bcs = d3d % bcs(1:2)
    case(2)
      self % nx = d3d % ny
      self % gnx = d3d % gny
      self % offx = d3d % offy
      self % bcs = d3d % bcs(3:4)
    case(3)
      self % nx = d3d % nz
      self % gnx = d3d % gnz
      self % offx = d3d % offz
      self % bcs = d3d % bcs(5:6)
   end select

   self % nxyz = [self % nx]
   self % cnt = self % nx
   self % gcnt = int(self % gnx, kind(self % gcnt))
end function

function poisfft_solver2d_from3d(d3d, direction) result(self)
   type(poisfft_solver2d) :: self
   type(poisfft_solver3d), intent(in) :: d3d
   integer, intent(in) :: direction
#ifdef MPI
   integer :: ie, dims

   call mpi_cartdim_get(d3d % mpi % comm, dims, ie)
   if (ie /= 0) stop 'error executing mpi_cartdim_get.'

   if (dims == 2) then
      ! we see the dimensions reversed in fortran!
      select case(direction)
       case(1)
         self % mpi % comm = d3d % mpi % comm
         self % mpi % comm_dim = 2
       case(2)
         call mpi_cart_sub(d3d % mpi % comm, [.true., .false.], self % mpi % comm, ie)
         if (ie /= 0) stop 'error executing mpi_cart_sub.'
         self % mpi % comm_dim = 1
       case(3)
         call mpi_cart_sub(d3d % mpi % comm, [.false., .true.], self % mpi % comm, ie)
         if (ie /= 0) stop 'error executing mpi_cart_sub.'
         self % mpi % comm_dim = 1
      end select
   else
      stop 'not implemented.'
   end if

#endif
   select case(direction)
    case(1)
      self % nx = d3d % ny
      self % gnx = d3d % gny
      self % offx = d3d % offy
      self % ny = d3d % nz
      self % gny = d3d % gnz
      self % offy = d3d % offz
      self % bcs = d3d % bcs(3:6)
    case(2)
      self % nx = d3d % nx
      self % gnx = d3d % gnx
      self % offx = d3d % offx
      self % ny = d3d % nz
      self % gny = d3d % gnz
      self % offy = d3d % offz
      self % bcs = d3d % bcs([1, 2, 5, 6])
    case(3)
      self % nx = d3d % nx
      self % gnx = d3d % gnx
      self % offx = d3d % offx
      self % ny = d3d % ny
      self % gny = d3d % gny
      self % offy = d3d % offy
      self % bcs = d3d % bcs(1:4)
   end select

   self % nxyz = [self % nx, self % ny]
   self % cnt = self % nx * self % ny
   self % gcnt = int(self % gnx, kind(self % gcnt)) * int(self % gny, kind(self % gcnt))
end function

function poisfft_solver2d__new(nxyz, lxyz, bcs, approximation, gnxyz, offs, mpi_comm, nthreads) result(self)
   type(poisfft_solver2d) :: self

   integer, intent(in) :: nxyz(2)
   real(RP), intent(in) :: lxyz(2)
   integer, intent(in) :: bcs(4)
   integer, intent(in), optional :: approximation, gnxyz(2), offs(2), mpi_comm, nthreads

   self % nxyz = nxyz
   self % lx = lxyz(1)
   self % ly = lxyz(2)
   self % nx = nxyz(1)
   self % ny = nxyz(2)

   if (present(gnxyz)) then
      self % gnx = gnxyz(1)
      self % gny = gnxyz(2)
   else
      self % gnx = self % nx
      self % gny = self % ny
   end if

   if (present(offs)) then
      self % offx = offs(1)
      self % offy = offs(2)
   end if

   self % cnt = product(self % nxyz)
   self % gcnt = product(int([self % gnx, self % gny], kind(self % gcnt)))
   self % bcs = bcs

   if (present(approximation)) self % approximation = approximation

   if (present(mpi_comm)) then
      self % mpi % comm = mpi_comm
#ifdef MPI
   else
      stop 'no pfft comm present in poisfft_solver2d__new.'
#endif
   end if

   self % nthreads = merge(nthreads, 1, present(nthreads))
   call init(self) ! create fftw plans and allocate working array
end function

subroutine poisfft_solver2d_init(self)
   type(poisfft_solver2d), intent(inout) :: self
   integer :: i

   self % norm_factor = norm_factor(self % gnx, self % bcs(1:2)) * norm_factor(self % gny, self % bcs(3:4))
   allocate(self % denomx(self % nx), self % denomy(self % ny))

   if (all(self % bcs == poisfft_periodic)) then
      call allocate_fftw_complex(self)
      self % forward = poisfft_plan2d(self, [fft_complex, fftw_forward])
      self % backward = poisfft_plan2d(self, [fft_complex, fftw_backward])
   else if (all(self % bcs == poisfft_dirichlet)) then
      call allocate_fftw_real(self)
      self % forward = poisfft_plan2d(self, [(fft_realodd00, i=1, 2)])
      self % backward = poisfft_plan2d(self, [(fft_realodd00, i=1, 2)])
   else if (all(self % bcs == poisfft_dirichletstag)) then
      call allocate_fftw_real(self)
      self % forward = poisfft_plan2d(self, [(fft_realodd10, i=1, 2)])
      self % backward = poisfft_plan2d(self, [(fft_realodd01, i=1, 2)])
   else if (all(self % bcs == poisfft_neumann)) then
      call allocate_fftw_real(self)
      self % forward = poisfft_plan2d(self, [(fft_realeven00, i=1, 2)])
      self % backward = poisfft_plan2d(self, [(fft_realeven00, i=1, 2)])
   else if (all(self % bcs == poisfft_neumannstag)) then
      call allocate_fftw_real(self)
      self % forward = poisfft_plan2d(self, [(fft_realeven10, i=1, 2)])
      self % backward = poisfft_plan2d(self, [(fft_realeven01, i=1, 2)])
   endif

   if (self % approximation == poisfft_finitedifference2) then
      self % denomx = eigenvalues(eig_fn_fd2, self % bcs(1:2), self % lx, self % nx, self % gnx, self % offx)
      self % denomy = eigenvalues(eig_fn_fd2, self % bcs(2:4), self % ly, self % ny, self % gny, self % offy)
   else if (self % approximation == poisfft_finitedifference4) then
      self % denomx = eigenvalues(eig_fn_fd4, self % bcs(1:2), self % lx, self % nx, self % gnx, self % offx)
      self % denomy = eigenvalues(eig_fn_fd4, self % bcs(2:4), self % ly, self % ny, self % gny, self % offy)
   else
      self % denomx = eigenvalues(eig_fn_spectral, self % bcs(1:2), self % lx, self % nx, self % gnx, self % offx)
      self % denomy = eigenvalues(eig_fn_spectral, self % bcs(2:4), self % ly, self % ny, self % gny, self % offy)
   end if
end subroutine

subroutine poisfft_solver2d__execute(self, phi, rhs)
   type(poisfft_solver2d), intent(inout) :: self
   real(RP), intent(out) :: phi(:, :)
   real(RP), intent(in) :: rhs(:, :)
   integer :: ngphi(2), ngrhs(2)

   ngphi = (ubound(phi) - [self % nx, self % ny]) / 2
   ngrhs = (ubound(rhs) - [self % nx, self % ny]) / 2
   if (all(self % bcs == poisfft_periodic)) then
      call poisfft_solver2d_fullperiodic(self, &
         phi(ngphi(1) + 1:ngphi(1) + self % nx, ngphi(2) + 1:ngphi(2) + self % ny), &
         rhs(ngrhs(1) + 1:ngrhs(1) + self % nx, ngrhs(2) + 1:ngrhs(2) + self % ny))
   else if (all(self % bcs == poisfft_dirichlet) .or. all(self % bcs == poisfft_dirichletstag)) then
      call poisfft_solver2d_fulldirichlet(self, &
         phi(ngphi(1) + 1:ngphi(1) + self % nx, ngphi(2) + 1:ngphi(2) + self % ny), &
         rhs(ngrhs(1) + 1:ngrhs(1) + self % nx, ngrhs(2) + 1:ngrhs(2) + self % ny))
   else if (all(self % bcs == poisfft_neumann) .or. all(self % bcs == poisfft_neumannstag)) then
      call poisfft_solver2d_fullneumann(self, &
         phi(ngphi(1) + 1:ngphi(1) + self % nx, ngphi(2) + 1:ngphi(2) + self % ny), &
         rhs(ngrhs(1) + 1:ngrhs(1) + self % nx, ngrhs(2) + 1:ngrhs(2) + self % ny))
   endif
end subroutine

subroutine poisfft_solver2d__finalize(self)
   type(poisfft_solver2d), intent(inout) :: self
   call finalize(self % forward)
   call finalize(self % backward)
   if (associated(self % rwork)) call deallocate_fftw(self % rwork)
   if (associated(self % cwork)) call deallocate_fftw(self % cwork)
end subroutine

function poisfft_solver1d__new(nxyz, lxyz, bcs, approximation, gnxyz, offs, mpi_comm, nthreads) result(self)
   type(poisfft_solver1d) :: self

   integer, intent(in) :: nxyz(1)
   real(RP), intent(in) :: lxyz(1)
   integer, intent(in) :: bcs(2)
   integer, intent(in), optional :: approximation, gnxyz(1), offs(1), mpi_comm, nthreads
#ifdef MPI
   integer :: ie, dims
#endif
      end if
      
      
      
     if (present(nthreads)) then
        D%nthreads = nthreads
      else
        D%nthreads = 1
      endif

       !create fftw plans and allocate working array
      call Init(D)
    end function PoisFFT_Solver1D__New




    subroutine Poisfft_Solver1D_Init(D)
      type(PoisFFT_Solver1D), intent(inout) :: D
      
      D%norm_factor = norm_factor(D%gnx,D%BCs(1:2))
      
      allocate(D%denomx(D%gnx))

      if (all(D%BCs==PoisFFT_Periodic)) then

        call allocate_fftw_complex(D)

        D%forward = PoisFFT_Plan1D(D, [FFT_Complex, FFTW_FORWARD])
        D%backward = PoisFFT_Plan1D(D, [FFT_Complex, FFTW_BACKWARD])

      else if (all(D%BCs==PoisFFT_Dirichlet)) then

        call allocate_fftw_real(D)

        D%forward = PoisFFT_Plan1D(D, [FFT_RealOdd00])
        D%backward = PoisFFT_Plan1D(D, [FFT_RealOdd00])
        
      else if (all(D%BCs==PoisFFT_DirichletStag)) then

        call allocate_fftw_real(D)

        D%forward = PoisFFT_Plan1D(D, [FFT_RealOdd10])
        D%backward = PoisFFT_Plan1D(D, [FFT_RealOdd01])
        
      else if (all(D%BCs==PoisFFT_Neumann)) then

        call allocate_fftw_real(D)

        D%forward = PoisFFT_Plan1D(D, [FFT_RealEven00])
        D%backward = PoisFFT_Plan1D(D, [FFT_RealEven00])
        
      else if (all(D%BCs==PoisFFT_NeumannStag)) then

        call allocate_fftw_real(D)

        D%forward = PoisFFT_Plan1D(D, [FFT_RealEven10])
        D%backward = PoisFFT_Plan1D(D, [FFT_RealEven01])
        
      else if (all(D%BCs==[PoisFFT_NeumannStag, PoisFFT_DirichletStag])) then
      
        call allocate_fftw_real(D)

        D%forward = PoisFFT_Plan1D(D, [FFT_RealEven11])
        D%backward = PoisFFT_Plan1D(D, [FFT_RealEven11])
        
      else if (all(D%BCs==[PoisFFT_DirichletStag, PoisFFT_NeumannStag])) then
      
        call allocate_fftw_real(D)

        D%forward = PoisFFT_Plan1D(D, [FFT_RealOdd11])
        D%backward = PoisFFT_Plan1D(D, [FFT_RealOdd11])
        
      endif
      
      if (D%approximation==PoisFFT_FiniteDifference2) then
        D%denomx = eigenvalues(eig_fn_FD2, D%BCs(1:2), D%Lx, D%gnx, D%gnx, D%offx)
      else if (D%approximation==PoisFFT_FiniteDifference4) then
        D%denomx = eigenvalues(eig_fn_FD4, D%BCs(1:2), D%Lx, D%gnx, D%gnx, D%offx)
      else
        D%denomx = eigenvalues(eig_fn_spectral, D%BCs(1:2), D%Lx, D%gnx, D%gnx, D%offx)
      end if

    end subroutine PoisFFT_Solver1D_Init
    
    
    
    
    subroutine PoisFFT_Solver1D__Execute(D,Phi,RHS)
      type(PoisFFT_Solver1D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:)
      real(RP), intent(in)  :: RHS(:)

      integer   :: ngPhi(1), ngRHS(1)


      ngPhi = (ubound(Phi)-D%nx)/2
      ngRHS = (ubound(RHS)-D%nx)/2


      if (all(D%BCs==PoisFFT_Periodic)) then

        call PoisFFT_Solver1D_FullPeriodic(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx))

      else if (all(D%BCs==PoisFFT_Dirichlet) .or. &
               all(D%BCs==PoisFFT_DirichletStag)) then

        call PoisFFT_Solver1D_FullDirichlet(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx))

      else if (all(D%BCs==PoisFFT_Neumann) .or. &
               all(D%BCs==PoisFFT_NeumannStag)) then

        call PoisFFT_Solver1D_FullNeumann(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx))
                 
      else if (all(D%BCs==PoisFFT_Dirichlet .or. &
                   D%BCs==PoisFFT_Neumann .or. &
                   D%BCs==PoisFFT_DirichletStag .or. &
                   D%BCs==PoisFFT_NeumannStag) &
               .and. &
               any(D%BCs==PoisFFT_DirichletStag .or. &
                   D%BCs==PoisFFT_Dirichlet)) then
   
        call PoisFFT_Solver1D_FullDirichlet(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx))
                 
      endif

    end subroutine PoisFFT_Solver1D__Execute




    subroutine PoisFFT_Solver1D__Finalize(D)
      type(PoisFFT_Solver1D), intent(inout) :: D

      call Finalize(D%forward)
      call Finalize(D%backward)

      if (associated(D%rwork)) call deallocate_fftw(D%rwork)
      if (associated(D%cwork)) call deallocate_fftw(D%cwork)

    endsubroutine PoisFFT_Solver1D__Finalize

   if (present(gnxyz)) then
      self % gnx = gnxyz(1)
   else
      self % gnx = self % nx
   end if

   if (present(offs)) self % offx = offs(1)
   self % cnt = self % nx
   self % gcnt = int(self % gnx, kind(self % gcnt))
   self % bcs = bcs

   if (present(approximation)) self % approximation = approximation

   if (present(mpi_comm)) then
      self % mpi % comm = mpi_comm
#ifdef MPI
      call MPI_Cartdim_get(self % mpi % comm, dims, ie)
      if (ie /= 0) stop 'error executing mpi_cartdim_get.'
      self % mpi_transpose_needed = dims > 0
   else
      stop 'no pfft comm present in poisfft_solver1d__new.'
#endif
   end if

   self % nthreads = merge(nthreads, 1, present(nthreads))
   call init(self) ! create fftw plans and allocate working array
end function

subroutine poisfft_solver1d_init(self)
   type(poisfft_solver1d), intent(inout) :: self

   self % norm_factor = norm_factor(self % gnx, self % bcs(1:2))
   allocate(self % denomx(self % gnx))

   if (all(self % bcs == poisfft_periodic)) then
      call allocate_fftw_complex(self)
      self % forward = poisfft_plan1d(self, [fft_complex, fftw_forward])
      self % backward = poisfft_plan1d(self, [fft_complex, fftw_backward])
   else if (all(self % bcs == poisfft_dirichlet)) then
      call allocate_fftw_real(self)
      self % forward = poisfft_plan1d(self, [fft_realodd00])
      self % backward = poisfft_plan1d(self, [fft_realodd00])
   else if (all(self % bcs == poisfft_dirichletstag)) then
      call allocate_fftw_real(self)
      self % forward = poisfft_plan1d(self, [fft_realodd10])
      self % backward = poisfft_plan1d(self, [fft_realodd01])
   else if (all(self % bcs == poisfft_neumann)) then
      call allocate_fftw_real(self)
      self % forward = poisfft_plan1d(self, [fft_realeven00])
      self % backward = poisfft_plan1d(self, [fft_realeven00])
   else if (all(self % bcs == poisfft_neumannstag)) then
      call allocate_fftw_real(self)
      self % forward = poisfft_plan1d(self, [fft_realeven10])
      self % backward = poisfft_plan1d(self, [fft_realeven01])
   endif

   if (self % approximation == poisfft_finitedifference2) then
      self % denomx = eigenvalues(eig_fn_fd2, self % bcs(1:2), self % lx, self % gnx, self % gnx, self % offx)
   else if (self % approximation == poisfft_finitedifference4) then
      self % denomx = eigenvalues(eig_fn_fd4, self % bcs(1:2), self % lx, self % gnx, self % gnx, self % offx)
   else
      self % denomx = eigenvalues(eig_fn_spectral, self % bcs(1:2), self % lx, self % gnx, self % gnx, self % offx)
   end if
end subroutine

subroutine poisfft_solver1d__execute(self, phi, rhs)
   type(poisfft_solver1d), intent(inout) :: self
   real(RP), intent(out) :: phi(:)
   real(RP), intent(in) :: rhs(:)

   integer :: ngphi(1), ngrhs(1)

   ngphi = (ubound(phi) - self % nx) / 2
   ngrhs = (ubound(rhs) - self % nx) / 2
   if (all(self % bcs == poisfft_periodic)) then
      call poisfft_solver1d_fullperiodic(self, &
         phi(ngphi(1) + 1:ngphi(1) + self % nx), &
         rhs(ngrhs(1) + 1:ngrhs(1) + self % nx))
   else if (all(self % bcs == poisfft_dirichlet) .or. all(self % bcs == poisfft_dirichletstag)) then
      call poisfft_solver1d_fulldirichlet(self, &
         phi(ngphi(1) + 1:ngphi(1) + self % nx), &
         rhs(ngrhs(1) + 1:ngrhs(1) + self % nx))
   else if (all(self % bcs == poisfft_neumann) .or. all(self % bcs == poisfft_neumannstag)) then
      call poisfft_solver1d_fullneumann(self, &
         phi(ngphi(1) + 1:ngphi(1) + self % nx), &
         rhs(ngrhs(1) + 1:ngrhs(1) + self % nx))
   endif
end subroutine

subroutine poisfft_solver1d__finalize(self)
   type(poisfft_solver1d), intent(inout) :: self
   call finalize(self % forward)
   call finalize(self % backward)
   if (associated(self % rwork)) call deallocate_fftw(self % rwork)
   if (associated(self % cwork)) call deallocate_fftw(self % cwork)
end subroutine

    
    
    
    
    
    
    
    
    function PoisFFT_Solver3D_nonuniform_z__New(nxyz,Lxy,z,z_u, &
                                  BCs,approximation, &
                                  gnxyz,offs,mpi_comm,nthreads,ierr) result(D)
      type(PoisFFT_Solver3D_nonuniform_z) :: D

      integer, intent(in)   :: nxyz(3)
      real(RP), intent(in)  :: Lxy(2)
      real(RP), intent(in)  :: z(:), z_u(0:)
      integer, intent(in)   :: bcs(6)
      integer, intent(in), optional :: approximation
      integer, intent(in), optional :: gnxyz(3)
      integer, intent(in), optional :: offs(3)
      integer, intent(in), optional :: mpi_comm
      integer, intent(in), optional :: nthreads
      integer, intent(out),optional :: ierr
      
      if (present(ierr)) ierr = 0

      D%nxyz = nxyz

      D%Lx = Lxy(1)
      D%Ly = Lxy(2)
      D%Lz = 0

      D%nx = nxyz(1)
      D%ny = nxyz(2)
      D%nz = nxyz(3)
      
      allocate(D%z(1:ubound(z,1)))
      D%z = z
      allocate(D%z_u(0:ubound(z_u,1)))
      D%z_u(:) = z_u
      
      if (present(gnxyz)) then
        D%gnx = gnxyz(1)
        D%gny = gnxyz(2)
        D%gnz = gnxyz(3)
      else
        D%gnx = D%nx
        D%gny = D%ny
        D%gnz = D%nz
      end if

      if (present(offs)) then
        D%offx = offs(1)
        D%offy = offs(2)
        D%offz = offs(3)
      end if

      D%cnt = product(D%nxyz)
      
      D%gcnt = product(int([D%gnx, D%gny, D%gnz], kind(D%gcnt)))

      D%BCs = BCs
      
      if (present(approximation)) D%approximation = approximation
      
      if (present(mpi_comm)) then
        D%mpi%comm = mpi_comm
#ifdef MPI      
      else
        stop "No PFFT comm present in PoisFFT_Solver3D__New."
#endif
      end if

      if (present(nthreads)) then
        D%nthreads = nthreads
      else
        D%nthreads = 1
      endif

      !create fftw plans and allocate working arrays
      call Init(D)
    end function PoisFFT_Solver3D_nonuniform_z__New
    
    
    
    
    subroutine PoisFFT_Solver3D_nonuniform_z_Init(D, ierr, errmsg)
      type(PoisFFT_Solver3D_nonuniform_z), intent(inout) :: D
      integer,   intent(out), optional :: ierr
      character, intent(out), optional :: errmsg
      integer :: real_forw, real_back
      integer :: real_forw_xy, real_back_xy
      integer :: real_forw_z, real_back_z
      integer :: i
      !$omp parallel
      !$omp single
      !$ D%nthreads = omp_get_num_threads()
      !$omp end single
      !$omp end parallel

      D%norm_factor = norm_factor(D%gnx,D%BCs(1:2)) * &
                      norm_factor(D%gny,D%BCs(3:4))
      
      allocate(D%denomx(D%nx))
      allocate(D%denomy(D%ny))
      
      
#ifdef MPI
      call PoisFFT_PFFT_init()
#else
      if (D%nthreads>1) call PoisFFT_InitThreads(1)
#endif

      if (all(D%BCs==PoisFFT_Neumann) .or. &
               all(D%BCs==PoisFFT_NeumannStag)) then
#ifdef MPI               
        if (D%nthreads>1) call PoisFFT_PFFT_InitThreads(D%nthreads)
#endif
        real_forw = real_transform_type_forward(D%BCs(1:2))
        real_back = real_transform_type_backward(D%BCs(1:2))        

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

        call allocate_fftw_real(D%Solvers2D(1))

        D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [real_forw, real_forw])
        D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [real_back, real_back])

        do i=2,D%nthreads
          D%Solvers2D(i) = D%Solvers2D(1)
          D%Solvers2D(i)%forward%planowner = .false.
          D%Solvers2D(i)%backward%planowner = .false.

          call allocate_fftw_real(D%Solvers2D(i))

        end do
        
        allocate(D%Solvers1D(3))
        D%Solvers1D(3) = PoisFFT_Solver1D_From3D(D,3)


      else if (all(D%BCs(1:4)==PoisFFT_Periodic) .and. &
               (all(D%BCs(5:6)==PoisFFT_Neumann) .or. &
                all(D%BCs(5:6)==PoisFFT_NeumannStag) )) then

#ifdef MPI
        call allocate_fftw_complex(D)
        
        allocate(D%Solvers1D(3:2+D%nthreads))

        D%Solvers1D(3) = PoisFFT_Solver1D_From3D(D,3)

        call allocate_fftw_real(D%Solvers1D(3))

        real_forw = real_transform_type_forward(D%BCs(5:6))
        real_back = real_transform_type_backward(D%BCs(5:6))
        D%Solvers1D(3)%forward = PoisFFT_Plan1D(D%Solvers1D(3), [real_forw, real_forw])
        D%Solvers1D(3)%backward = PoisFFT_Plan1D(D%Solvers1D(3), [real_back, real_back])

        do i = 4, 2+D%nthreads
          D%Solvers1D(i) = D%Solvers1D(3)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_real(D%Solvers1D(i))
        end do

        if (D%ny<D%gny) then
          if (D%nthreads>1) call PoisFFT_PFFT_InitThreads(D%nthreads)

          allocate(D%Solvers2D(1))

          D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D, 3)

          call allocate_fftw_complex(D%Solvers2D(1))

          D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_FORWARD])
          D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_BACKWARD])

          
        else
          allocate(D%Solvers2D(D%nthreads))

          D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

          call allocate_fftw_complex(D%Solvers2D(1))

          D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), &
                                                  [FFT_Complex, FFTW_FORWARD], &
                                                  distributed=.false.)
          D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), &
                                                   [FFT_Complex, FFTW_BACKWARD], &
                                                   distributed=.false.)

          do i = 2, D%nthreads
            D%Solvers2D(i) = D%Solvers2D(1)
            D%Solvers2D(i)%forward%planowner = .false.
            D%Solvers2D(i)%backward%planowner = .false.

            call allocate_fftw_complex(D%Solvers2D(i))
          end do
        end if

        call Init_MPI_Buffers(D, 3)

#else
        call allocate_fftw_complex(D)

        allocate(D%Solvers1D(D%nthreads))

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

        call allocate_fftw_complex(D%Solvers2D(1))

        D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_FORWARD])
        D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_BACKWARD])

        do i = 2, D%nthreads
          D%Solvers2D(i) = D%Solvers2D(1)
          D%Solvers2D(i)%forward%planowner = .false.
          D%Solvers2D(i)%backward%planowner = .false.

          call allocate_fftw_complex(D%Solvers2D(i))
        end do
        
!MPI
#endif


      else
        stop "Unknown combination of boundary conditions."
      endif

#ifdef MPI
      call Init_MPI_Buffers(D, 3)
#endif   


      if (D%approximation==PoisFFT_FiniteDifference2) then
        D%denomx = eigenvalues(eig_fn_FD2, D%BCs(1:2), D%Lx, D%nx, D%gnx, D%offx)
        D%denomy = eigenvalues(eig_fn_FD2, D%BCs(3:4), D%Ly, D%ny, D%gny, D%offy)
      else if (D%approximation==PoisFFT_FiniteDifference4) then
        D%denomx = eigenvalues(eig_fn_FD4, D%BCs(1:2), D%Lx, D%gnx, D%gnx, D%offx)
        D%denomy = eigenvalues(eig_fn_FD4, D%BCs(3:4), D%Ly, D%gny, D%gny, D%offy)
      else
        D%denomx = eigenvalues(eig_fn_spectral, D%BCs(1:2), D%Lx, D%nx, D%gnx, D%offx)
        D%denomy = eigenvalues(eig_fn_spectral, D%BCs(3:4), D%Ly, D%ny, D%gny, D%offy)
      end if
      
      allocate(D%mat_a(2:D%gnz))
      allocate(D%mat_b(1:D%gnz))
      allocate(D%mat_c(1:max(1,D%gnz-1)))
    
      do i = 2, D%gnz
        D%mat_a(i) = 1 / (D%z(i)-D%z(i-1)) / (D%z_u(i) - D%z_u(i-1))
      end do
      i = 1
      if (D%gnz==1) then
        D%mat_b(1) = 0
      else
        D%mat_b(1) = -1 / (D%z(i+1)-D%z(i)) / (D%z_u(i) - D%z_u(i-1))
      end if
      do i = 2, D%gnz-1
        D%mat_b(i) = -1 / (D%z(i+1)-D%z(i)) / (D%z_u(i) - D%z_u(i-1))
        D%mat_b(i) = D%mat_b(i) - 1 / (D%z(i)-D%z(i-1)) / (D%z_u(i) - D%z_u(i-1))   
      end do
      if (D%gnz>1) then
        i = D%gnz
        D%mat_b(D%gnz) = -1 / (D%z(i)-D%z(i-1)) / (D%z_u(i) - D%z_u(i-1))
      end if
      
      do i = 1, D%gnz-1
        D%mat_c(i) = 1 / (D%z(i+1)-D%z(i)) / (D%z_u(i) - D%z_u(i-1))
      end do
    end subroutine PoisFFT_Solver3D_nonuniform_z_Init
    
    


    subroutine PoisFFT_Solver3D_nonuniform_z__Execute(D,Phi,RHS)
      type(PoisFFT_Solver3D_nonuniform_z), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)

      integer   :: ngPhi(3), ngRHS(3)

      ngPhi = (ubound(Phi)-[D%nx,D%ny,D%nz])/2
      ngRHS = (ubound(RHS)-[D%nx,D%ny,D%nz])/2

      if (all(D%BCs==PoisFFT_Neumann) .or. &
               all(D%BCs==PoisFFT_NeumannStag)) then

        call PoisFFT_Solver3D_nonuniform_z_FullNeumann(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))

      else if (all(D%BCs(1:4)==PoisFFT_Periodic) .and. &
               (all(D%BCs(5:6)==PoisFFT_Neumann) .or. &
                all(D%BCs(5:6)==PoisFFT_NeumannStag))) then
                
        call PoisFFT_Solver3D_nonuniform_z_PPNs(D,&
                Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                    ngPhi(2)+1:ngPhi(2)+D%ny,&
                    ngPhi(3)+1:ngPhi(3)+D%nz),&
                RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                    ngRHS(2)+1:ngRHS(2)+D%ny,&
                    ngRHS(3)+1:ngRHS(3)+D%nz))

      endif

    end subroutine PoisFFT_Solver3D_nonuniform_z__Execute
    
    
    
    
    subroutine PoisFFT_Solver3D_nonuniform_z__Finalize(D)
      type(PoisFFT_Solver3D_nonuniform_z), intent(inout) :: D
      integer :: i

      call Finalize(D%PoisFFT_Solver3D)
      
      deallocate(D%z, D%z_u)
      
      deallocate(D%mat_a, D%mat_b, D%mat_c)

    endsubroutine PoisFFT_Solver3D_nonuniform_z__Finalize
    
    



    
    
    
    
    
    
    
    
    
#include "poisfft-solvers-inc.f90"

#include "poisfft-solvers-nonuniform_z-inc.f90"


#ifdef MPI
subroutine init_mpi_buffers(self, dir)
   interface
      subroutine mpi_alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierror)
         integer :: sendbuf(*), recvbuf(*)
         integer :: sendcount, sendtype, recvcount, recvtype
         integer :: comm, ierror
      end subroutine
   end interface

   type(poisfft_solver3d), intent(inout), target :: self
   integer, intent(in) :: dir
   type(mpi_vars_1d), pointer :: mpi
   integer :: i, ie

   mpi => self % solvers1d(dir) % mpi

   allocate(mpi % snxs(mpi % np), mpi % snzs(mpi % np), mpi % rnxs(mpi % np), mpi % rnzs(mpi % np))
   allocate(mpi % sumrnzs(mpi % np))
   allocate(mpi % sdispls(mpi % np), mpi % scounts(mpi % np), mpi % rdispls(mpi % np), mpi % rcounts(mpi % np))

   if (dir == 2) then
      mpi % snxs(1:mpi % np - 1) = (self % nx / mpi % np)
      mpi % snxs(mpi % np) = self % nx - sum(mpi % snxs(1:mpi % np - 1))
      mpi % snzs = self % ny
      mpi % scounts = mpi % snxs * self % ny * self % nz
      mpi % sdispls(1) = 0
      do i = 2, mpi % np
         mpi % sdispls(i) = mpi % sdispls(i - 1) + mpi % scounts(i - 1)
      end do

      call mpi_alltoall(mpi % snxs, 1, mpi_integer, mpi % rnxs, 1, mpi_integer, mpi % comm, ie)
      if (.not. all(mpi % rnxs(2:mpi % np) == mpi % rnxs(1))) stop 'poisfft internal error: .not.all(mpi%rnxs(2:mpi%np)==mpi%rnxs(1))'
      call mpi_alltoall(mpi % snzs, 1, mpi_integer, mpi % rnzs, 1, mpi_integer, mpi % comm, ie)
      call mpi_alltoall(mpi % scounts, 1, mpi_integer, mpi % rcounts, 1, mpi_integer, mpi % comm, ie)

      mpi % rdispls(1) = 0
      do i = 2, mpi % np; mpi % rdispls(i) = mpi % rdispls(i - 1) + mpi % rcounts(i - 1); end do
      do i = 1, mpi % np; mpi % sumrnzs(i) = sum(mpi % rnzs(1:i - 1)); end do

      allocate(mpi % tmp1(1:self % ny, 1:self % nz, 1:self % nx))
      allocate(mpi % tmp2(0:sum(mpi % rcounts) - 1))
      allocate(mpi % rwork(sum(mpi % rnzs), self % nz, mpi % rnxs(1)))
   else if (dir == 3) then
      mpi % snxs(1:mpi % np - 1) = (self % nx / mpi % np)
      mpi % snxs(mpi % np) = self % nx - sum(mpi % snxs(1:mpi % np - 1))
      mpi % snzs = self % nz
      mpi % scounts = mpi % snxs * self % ny * self % nz
      mpi % sdispls(1) = 0
      do i = 2, mpi % np
         mpi % sdispls(i) = mpi % sdispls(i - 1) + mpi % scounts(i - 1)
      end do

      call mpi_alltoall(mpi % snxs, 1, mpi_integer, mpi % rnxs, 1, mpi_integer, mpi % comm, ie)
      if (.not. all(mpi % rnxs(2:mpi % np) == mpi % rnxs(1))) stop 'poisfft internal error: .not.all(mpi%rnxs(2:mpi%np)==mpi%rnxs(1))'
      call mpi_alltoall(mpi % snzs, 1, mpi_integer, mpi % rnzs, 1, mpi_integer, mpi % comm, ie)
      call mpi_alltoall(mpi % scounts, 1, mpi_integer, mpi % rcounts, 1, mpi_integer, mpi % comm, ie)

      mpi % rdispls(1) = 0
      do i = 2, mpi % np
         mpi % rdispls(i) = mpi % rdispls(i - 1) + mpi % rcounts(i - 1)
      end do

      do i = 1, mpi % np
         mpi % sumrnzs(i) = sum(mpi % rnzs(1:i - 1))
      end do

      allocate(mpi % tmp1(1:self % nz, 1:self % ny, 1:self % nx))
      allocate(mpi % tmp2(0:sum(mpi % rcounts) - 1))
      allocate(mpi % rwork(sum(mpi % rnzs), self % ny, mpi % rnxs(1)))
   else
      stop 'not implemented.'
   end if
end subroutine




    function eigenvalues(f, BCs, L, n, gn, off) result(res)
      procedure(eig_fn_spectral) :: f
      integer, intent(in) :: BCs(2)
      real(RP), intent(in) :: L
      integer, intent(in) :: n, gn, off
      real(RP) :: res(n)
      integer :: i
      real(RP) :: dkx, dkx_h
      
      dkx = pi * grid_dx(gn, L, BCs) / L
      dkx_h = dkx / 2
      
      if (all(BCs==PoisFFT_Periodic)) then
        do i = 1, n
          if (i+off<gn/2) then
            res(i) = f((i-1+off)*dkx)
          else
            res(i) = f((gn-(i-1+off))*dkx)
          end if
        end do
      else if (all(BCs==PoisFFT_Dirichlet)) then
        forall(i=1:n) res(i) = f((i+off)*dkx_h)
      else if (all(BCs==PoisFFT_DirichletStag)) then
        forall(i=1:n) res(i) = f((i+off)*dkx_h)
      else if (all(BCs==PoisFFT_Neumann)) then
        forall(i=1:n) res(i) = f((i-1+off)*dkx_h)
      else if (all(BCs==PoisFFT_NeumannStag)) then
        forall(i=1:n) res(i) = f((i-1+off)*dkx_h)
      else if (all(BCs==[PoisFFT_NeumannStag,PoisFFT_DirichletStag])) then
        forall(i=1:n) res(i) = f((i+off-0.5_RP)*dkx_h)
      else if (all(BCs==[PoisFFT_DirichletStag,PoisFFT_NeumannStag])) then
        forall(i=1:n) res(i) = f((i+off-0.5_RP)*dkx_h)
      end if

      res = res / (grid_dx(gn, L, BCs))**2
      
    end function


    pure real(RP) function eig_fn_spectral(x) result(f)
      real(RP), intent(in) :: x
      f = -(2 * x)**2
    end function
  
    pure real(RP) function eig_fn_FD2(x) result(f)
      real(RP), intent(in) :: x
      f = -(2 * sin(x))**2
    end function
  
    pure real(RP) function eig_fn_FD4(x) result(f)
      real(RP), intent(in) :: x
      f = -((sin(3*x) - 27*sin(x)) / 12)**2
    end function
    
      
    function grid_dx(gn, L, BCs) result(res)
      real(RP) :: res
      integer, intent(in) :: gn
      real(RP), intent(in) :: L
      integer, intent(in) :: BCs(2)
      
      if (all(BCs==PoisFFT_Periodic)) then
        res = L / gn
      else if (all(BCs==PoisFFT_Dirichlet)) then
        res = L / (gn+1)
      else if (all(BCs==PoisFFT_DirichletStag .or. &
                   BCs==PoisFFT_NeumannStag)) then
        res = L / gn
      else if (all(BCs==PoisFFT_Neumann)) then
        res = L / (gn-1)
      else
        res = 0
      end if
    end function

    function norm_factor(gn,BCs) result(res)
      real(RP) :: res
      integer, intent(in) :: gn
      integer, intent(in) :: BCs(2)
      
      if (all(BCs==PoisFFT_Periodic)) then
        res = gn
      else if (all(BCs==PoisFFT_Dirichlet)) then
        res = 2 * (gn+1)
      else if (all(BCs==PoisFFT_Neumann)) then
        res = 2 * (gn-1)
      else if (all(BCs==PoisFFT_DirichletStag .or. &
                   BCs==PoisFFT_NeumannStag)) then
        res = 2 * gn
      else
        res = 0
      end if
    end function
    
    function real_transform_type_forward(BCs) result(res)
      integer :: res
      integer, intent(in) :: BCs(2)
      
      if (all(BCs==PoisFFT_Dirichlet)) then
          res = FFT_RealOdd00
      else if (all(BCs==PoisFFT_DirichletStag)) then
          res = FFT_RealOdd10
      else if (all(BCs==PoisFFT_Neumann)) then
          res = FFT_RealEven00
      else if (all(BCs==PoisFFT_NeumannStag)) then
          res = FFT_RealEven10
      else if (all(BCs==[PoisFFT_NeumannStag,PoisFFT_DirichletStag])) then
          res = FFT_RealEven11
      else if (all(BCs==[PoisFFT_DirichletStag,PoisFFT_NeumannStag])) then
          res = FFT_RealOdd11
      else
          res = -1
      end if
    end function
    
    function real_transform_type_backward(BCs) result(res)
      integer :: res
      integer, intent(in) :: BCs(2)
      
      if (all(BCs==PoisFFT_Dirichlet)) then
          res = FFT_RealOdd00
      else if (all(BCs==PoisFFT_DirichletStag)) then
          res = FFT_RealOdd01
      else if (all(BCs==PoisFFT_Neumann)) then
          res = FFT_RealEven00
      else if (all(BCs==PoisFFT_NeumannStag)) then
          res = FFT_RealEven01
      else if (all(BCs==[PoisFFT_NeumannStag,PoisFFT_DirichletStag])) then
          res = FFT_RealEven11
      else if (all(BCs==[PoisFFT_DirichletStag,PoisFFT_NeumannStag])) then
          res = FFT_RealOdd11
      else
          res = -1
      end if
    end function
    
    

    
#ifdef MPI     
    subroutine Init_MPI_Buffers(D, dir)
      interface
        subroutine MPI_ALLTOALL(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT,&
                                RECVTYPE, COMM, IERROR)
            INTEGER   SENDBUF(*), RECVBUF(*)
            INTEGER   SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE
            INTEGER   COMM, IERROR
        end subroutine
      end interface
      
      class(PoisFFT_Solver3D), intent(inout), target :: D
      integer, intent(in) :: dir
      type(mpi_vars_1d), pointer :: mpi
      integer :: i, ie
      
      mpi => D%Solvers1D(dir)%mpi
    
      allocate(mpi%snxs(mpi%np))
      allocate(mpi%snzs(mpi%np))

      allocate(mpi%rnxs(mpi%np))
      allocate(mpi%rnzs(mpi%np))

      allocate(mpi%sumrnzs(mpi%np))

      allocate(mpi%sdispls(mpi%np))
      allocate(mpi%scounts(mpi%np))

      allocate(mpi%rdispls(mpi%np))
      allocate(mpi%rcounts(mpi%np))
        
      if (dir==2) then
        mpi%snxs(1:mpi%np-1) = (D%nx / mpi%np)
        mpi%snxs(mpi%np) = D%nx - sum(mpi%snxs(1:mpi%np-1))
        mpi%snzs = D%ny
        mpi%scounts = mpi%snxs*D%ny*D%nz
        mpi%sdispls(1) = 0
        do i = 2, mpi%np
          mpi%sdispls(i) = mpi%sdispls(i-1) + mpi%scounts(i-1)
        end do

        call MPI_AllToAll(mpi%snxs, 1, MPI_INTEGER, &
                          mpi%rnxs, 1, MPI_INTEGER, &
                          mpi%comm, ie)
        if (.not.all(mpi%rnxs(2:mpi%np)==mpi%rnxs(1))) &
          stop "PoisFFT internal error: .not.all(mpi%rnxs(2:mpi%np)==mpi%rnxs(1))"
        call MPI_AllToAll(mpi%snzs, 1, MPI_INTEGER, &
                          mpi%rnzs, 1, MPI_INTEGER, &
                          mpi%comm, ie)
        call MPI_AllToAll(mpi%scounts, 1, MPI_INTEGER, &
                          mpi%rcounts, 1, MPI_INTEGER, &
                          mpi%comm, ie)

        mpi%rdispls(1) = 0
        do i = 2, mpi%np
          mpi%rdispls(i) = mpi%rdispls(i-1) + mpi%rcounts(i-1)
        end do

        do i = 1, mpi%np
          mpi%sumrnzs(i) = sum(mpi%rnzs(1:i-1))
        end do

        allocate(mpi%tmp1(1:D%ny,1:D%nz,1:D%nx))
        allocate(mpi%tmp2(0:sum(mpi%rcounts)-1))
        allocate(mpi%rwork(sum(mpi%rnzs), D%nz, mpi%rnxs(1)))
        
      else if (dir==3) then
      
        mpi%snxs(1:mpi%np-1) = (D%nx / mpi%np)
        mpi%snxs(mpi%np) = D%nx - sum(mpi%snxs(1:mpi%np-1))
        mpi%snzs = D%nz
        mpi%scounts = mpi%snxs*D%ny*D%nz
        mpi%sdispls(1) = 0
        do i = 2, mpi%np
          mpi%sdispls(i) = mpi%sdispls(i-1) + mpi%scounts(i-1)
        end do

        call MPI_AllToAll(mpi%snxs, 1, MPI_INTEGER, &
                          mpi%rnxs, 1, MPI_INTEGER, &
                          mpi%comm, ie)
        if (.not.all(mpi%rnxs(2:mpi%np)==mpi%rnxs(1))) &
          stop "PoisFFT internal error: .not.all(mpi%rnxs(2:mpi%np)==mpi%rnxs(1))"
        call MPI_AllToAll(mpi%snzs, 1, MPI_INTEGER, &
                          mpi%rnzs, 1, MPI_INTEGER, &
                          mpi%comm, ie)
        call MPI_AllToAll(mpi%scounts, 1, MPI_INTEGER, &
                          mpi%rcounts, 1, MPI_INTEGER, &
                          mpi%comm, ie)

        mpi%rdispls(1) = 0
        do i = 2, mpi%np
          mpi%rdispls(i) = mpi%rdispls(i-1) + mpi%rcounts(i-1)
        end do

        do i = 1, mpi%np
          mpi%sumrnzs(i) = sum(mpi%rnzs(1:i-1))
        end do

        allocate(mpi%tmp1(1:D%nz,1:D%ny,1:D%nx))
        allocate(mpi%tmp2(0:sum(mpi%rcounts)-1))
        allocate(mpi%rwork(sum(mpi%rnzs), D%ny, mpi%rnxs(1)))        
      else
        stop "Not implemented."
      end if
    end subroutine Init_MPI_Buffers
#endif

#undef RP
#undef CP
#undef MPI_RP
