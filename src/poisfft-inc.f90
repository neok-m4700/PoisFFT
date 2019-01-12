#if (PREC==2)

#define RP DRP
#define CP DCP
#define MPI_RP MPI_DOUBLE_PRECISION

#define _RP _DRP

#else

#define RP SRP
#define CP SCP
#define MPI_RP MPI_REAL

#define _RP _SRP
#define _CP _SCP

#endif

use iso_c_binding
!$ use omp_lib
use poisfft_precisions
use poisfft_parameters
#ifdef MPI
use mpi
#endif
implicit none

private
public :: poisfft_solver1d, poisfft_solver2d, poisfft_solver3d, finalize, execute

real(RP), parameter, private :: pi = 4 * atan(1._RP) ! 3.141592653589793238462_RP

interface poisfft_solver1d
   module procedure poisfft_solver1d__new
end interface

interface poisfft_solver2d
   module procedure poisfft_solver2d__new
end interface

interface poisfft_solver3d
   module procedure poisfft_solver3d__new
end interface

interface finalize
   module procedure poisfft_solver1d__finalize
   module procedure poisfft_solver2d__finalize
   module procedure poisfft_solver3d__finalize
end interface

interface init
   module procedure poisfft_solver1d_init
   module procedure poisfft_solver2d_init
   module procedure poisfft_solver3d_init
end interface

interface execute
   module procedure poisfft_solver1d__execute
   module procedure poisfft_solver2d__execute
   module procedure poisfft_solver3d__execute
end interface

contains

function poisfft_solver3d__new(nxyz, lxyz, bcs, approximation, gnxyz, offs, mpi_comm, nthreads) result(self)
   type(poisfft_solver3d) :: self

   integer, intent(in) :: nxyz(3)
   real(RP), intent(in) :: lxyz(3)
   integer, intent(in) :: bcs(6)
   integer, intent(in), optional :: approximation
   integer, intent(in), optional :: gnxyz(3)
   integer, intent(in), optional :: offs(3)
   integer, intent(in), optional :: mpi_comm
   integer, intent(in), optional :: nthreads

   self % nxyz = nxyz

   self % lx = lxyz(1)
   self % ly = lxyz(2)
   self % lz = lxyz(3)

   self % nx = nxyz(1)
   self % ny = nxyz(2)
   self % nz = nxyz(3)

   if (present(gnxyz)) then
      self % gnx = gnxyz(1)
      self % gny = gnxyz(2)
      self % gnz = gnxyz(3)
   else
      self % gnx = self % nx
      self % gny = self % ny
      self % gnz = self % nz
   end if

   if (present(offs)) then
      self % offx = offs(1)
      self % offy = offs(2)
      self % offz = offs(3)
   end if

   self % cnt = product(self % nxyz)
   self % gcnt = product(int([self % gnx, self % gny, self % gnz], kind(self % gcnt)))
   self % bcs = bcs

   if (self % bcs(1) /= self % bcs(2) .or. self % bcs(3) /= self % bcs(4) .or. self % bcs(5) /= self % bcs(6)) then
      stop 'both boundary consitions in one direction must be identical.'
   end if

   if (present(approximation)) self % approximation = approximation

   if (present(mpi_comm)) then
      self % mpi % comm = mpi_comm
#ifdef MPI
   else
      stop 'No PFFT comm present in PoisFFT_Solver3D__New.'
#endif
   end if

   self % nthreads = merge(nthreads, 1, present(nthreads))

   !create fftw plans and allocate working arrays
   call init(self)
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

   if (all(self % bcs == poisfft_periodic)) then

#ifdef MPI
      if (self % nthreads > 1) call poisfft_pfft_initthreads(self % nthreads)
#else
      if (self % nthreads > 1) call poisfft_initthreads(self % nthreads)
#endif

      call allocate_fftw_complex(self)

      self % forward = poisfft_plan3d(self, [fft_complex, fftw_forward])
      self % backward = poisfft_plan3d(self, [fft_complex, fftw_backward])

   else if ( &
      all(self % bcs == poisfft_dirichlet) .or. &
      all(self % bcs == poisfft_dirichletstag) .or. &
      all(self % bcs == poisfft_neumann) .or. &
      all(self % bcs == poisfft_neumannstag)) then

#ifdef MPI
      if (self % nthreads > 1) call poisfft_pfft_initthreads(self % nthreads)
#else
      if (self % nthreads > 1) call poisfft_initthreads(self % nthreads)
#endif

      call allocate_fftw_real(self)

      real_forw = real_transform_type_forward(self % bcs(1))
      real_back = real_transform_type_backward(self % bcs(1))

      self % forward = poisfft_plan3d(self, [(real_forw, i=1, 3)])
      self % backward = poisfft_plan3d(self, [(real_back, i=1, 3)])

   else if ( &
      all(self % bcs(1:4) == poisfft_periodic) .and. &
      (all(self % bcs(5:6) == poisfft_neumann) .or. all(self % bcs(5:6) == poisfft_neumannstag))) then

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

#ifdef MPI
      !1..x, 2..y, 3..z, not different threads, but directions
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
       case(PoisFFT_Dirichlet)
         res = FFT_RealOdd00
       case(PoisFFT_DirichletStag)
         res = FFT_RealOdd10
       case(PoisFFT_Neumann)
         res = FFT_RealEven00
       case(PoisFFT_NeumannStag)
         res = FFT_RealEven10
       case default
         res = -1
      end select
   end function

   function real_transform_type_backward(bc) result(res)
      integer :: res
      integer, intent(in) :: bc
      select case(bc)
       case(PoisFFT_Dirichlet)
         res = FFT_RealOdd00
       case(PoisFFT_DirichletStag)
         res = FFT_RealOdd01
       case(PoisFFT_Neumann)
         res = FFT_RealEven00
       case(PoisFFT_NeumannStag)
         res = FFT_RealEven01
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
         phi(ngphi(1) + 1:ngphi(1) + self % nx, &
         ngphi(2) + 1:ngphi(2) + self % ny, &
         ngphi(3) + 1:ngphi(3) + self % nz), &
         rhs(ngrhs(1) + 1:ngrhs(1) + self % nx, &
         ngrhs(2) + 1:ngrhs(2) + self % ny, &
         ngrhs(3) + 1:ngrhs(3) + self % nz))

   else if (all(self % bcs == poisfft_dirichlet) .or. all(self % bcs == poisfft_dirichletstag)) then
      call poisfft_solver3d_fulldirichlet(self, &
         phi(ngphi(1) + 1:ngphi(1) + self % nx, &
         ngphi(2) + 1:ngphi(2) + self % ny, &
         ngphi(3) + 1:ngphi(3) + self % nz), &
         rhs(ngrhs(1) + 1:ngrhs(1) + self % nx, &
         ngrhs(2) + 1:ngrhs(2) + self % ny, &
         ngrhs(3) + 1:ngrhs(3) + self % nz))

   else if (all(self % bcs == poisfft_neumann) .or. all(self % bcs == poisfft_neumannstag)) then
      call poisfft_solver3d_fullneumann(self, &
         phi(ngphi(1) + 1:ngphi(1) + self % nx, &
         ngphi(2) + 1:ngphi(2) + self % ny, &
         ngphi(3) + 1:ngphi(3) + self % nz), &
         rhs(ngrhs(1) + 1:ngrhs(1) + self % nx, &
         ngrhs(2) + 1:ngrhs(2) + self % ny, &
         ngrhs(3) + 1:ngrhs(3) + self % nz))

   else if (all(self % bcs(1:4) == poisfft_periodic) .and. (all(self % bcs(5:6) == poisfft_neumann) .or. all(self % bcs(5:6) == poisfft_neumannstag))) then
      call poisfft_solver3d_ppns(self, &
         phi(ngphi(1) + 1:ngphi(1) + self % nx, &
         ngphi(2) + 1:ngphi(2) + self % ny, &
         ngphi(3) + 1:ngphi(3) + self % nz), &
         rhs(ngrhs(1) + 1:ngrhs(1) + self % nx, &
         ngrhs(2) + 1:ngrhs(2) + self % ny, &
         ngrhs(3) + 1:ngrhs(3) + self % nz))

   else if (all(self % bcs(1:2) == poisfft_periodic) .and. (all(self % bcs(3:6) == poisfft_neumann) .or. all(self % bcs(3:6) == poisfft_neumannstag))) then
      call poisfft_solver3d_pnsns(self, &
         phi(ngphi(1) + 1:ngphi(1) + self % nx, &
         ngphi(2) + 1:ngphi(2) + self % ny, &
         ngphi(3) + 1:ngphi(3) + self % nz), &
         rhs(ngrhs(1) + 1:ngrhs(1) + self % nx, &
         ngrhs(2) + 1:ngrhs(2) + self % ny, &
         ngrhs(3) + 1:ngrhs(3) + self % nz))
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
      !we see the dimensions reversed in fortran!
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
      !we see the dimensions reversed in fortran!
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
   integer, intent(in), optional :: approximation
   integer, intent(in), optional :: gnxyz(2)
   integer, intent(in), optional :: offs(2)
   integer, intent(in), optional :: mpi_comm
   integer, intent(in), optional :: nthreads

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
      stop 'No PFFT comm present in PoisFFT_Solver2D__New.'
#endif
   end if

   self % nthreads = merge(nthreads, 1, present(nthreads))

   !create fftw plans and allocate working array
   call init(self)
end function

subroutine poisfft_solver2d_init(self)
   type(poisfft_solver2d), intent(inout) :: self
   integer :: i

   self % norm_factor = norm_factor(self % gnx, self % bcs(1:2)) * &
      norm_factor(self % gny, self % bcs(3:4))

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
         phi(ngphi(1) + 1:ngphi(1) + self % nx, &
         ngphi(2) + 1:ngphi(2) + self % ny), &
         rhs(ngrhs(1) + 1:ngrhs(1) + self % nx, &
         ngrhs(2) + 1:ngrhs(2) + self % ny))
   else if (all(self % bcs == poisfft_dirichlet) .or. &
      all(self % bcs == poisfft_dirichletstag)) then
      call poisfft_solver2d_fulldirichlet(self, &
         phi(ngphi(1) + 1:ngphi(1) + self % nx, &
         ngphi(2) + 1:ngphi(2) + self % ny), &
         rhs(ngrhs(1) + 1:ngrhs(1) + self % nx, &
         ngrhs(2) + 1:ngrhs(2) + self % ny))
   else if (all(self % bcs == poisfft_neumann) .or. &
      all(self % bcs == poisfft_neumannstag)) then
      call poisfft_solver2d_fullneumann(self, &
         phi(ngphi(1) + 1:ngphi(1) + self % nx, &
         ngphi(2) + 1:ngphi(2) + self % ny), &
         rhs(ngrhs(1) + 1:ngrhs(1) + self % nx, &
         ngrhs(2) + 1:ngrhs(2) + self % ny))
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
   integer, intent(in), optional :: approximation
   integer, intent(in), optional :: gnxyz(1)
   integer, intent(in), optional :: offs(1)
   integer, intent(in), optional :: mpi_comm
   integer, intent(in), optional :: nthreads
#ifdef MPI
   integer :: ie, dims
#endif
   self % nxyz = nxyz
   self % lx = lxyz(1)
   self % nx = nxyz(1)

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
      if (ie /= 0) stop 'Error executing MPI_Cartdim_get.'
      self % mpi_transpose_needed = dims > 0
   else
      stop 'No PFFT comm present in PoisFFT_Solver1D__New.'
#endif
   end if

   self % nthreads = merge(nthreads, 1, present(nthreads))

   !create fftw plans and allocate working array
   call init(self)
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

#include "poisfft-solvers-inc.f90"

function eigenvalues(f, bcs, l, n, gn, off) result(res)
   procedure(eig_fn_spectral) :: f
   integer, intent(in) :: bcs(2)
   real(RP), intent(in) :: l
   integer, intent(in) :: n, gn, off
   real(RP) :: res(n)
   integer :: i
   real(RP) :: dkx, dkx_h

   dkx = pi * grid_dx(gn, l, bcs) / l
   dkx_h = dkx / 2

   if (all(bcs == poisfft_periodic)) then
      do i = 1, n
         if (i + off < gn / 2) then
            res(i) = f((i - 1 + off) * dkx)
         else
            res(i) = f((gn - (i - 1 + off)) * dkx)
         end if
      end do
   else if (all(bcs == poisfft_dirichlet)) then
      forall(i=1:n) res(i) = f((i + off) * dkx_h)
   else if (all(bcs == poisfft_dirichletstag)) then
      forall(i=1:n) res(i) = f((i + off) * dkx_h)
   else if (all(bcs == poisfft_neumann)) then
      forall(i=1:n) res(i) = f((i - 1 + off) * dkx_h)
   else if (all(bcs == poisfft_neumannstag)) then
      forall(i=1:n) res(i) = f((i - 1 + off) * dkx_h)
   end if

   res = res / (grid_dx(gn, l, bcs))**2
end function


pure real(RP) function eig_fn_spectral(x) result(f)
   real(RP), intent(in) :: x
   f = - (2 * x)**2
end function

pure real(RP) function eig_fn_FD2(x) result(f)
   real(RP), intent(in) :: x
   f = - (2 * sin(x))**2
end function

pure real(RP) function eig_fn_FD4(x) result(f)
   real(RP), intent(in) :: x
   f = - ((sin(3 * x) - 27 * sin(x)) / 12)**2
end function

function grid_dx(gn, l, bcs) result(res)
   real(RP) :: res
   integer, intent(in) :: gn
   real(RP), intent(in) :: l
   integer, intent(in) :: bcs(2)

   if (all(bcs == poisfft_periodic)) then
      res = l / gn
   else if (all(bcs == poisfft_dirichlet)) then
      res = l / (gn + 1)
   else if (all(bcs == poisfft_dirichletstag)) then
      res = l / gn
   else if (all(bcs == poisfft_neumann)) then
      res = l / (gn - 1)
   else if (all(bcs == poisfft_neumannstag)) then
      res = l / gn
   else
      res = 0
   end if
end function

function norm_factor(gn, bcs) result(res)
   real(RP) :: res
   integer, intent(in) :: gn
   integer, intent(in) :: bcs(2)

   if (all(bcs == poisfft_periodic)) then
      res = gn
   else if (all(bcs == poisfft_dirichlet)) then
      res = 2 * (gn + 1)
   else if (all(bcs == poisfft_dirichletstag)) then
      res = 2 * gn
   else if (all(bcs == poisfft_neumann)) then
      res = 2 * (gn - 1)
   else if (all(bcs == poisfft_neumannstag)) then
      res = 2 * gn
   else
      res = 0
   end if
end function

#ifdef MPI
subroutine init_mpi_buffers(self, dir)
   interface
      subroutine mpi_alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, &
         recvtype, comm, ierror)
         integer sendbuf(*), recvbuf(*)
         integer sendcount, sendtype, recvcount, recvtype
         integer comm, ierror
      end subroutine
   end interface

   type(poisfft_solver3d), intent(inout), target :: self
   integer, intent(in) :: dir
   type(mpi_vars_1d), pointer :: mpi
   integer :: i, ie

   mpi => self % Solvers1D(dir) % mpi

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
      do i = 2, mpi % np
         mpi % rdispls(i) = mpi % rdispls(i - 1) + mpi % rcounts(i - 1)
      end do

      do i = 1, mpi % np
         mpi % sumrnzs(i) = sum(mpi % rnzs(1:i - 1))
      end do

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
      stop 'Not implemented.'
   end if
end subroutine

#endif

#undef RP
#undef CP
#undef _RP
#undef _CP
#undef MPI_RP
