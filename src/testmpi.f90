module endianness
   use iso_fortran_env
   implicit none

   private
   public :: getendianness, bigend

   logical, save :: littleendian = .true.

   interface bigend
      module procedure bigend32
      module procedure bigend64
   end interface

contains

   subroutine getendianness
      integer(int8), dimension(4) :: bytes !may not work on some processors

      bytes = transfer(1_int32, bytes, 4)
      if (bytes(4) == 1) then
         littleendian = .false.
      else
         littleendian = .true.
      endif
   end subroutine

   elemental function bigend32(x) result(res)
      real(real32) :: res
      real(real32), intent(in) :: x
      integer(int8), dimension(4) :: bytes !may not work on some processors

      if (.not. littleendian) then
         res = x
      else
         bytes = transfer(x, bytes, 4)
         res = transfer(bytes(4:1:-1), res)
      endif
   end function

   elemental function bigend64(x) result(res)
      real(real64) :: res
      real(real64), intent(in) :: x
      integer(int8), dimension(8) :: bytes !may not work on some processors

      if (.not. littleendian) then
         res = x
      else
         bytes = transfer(x, bytes, 8)
         res = transfer(bytes(8:1:-1), res)
      endif
   end function

end module

module parameters
   use poisfft_constants
   use iso_c_binding
   integer, parameter :: rp = drp
   real(rp), parameter :: pi = 4 * atan(1._rp) ! 3.141592653589793238462_rp
end module

module my_mpi
   use parameters
   use mpi

   integer :: nimg, npxyz(3), pxyz(3)
   integer :: nxims, nyims, nzims
   integer :: myim, myrank, iim, jim, kim
   integer :: w_im, e_im, s_im, n_im, b_im, t_im
   integer :: w_rk, e_rk, s_rk, n_rk, b_rk, t_rk
   integer :: glob_comm, cart_comm, cart_comm_dim = -1
   logical :: rt = .true.
   integer :: mpi_rp = -huge(1)

   interface error_stop
      module procedure error_stop_int
      module procedure error_stop_char
   end interface

contains
   integer function caf_this_image() result(res)
      integer :: ie
      call mpi_comm_rank(glob_comm, res, ie)
      res = res + 1
      if (ie /= 0) call error_stop('mpi_comm_rank error')
   end function

   integer function caf_num_images() result(res)
      integer :: ie
      call mpi_comm_size(glob_comm, res, ie)
      if (ie /= 0) call error_stop('mpi_comm_size error')
   end function

   integer function caf_image_index(sub) result(res)
      integer, intent(in) :: sub(3)
      integer :: ie

      if (cart_comm_dim == -1) then
         call mpi_cartdim_get(cart_comm, cart_comm_dim, ie)
         if (ie /= 0) call error_stop('mpi_cartdim_get')
      end if
      call mpi_cart_rank(cart_comm, sub(3:4 - cart_comm_dim:-1) - 1, res, ie)
      if (ie /= 0) call error_stop('mpi_cart_rank')
      res = res + 1
   end function

   subroutine get_image_coords()
      integer :: ie
      call mpi_cart_coords(cart_comm, myrank, 3, pxyz, ie)
      if (ie /= 0) call error_stop('mpi_cart_coords')

      pxyz = pxyz(3:1:-1)

      iim = pxyz(1) + 1
      jim = pxyz(2) + 1
      kim = pxyz(3) + 1

      w_im = caf_image_index([merge(iim - 1, nxims, iim > 1), jim, kim])
      e_im = caf_image_index([merge(iim + 1, 1, iim < nxims), jim, kim])
      s_im = caf_image_index([iim, merge(jim - 1, nyims, jim > 1), kim])
      n_im = caf_image_index([iim, merge(jim + 1, 1, jim < nyims), kim])
      b_im = caf_image_index([iim, jim, merge(kim - 1, nzims, kim > 1)])
      t_im = caf_image_index([iim, jim, merge(kim + 1, 1, kim < nzims)])

      w_rk = w_im - 1; e_rk = e_im - 1
      s_rk = s_im - 1; n_rk = n_im - 1
      b_rk = b_im - 1; t_rk = t_im - 1
   end subroutine

   subroutine error_stop_int(n)
      integer, intent(in) :: n
      integer :: ie
      if (rt) write(ferr, *) 'error', n, ' ie=', ie
      call mpi_finalize(ie)
      error stop 11
   end subroutine

   subroutine error_stop_char(ch)
      character(*), intent(in) :: ch
      integer :: ie
      if (rt) write(ferr, *) 'error', ch, ' ie=', ie
      call mpi_finalize(ie)
      error stop 12
   end subroutine
end module

module subs
   use parameters
   use my_mpi
   use poisfft
   implicit none

contains
   subroutine exchange_boundaries_3d(comm, phi, nx, ny, nz, bcs)
      integer, intent(in) :: comm
      real(rp), intent(inout), contiguous :: phi(0:, 0:, 0:)
      integer, intent(in) :: nx, ny, nz, bcs(6)
      logical :: oddx, oddy, oddz, evenx, eveny, evenz
      integer :: ierr, tag, status(mpi_status_size)

      oddx = mod(iim, 2) == 1
      evenx = .not. oddx
      oddy = mod(jim, 2) == 1
      eveny = .not. oddy
      oddz = mod(kim, 2) == 1
      evenz = .not. oddz
      call mpi_barrier(comm, ierr)

#define SEND_W3 if (iim > 1) then; call send(phi(1, 1:ny, 1:nz), w_rk); end if
#define RECV_W3 if (iim > 1) then; call recv(phi(0, 1:ny, 1:nz), w_rk); end if
#define SEND_E3 if (iim < nxims) then; call send(phi(nx, 1:ny, 1:nz), e_rk); end if
#define RECV_E3 if (iim < nxims) then; call recv(phi(nx + 1, 1:ny, 1:nz), e_rk); end if
#define SEND_S3 if (jim > 1) then; call send(phi(1:nx, 1, 1:nz), s_rk); end if
#define RECV_S3 if (jim > 1) then; call recv(phi(1:nx, 0, 1:nz), s_rk); end if
#define SEND_N3 if (jim < nyims) then; call send(phi(1:nx, ny, 1:nz), n_rk); end if
#define RECV_N3 if (jim < nyims) then; call recv(phi(1:nx, ny + 1, 1:nz), n_rk); end if
#define SEND_B3 if (kim > 1) then; call send(phi(1:nx, 1:ny, 1), b_rk); end if
#define RECV_B3 if (kim > 1) then; call recv(phi(1:nx, 1:ny, 0), b_rk); end if
#define SEND_T3 if (kim < nzims) then; call send(phi(1:nx, 1:ny, nz), t_rk); end if
#define RECV_T3 if (kim < nzims) then; call recv(phi(1:nx, 1:ny, nz + 1), t_rk); end if

      ! internal boundaries
      if (oddx) then; SEND_W3; else; RECV_E3; end if
      if (evenx) then; SEND_W3; else; RECV_E3; end if

      if (oddx) then; SEND_E3; else; RECV_W3; end if
      if (evenx) then; SEND_E3; else; RECV_W3; end if

      if (oddy) then; SEND_S3; else; RECV_N3; end if
      if (eveny) then; SEND_S3; else; RECV_N3; end if

      if (oddy) then; SEND_N3; else; RECV_S3; end if
      if (eveny) then; SEND_N3; else; RECV_S3; end if

      if (oddz) then; SEND_B3; else; RECV_T3; end if
      if (evenz) then; SEND_B3; else; RECV_T3; end if

      if (oddz) then; SEND_T3; else; RECV_B3; end if
      if (evenz) then; SEND_T3; else; RECV_B3; end if

      ! global domain boundaries
      if (bcs(1) == poisfft_periodic) then
         if (nxims > 1) then
            if (iim == 1) then
               call send(phi(1, 1:ny, 1:nz), w_rk)
            else if (iim == nxims) then
               call recv(phi(nx + 1, 1:ny, 1:nz), e_rk)
            end if
            if (iim == nxims) then
               call send(phi(nx, 1:ny, 1:nz), e_rk)
            else if (iim == 1) then
               call recv(phi(0, 1:ny, 1:nz), w_rk)
            end if
         else
            phi(0, :, :) = phi(nx, :, :)
            phi(nx + 1, :, :) = phi(1, :, :)
         end if
      else
         if (bcs(1) == poisfft_neumannstag) then
            if (iim == 1) phi(0, :, :) = phi(1, :, :)
         elseif(bcs(1) == poisfft_dirichletstag) then
            if (iim == 1) phi(0, :, :) = -phi(1, :, :)
         end if
         if (bcs(2) == poisfft_neumannstag) then
            if (iim == nxims) phi(nx + 1, :, :) = phi(nx, :, :)
         elseif(bcs(2) == poisfft_dirichletstag) then
            if (iim == nxims) phi(nx + 1, :, :) = -phi(nx, :, :)
         end if
      end if

      if (bcs(3) == poisfft_periodic) then
         if (nyims > 1) then
            if (jim == 1) then
               call send(phi(1:nx, 1, 1:nz), s_rk)
            else if (jim == nyims) then
               call recv(phi(1:nx, ny + 1, 1:nz), n_rk)
            end if
            if (jim == nyims) then
               call send(phi(1:nx, ny, 1:nz), n_rk)
            else if (jim == 1) then
               call recv(phi(1:nx, 0, 1:nz), s_rk)
            end if
         else
            phi(:, 0, :) = phi(:, ny, :)
            phi(:, ny + 1, :) = phi(:, 1, :)
         end if
      else
         if (bcs(3) == poisfft_neumannstag) then
            if (jim == 1) phi(:, 0, :) = phi(:, 1, :)
         elseif(bcs(3) == poisfft_dirichletstag) then
            if (jim == 1) phi(:, 0, :) = -phi(:, 1, :)
         end if
         if (bcs(4) == poisfft_neumannstag) then
            if (jim == nyims) phi(:, ny + 1, :) = phi(:, ny, :)
         elseif(bcs(4) == poisfft_dirichletstag) then
            if (jim == nyims) phi(:, ny + 1, :) = -phi(:, ny, :)
         end if
      end if

      if (bcs(5) == poisfft_periodic) then
         if (nzims > 1) then
            if (kim == 1) then
               call send(phi(1:nx, 1:ny, 1), b_rk)
            else if (kim == nzims) then
               call recv(phi(1:nx, 1:ny, nz + 1), t_rk)
            end if
            if (kim == nzims) then
               call send(phi(1:nx, 1:ny, nz), t_rk)
            else if (kim == 1) then
               call recv(phi(1:nx, 1:ny, 0), b_rk)
            end if
         else
            phi(:, :, 0) = phi(:, :, nz)
            phi(:, :, nz + 1) = phi(:, :, 1)
         end if
      else
         if (bcs(5) == poisfft_neumannstag) then
            if (kim == 1) phi(:, :, 0) = phi(:, :, 1)
         elseif(bcs(5) == poisfft_dirichletstag) then
            if (kim == 1) phi(:, :, 0) = -phi(:, :, 1)
         end if
         if (bcs(6) == poisfft_neumannstag) then
            if (kim == nzims) phi(:, :, nz + 1) = phi(:, :, nz)
         elseif(bcs(6) == poisfft_dirichletstag) then
            if (kim == nzims) phi(:, :, nz + 1) = -phi(:, :, nz)
         end if
      end if

      call mpi_barrier(comm, ierr)
   contains
      subroutine send(a, to)
         real(rp), intent(in) :: a(:, :)
         integer, intent(in) :: to

         call mpi_send(a, size(a), mpi_rp, to, 1, comm, ierr)
         if (ierr /= 0) stop 'error sending mpi message.'
      end subroutine

      subroutine recv(a, from)
         real(rp), intent(out) :: a(:, :)
         integer, intent(in) :: from

         call mpi_recv(a, size(a), mpi_rp, from, 1, comm, status, ierr)
         if (ierr /= 0) stop 'error sending mpi message.'
      end subroutine
   end subroutine

   subroutine exchange_boundaries_2d(comm, phi, nx, ny, bcs)
      integer, intent(in) :: comm
      real(rp), intent(inout), contiguous :: phi(0:, 0:)
      integer, intent(in) :: nx, ny, bcs(4)
      logical :: oddx, oddy, evenx, eveny
      integer :: ierr, tag, status(mpi_status_size)

      oddx = mod(iim, 2) == 1
      evenx = .not. oddx
      oddy = mod(jim, 2) == 1
      eveny = .not. oddy
      call mpi_barrier(comm, ierr)

#define SEND_W2 if (iim > 1) then; call send(phi(1, 1:ny), w_rk); end if
#define RECV_W2 if (iim > 1) then; call recv(phi(0, 1:ny), w_rk); end if
#define SEND_E2 if (iim < nxims) then; call send(phi(nx, 1:ny), e_rk); end if
#define RECV_E2 if (iim < nxims) then; call recv(phi(nx + 1, 1:ny), e_rk); end if
#define SEND_S2 if (jim > 1) then; call send(phi(1:nx, 1), s_rk); end if
#define RECV_S2 if (jim > 1) then; call recv(phi(1:nx, 0), s_rk); end if
#define SEND_N2 if (jim < nyims) then; call send(phi(1:nx, ny), n_rk); end if
#define RECV_N2 if (jim < nyims) then; call recv(phi(1:nx, ny + 1), n_rk); end if

      ! internal boundaries
      if (oddx) then; SEND_W2; else; RECV_E2; end if
      if (evenx) then; SEND_W2; else; RECV_E2; end if

      if (oddx) then; SEND_E2; else; RECV_W2; end if
      if (evenx) then; SEND_E2; else; RECV_W2; end if

      if (oddy) then; SEND_S2; else; RECV_N2; end if
      if (eveny) then; SEND_S2; else; RECV_N2; end if

      if (oddy) then; SEND_N2; else; RECV_S2; end if
      if (eveny) then; SEND_N2; else; RECV_S2; end if

      !global domain boundaries
      if (bcs(1) == poisfft_periodic) then
         if (nxims > 1) then
            if (iim == 1) then
               call send(phi(1, 1:ny), w_rk)
            else if (iim == nxims) then
               call recv(phi(nx + 1, 1:ny), e_rk)
            end if
            if (iim == nxims) then
               call send(phi(nx, 1:ny), e_rk)
            else if (iim == 1) then
               call recv(phi(0, 1:ny), w_rk)
            end if
         else
            phi(0, :) = phi(nx, :)
            phi(nx + 1, :) = phi(1, :)
         end if
      else
         if (bcs(1) == poisfft_neumannstag) then
            if (iim == 1) phi(0, :) = phi(1, :)
         elseif(bcs(1) == poisfft_dirichletstag) then
            if (iim == 1) phi(0, :) = -phi(1, :)
         end if
         if (bcs(2) == poisfft_neumannstag) then
            if (iim == nxims) phi(nx + 1, :) = phi(nx, :)
         elseif(bcs(2) == poisfft_dirichletstag) then
            if (iim == nxims) phi(nx + 1, :) = -phi(nx, :)
         end if
      end if

      if (bcs(3) == poisfft_periodic) then
         if (nyims > 1) then
            if (jim == 1) then
               call send(phi(1:nx, 1), s_rk)
            else if (jim == nyims) then
               call recv(phi(1:nx, ny + 1), n_rk)
            end if
            if (jim == nyims) then
               call send(phi(1:nx, ny), n_rk)
            else if (jim == 1) then
               call recv(phi(1:nx, 0), s_rk)
            end if
         else
            phi(:, 0) = phi(:, ny)
            phi(:, ny + 1) = phi(:, 1)
         end if
      else
         if (bcs(3) == poisfft_neumannstag) then
            if (jim == 1) phi(:, 0) = phi(:, 1)
         elseif(bcs(3) == poisfft_dirichletstag) then
            if (jim == 1) phi(:, 0) = -phi(:, 1)
         end if
         if (bcs(4) == poisfft_neumannstag) then
            if (jim == nyims) phi(:, ny + 1) = phi(:, ny)
         elseif(bcs(4) == poisfft_dirichletstag) then
            if (jim == nyims) phi(:, ny + 1) = -phi(:, ny)
         end if
      end if

      call mpi_barrier(comm, ierr)

   contains
      subroutine send(a, to)
         real(rp), intent(in) :: a(:)
         integer, intent(in) :: to

         call mpi_send(a, size(a), mpi_rp, to, 1, comm, ierr)
         if (ierr /= 0) stop 'error sending mpi message.'
      end subroutine

      subroutine recv(a, from)
         real(rp), intent(out) :: a(:)
         integer, intent(in) :: from

         call mpi_recv(a, size(a), mpi_rp, from, 1, comm, status, ierr)
         if (ierr /= 0) stop 'error sending mpi message.'
      end subroutine
   end subroutine

   subroutine res3d(nx, ny, nz, phi, rhs, aw, ae, as, an, ab, at, r)
      implicit none
      intrinsic mod, abs, max
      integer, parameter :: knd = rp

      integer, intent(in) :: nx, ny, nz
      real(knd), dimension(0:nx + 1, 0:ny + 1, 0:nz + 1), intent(inout) :: phi
      real(knd), dimension(1:nx, 1:ny, 1:nz), intent(in) :: rhs
      real(knd), intent(in) :: aw, ae, as, an, ab, at
      real(knd), intent(out) :: r
      integer i, j, k, l
      real(knd) :: p, ap
      integer :: ie

      r = 0
      ap = aw + ae + as + an + ab + at
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               p = 0
               p = p + phi(i - 1, j, k) * aw
               p = p + phi(i + 1, j, k) * ae
               p = p + phi(i, j - 1, k) * as
               p = p + phi(i, j + 1, k) * an
               p = p + phi(i, j, k - 1) * ab
               p = p + phi(i, j, k + 1) * at
               p = p - rhs(i, j, k)
               p = abs(-p + ap * phi(i, j, k))
               r = max(r, abs(p))
            end do
         end do
      end do

      call mpi_allreduce(mpi_in_place, r, 1, mpi_rp, mpi_max, glob_comm, ie)
   end subroutine

   subroutine res2d(nx, ny, phi, rhs, aw, ae, as, an, r)
      implicit none
      intrinsic mod, abs, max
      integer, parameter :: knd = rp

      integer, intent(in) :: nx, ny
      real(knd), dimension(0:nx + 1, 0:ny + 1), intent(inout) :: phi
      real(knd), dimension(1:nx, 1:ny), intent(in) :: rhs
      real(knd), intent(in) :: aw, ae, as, an
      real(knd), intent(out) :: r
      integer i, j, k, l
      real(knd) :: p, ap
      integer :: ie

      r = 0
      ap = aw + ae + as + an

      do j = 1, ny
         do i = 1, nx
            p = 0
            p = p + phi(i - 1, j) * aw
            p = p + phi(i + 1, j) * ae
            p = p + phi(i, j - 1) * as
            p = p + phi(i, j + 1) * an
            p = p - rhs(i, j)
            p = abs(-p + ap * phi(i, j))
            r = max(r, abs(p))
         end do
      end do

      call mpi_allreduce(mpi_in_place, r, 1, mpi_rp, mpi_max, glob_comm, ie)
   end subroutine
end module subs

program testpoisson_mpi
   use iso_fortran_env
   use iso_c_binding
   use poisfft, &
      poisfft_solver1d => poisfft_solver1d_dp, &
      poisfft_solver2d => poisfft_solver2d_dp, &
      poisfft_solver3d => poisfft_solver3d_dp
   use my_mpi
   use subs
   implicit none

   integer(c_intptr_t) :: ng(3) = [21, 32, 43] ! [131, 123, 127]
   real(rp), dimension(:, :, :), allocatable :: phi, rhs
   real(rp) :: dx, dy, dz, ls(3)
   integer i, j, k
   real(rp) :: r, x, y, z
   integer(int64) :: t1, t2, trate
   type(poisfft_solver3d) :: solver3d
   type(poisfft_solver2d) :: solver2d
   type(poisfft_solver1d) :: solver1d
   character(50) :: ch50
   integer :: nx, ny, nz
   integer(c_intptr_t), dimension(3) :: nxyz, off, nxyz2, nsxyz2
   real(drp) :: s
   integer :: ie = 0
   character(200) :: fname
   logical :: flag

   integer :: required, provided

   ! MPI_THREAD_SINGLE=0 < MPI_THREAD_FUNNELED=1 < MPI_THREAD_SERIALIZED=2 < MPI_THREAD_MULTIPLE=3
   required = mpi_thread_serialized

   call mpi_initialized(flag, ie)
   if (.not. flag) then
      call mpi_init_thread(required, provided, ie)
      if (ie /= 0) call error_stop('Error mpi_init_thread.')
   else
      call mpi_query_thread(provided, ie)
      if (ie /= 0) call error_stop('Error mpi_query_thread.')
   end if

   glob_comm = mpi_comm_world
   myim = caf_this_image()
   myrank = myim - 1

   if (myrank /= 0) rt = .false.

   if (provided < required) then
      if (rt) then
         write(ferr, *) '------------------------------'
         write(ferr, *) 'error, the provided mpi threading support smaller than required!'
         write(ferr, *) 'required:', required
         write(ferr, *) 'provided:', provided
         write(ferr, *) 'trying to continue anyway, but a crash is likely and the results will be questionable.'
         write(ferr, *) '------------------------------'
      end if
   end if

   if (rp == kind(1.)) then
      mpi_rp = mpi_real
   else if (rp == kind(dble(1.))) then
      mpi_rp = mpi_double_precision
   end if

   call system_clock(count_rate=trate)

   nimg = caf_num_images()

   if (command_argument_count() >= 1) then
      call get_command_argument(1, value=ch50)
      read(ch50, *, iostat=ie) npxyz
      if (ie /= 0) then
         write(ferr, *) 'the process grid should be provided as "npx, npy, npz" where nxp,npy and npz are integers.'
         error stop 111
      end if
   else
      npxyz(1) = 1
      npxyz(2) = nint(sqrt(real(nimg)))
      npxyz(3) = nimg / npxyz(2)
      if (rt) write(*, *) 'trying to decompose in', npxyz, 'process grid.'
   end if

   if (command_argument_count() >= 2) then
      call get_command_argument(2, value=ch50)
      read(ch50, *) ng
   end if

   nxims = npxyz(1)
   nyims = npxyz(2)
   nzims = npxyz(3)

   if (product(npxyz) /= nimg) then
      if (rt) then
         write(ferr, *) 'could not decompose the processes to n x n grid.'
         write(ferr, *) 'try a perfect square for number of processes.'
      end if
      call error_stop(25)
   end if

   call poisfft_initmpigrid(glob_comm, npxyz(3:2:-1), cart_comm, ie)
   if (ie /= 0) call error_stop(30)

   call get_image_coords

   call poisfft_localgridsize(3, ng, cart_comm, nxyz, off, nxyz2, nsxyz2)
   if (any(nxyz /= nxyz2) .or. any(off /= nsxyz2)) call error_stop(40)

   if (any(nxyz == 0)) then
      write(ferr, *) 'process', pxyz, 'has grid dimensions', nxyz, '.'
      write(ferr, *) 'try different process grid distribution.'
      call error_stop(45)
   end if

   nx = nxyz(1)
   ny = nxyz(2)
   nz = nxyz(3)

   ls = [2 * pi, 2 * pi, 2 * pi]
   dx = ls(1) / ng(1)
   dy = ls(2) / ng(2)
   dz = ls(3) / ng(3)

   allocate(rhs(nx, ny, nz), stat=ie)
   if (ie /= 0) call error_stop(50)

   allocate(phi(0:nx + 1, 0:ny + 1, 0:nz + 1))
   if (ie /= 0) call error_stop(60)

   call random_seed(put=[(0, i=1, 100)])
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            call random_number(rhs(i, j, k))
            call random_number(phi(i, j, k))
            ! x = (i + off(1) - 1._rp / 2) * dx
            ! y = (j + off(2) - 1._rp / 2) * dy
            ! z = (k + off(3) - 1._rp / 2) * dz
            ! rhs(i,j,k) = x * sin(y) / (abs(z)+0.1)
            ! phi(i,j,k) = x * sin(y) / (abs(z)+0.1)
         end do
      end do
   end do

   call mpi_allreduce(sum(rhs), s, 1, mpi_rp, mpi_sum, glob_comm, ie)
   rhs = rhs - s / product(ng)

   call mpi_barrier(glob_comm, ie)
   if (rt) write(ferr, *) '3d pnsns'
   call compute3d([(poisfft_periodic, i=1, 2), (poisfft_neumannstag, i=3, 6)])
   call mpi_barrier(glob_comm, ie)
   if (rt) write(ferr, *) '3d ppns'
   call compute3d([(poisfft_periodic, i=1, 4), (poisfft_neumannstag, i=5, 6)])
   call mpi_barrier(glob_comm, ie)
   if (rt) write(ferr, *) '3d staggered dirichlet'
   call compute3d([(poisfft_dirichletstag, i=1, 6)])
   call mpi_barrier(glob_comm, ie)
   if (rt) write(ferr, *) '3d staggered neumann:'
   call compute3d([(poisfft_neumannstag, i=1, 6)])
   call mpi_barrier(glob_comm, ie)
   if (rt) write(ferr, *) '3d periodic'
   call compute3d([(poisfft_periodic, i=1, 6)])
   call mpi_barrier(glob_comm, ie)
   if (rt) write(ferr, *) '2d periodic'
   call compute2d([(poisfft_periodic, i=1, 4)])
   call mpi_barrier(glob_comm, ie)
   if (rt) write(ferr, *) '1d periodic'
   call compute1d([(poisfft_periodic, i=1, 2)])
   call mpi_barrier(glob_comm, ie)
   call save_vtk
   call mpi_finalize(ie)

contains
   subroutine compute3d(bcs)
      integer, intent(in) :: bcs(6)

      phi(1:nx, 1:ny, 1:nz) = rhs
      call exchange_boundaries_3d(glob_comm, phi, nx, ny, nz, bcs)
      call res3d(nx, ny, nz, phi, rhs, dx**(-2), dx**(-2), dy**(-2), dy**(-2), dz**(-2), dz**(-2), r)

      if (rt) write(*, *) 'r1:', r
      solver3d = poisfft_solver3d([nx, ny, nz], ls, bcs, poisfft_finitedifference2, int(ng), int(off), cart_comm)

      do i = 1, 2
         call system_clock(count=t1)
         call execute(solver3d, phi, rhs)
         call system_clock(count=t2)
         if (rt) write(ferr, *) 'solver cpu time', real(t2 - t1) / real(trate)
      end do

      call finalize(solver3d)
      call exchange_boundaries_3d(glob_comm, phi, nx, ny, nz, bcs)
      call res3d(nx, ny, nz, phi, rhs, dx**(-2), dx**(-2), dy**(-2), dy**(-2), dz**(-2), dz**(-2), r)

      if (rt) write(*, *) 'r2:', r
      if (rt) write(ferr, *) '--------'
   end subroutine

   subroutine compute2d(bcs)
      integer, intent(in) :: bcs(4)
      integer :: sub_comm, ie
      integer :: n

      call mpi_cart_sub(cart_comm, [.false., .true.], sub_comm, ie)
      call mpi_allreduce(sum(rhs(:, :, 1)), s, 1, mpi_rp, mpi_sum, sub_comm, ie)
      rhs(:, :, 1) = rhs(:, :, 1) - s / product(ng(1:2))

      call exchange_boundaries_2d(glob_comm, phi(:, :, 1), nx, ny, bcs)
      call res2d(nx, ny, phi(:, :, 1), rhs(:, :, 1), dx**(-2), dx**(-2), dy**(-2), dy**(-2), r)

      if (rt) write(*, *) 'r1:', r
      solver2d = poisfft_solver2d([nx, ny], ls(1:2), bcs, poisfft_finitedifference2, int(ng(1:2)), int(off(1:2)), sub_comm)

      do i = 1, 1
         call system_clock(count=t1)
         call execute(solver2d, phi(:, :, i), rhs(:, :, i))
         call system_clock(count=t2)
         if (rt) write(ferr, *) 'solver cpu time', real(t2 - t1) / real(trate)
      end do

      call finalize(solver2d)
      call exchange_boundaries_2d(glob_comm, phi(:, :, 1), nx, ny, bcs)
      call res2d(nx, ny, phi(:, :, 1), rhs(:, :, 1), dx**(-2), dx**(-2), dy**(-2), dy**(-2), r)

      if (rt) write(*, *) 'r2:', r
      if (rt) write(ferr, *) '--------'
   end subroutine

   subroutine compute1d(bcs)
      integer, intent(in) :: bcs(2)
      integer :: sub_comm = -1, ie

      call mpi_cart_sub(cart_comm, [.false., .true.], sub_comm, ie)
      call mpi_allreduce(sum(rhs(1, :, 1)), s, 1, mpi_rp, mpi_sum, sub_comm, ie)
      rhs(1, :, 1) = rhs(1, :, 1) - s / ng(2)

      ! call exchange_boundaries_1d(glob_comm, phi(:,:,1), nx, bcs)
      ! call res1d(nx,ny,phi(:,1,1),rhs(:,1,1), dx**(-2),dx**(-2), r)

      if (rt) write(*, *) 'r1:', r
      solver1d = poisfft_solver1d([ny], ls(1:1), bcs, poisfft_finitedifference2, int(ng(2:2)), int(off(2:2)), sub_comm)

      do i = 1, 1
         call system_clock(count=t1)
         call execute(solver1d, phi(1, :, i), rhs(1, :, i))
         call system_clock(count=t2)
         if (rt) write(ferr, *) 'solver cpu time', real(t2 - t1) / real(trate)
      end do
      call finalize(solver1d)

      ! call exchange_boundaries_1d(glob_comm, phi(:,1,1), nx, ny, bcs)
      ! call res1d(nx,ny,phi(:,1,1),rhs(:,1,1), dx**(-2),dx**(-2), r)
      if (rt) write(*, *) 'r2:', r
      if (rt) write(ferr, *) '--------'
   end subroutine

   subroutine save_vtk
      use endianness
      integer filetype, fh, unit
      integer(mpi_offset_kind) pos
      real(rp), allocatable :: buffer(:, :, :), buf(:)
      character :: lf = achar(10)
      character(70) :: str
      character(10) :: fm = '(*(1x,g0))'

      call getendianness

      if (rt) then
         open(newunit=unit, file='out.vtk', access='stream', status='replace', form='unformatted', action='write')

         write(unit) '# vtk DataFile Version 2.0', lf
         write(unit) 'CLMM output file', lf
         write(unit) 'BINARY', lf
         write(unit) 'DATASET RECTILINEAR_GRID', lf
         str = 'DIMENSIONS'
         write(str(12:), fm) ng(1), ng(2), ng(3)
         write(unit) str, lf
         str = 'X_COORDINATES'
         write(str(15:), fm) ng(1), 'double'
         write(unit) str, lf
         write(unit) bigend([((i - 0.5) * dx, i=1, ng(1))]), lf
         str = 'Y_COORDINATES'
         write(str(15:), fm) ng(2), 'double'
         write(unit) str, lf
         write(unit) bigend([((j - 0.5) * dy, j=1, ng(2))]), lf
         str = 'Z_COORDINATES'
         write(str(15:), fm) ng(3), 'double'
         write(unit) str, lf
         write(unit) bigend([((k - 0.5) * dz, k=1, ng(3))]), lf
         str = 'POINT_DATA'
         write(str(12:), fm) product(ng)
         write(unit) str, lf
         write(unit) lf
         write(unit) 'SCALARS Phi double', lf
         write(unit) 'LOOKUP_TABLE default', lf
         close(unit)
      end if

      call mpi_barrier(glob_comm, ie)

      call mpi_file_open(glob_comm, 'out.vtk', mpi_mode_append + mpi_mode_wronly, mpi_info_null, fh, ie)
      if (ie /= 0) call error_stop('open')

      call mpi_type_create_subarray(3, int(ng), int(nxyz), int(off), mpi_order_fortran, mpi_rp, filetype, ie)
      if (ie /= 0) call error_stop('create_subarray')

      call mpi_type_commit(filetype, ie)
      if (ie /= 0) call error_stop('type_commit')

      call mpi_barrier(glob_comm, ie)
      call mpi_file_get_position(fh, pos, ie)
      call mpi_barrier(glob_comm, ie)

      call mpi_file_set_view(fh, pos, mpi_rp, filetype, 'native', mpi_info_null, ie)
      if (ie /= 0) call error_stop('set_view')

      allocate(buffer(1:nx, 1:ny, 1:nz))
      buffer = bigend(phi(1:nx, 1:ny, 1:nz))

      call mpi_file_write_all(fh, buffer, nx * ny * nz, mpi_rp, mpi_status_ignore, ie)
      if (ie /= 0) call error_stop('write_all')

      call mpi_file_close(fh, ie)
      if (ie /= 0) call error_stop('close')

   end subroutine

end program
