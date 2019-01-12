module kinds
   use poisfft_precisions
   use iso_c_binding
   use iso_fortran_env

   integer, parameter :: rp = drp
   integer, parameter :: size_kind = c_size_t
end module

module globals
   use kinds
   implicit none
   real(rp), parameter :: pi = 4 * atan(1._rp) ! pi = 3.141592653589793238462_rp
   real(rp), parameter :: lx = 2 * pi, ly = lx * 1.1, lz = lx / 1.1
   integer :: nx = 21, ny = 32, nz = 25
   real(rp) :: dx, dy, dz
   real(rp), dimension(:, :, :), allocatable :: phi3d, rhs3d
   real(rp), dimension(:, :), allocatable :: phi2d, rhs2d
   real(rp), dimension(:), allocatable :: phi1d, rhs1d
end module

module residues
   use kinds
   use globals
   use poisfft
   implicit none

contains
   subroutine resexact1d_dir(phi, r)
      real(rp), intent(in) :: phi(0:)
      real(rp), intent(out) :: r
      integer :: i
      real(rp) :: x, p
      real(rp) :: xx, f
      xx(i) = dx * i
      !       f(x) = (x**3 - x*(nx*dx)**2)/6
      f(x) = -lx**2 / (5 * pi)**2 * sin(5 * pi * x / lx)
      !     f(x) = -(x/(3*lx**2))*(lx**3-2*lx*x**2+x**3)
      r = 0
      do i = 1, nx
         x = xx(i)
         p = abs(phi(i) - f(x))
         !         print *, x, f(x), phi(i)
         r = r + p**2
      end do
      r = sqrt(r) / nx
   end subroutine

   subroutine resexact1d_dirstag(phi, r)
      real(rp), intent(in) :: phi(0:)
      real(rp), intent(out) :: r
      integer :: i
      real(rp) :: x, p
      real(rp) :: xx, f
      xx(i) = dx / 2 + dx * (i - 1)
      !       f(x) = (x**3 - x*(nx*dx)**2)/6
      f(x) = -lx**2 / (5 * pi)**2 * sin(5 * pi * x / lx)
      !     f(x) = -(x/(3*lx**2))*(lx**3-2*lx*x**2+x**3)
      r = 0
      do i = 1, nx
         x = xx(i)
         p = abs(phi(i) - f(x))
         r = r + p**2
      end do
      r = sqrt(r) / nx
   end subroutine

   subroutine resexact1d_neum(phi, r)
      real(rp), intent(in) :: phi(0:)
      real(rp), intent(out) :: r
      integer :: i
      real(rp) :: x, p
      real(rp) :: xx, f, g
      !       xx(i) = dx/2 + dx*(i-1) - 0.5*lx
      !       g(x) = (x**3/6 - x*(nx*dx)**2/8)
      !       g(x) = (5._rp/16._rp)*lx**2*x+x**5/(5*lx**2)-x**3/2
      xx(i) = dx * (i - 1)
      g(x) = - (lx / (6 * pi))**2 * cos(6 * pi * (x / lx))
      f(x) = g(x) + (phi(1) - g(xx(1)))

      r = 0
      do i = 1, nx
         x = xx(i)
         p = abs(phi(i) - f(x))
         !         print *, x, f(x), phi(i)
         r = r + p**2
      end do
      r = sqrt(r) / nx
   end subroutine

   subroutine resexact1d_neumstag(phi, r)
      real(rp), intent(in) :: phi(0:)
      real(rp), intent(out) :: r
      integer :: i
      real(rp) :: x, p
      real(rp) :: xx, f, g
      !       xx(i) = dx/2 + dx*(i-1) - 0.5*lx
      !       g(x) = (x**3/6 - x*(nx*dx)**2/8)
      !       g(x) = (5._rp/16._rp)*lx**2*x+x**5/(5*lx**2)-x**3/2
      xx(i) = dx / 2 + dx * (i - 1)
      g(x) = - (lx / (6 * pi))**2 * cos(6 * pi * (x / lx))
      f(x) = g(x) + (phi(1) - g(xx(1)))

      r = 0
      do i = 1, nx
         x = xx(i)
         p = abs(phi(i) - f(x))
         r = r + p**2
      end do
      r = sqrt(r) / nx
   end subroutine

   subroutine resexact1d_per(phi, r)
      real(rp), intent(in) :: phi(0:)
      real(rp), intent(out) :: r
      integer :: i
      real(rp) :: x, p
      real(rp) :: xx, f, g
      xx(i) = dx * (i - 1)
      g(x) = - (lx / (6 * pi))**2 * cos(6 * pi * (x / lx))
      f(x) = g(x) + (phi(1) - g(xx(1)))

      r = 0
      do i = 1, nx
         x = xx(i)
         p = abs(phi(i) - f(x))
         !         print *, x, f(x), phi(i)
         r = r + p**2
      end do
      r = sqrt(r) / nx
   end subroutine

   subroutine resexact3d_dir(phi, r)
      real(rp), intent(in) :: phi(0:, 0:, 0:)
      real(rp), intent(out) :: r
      integer :: i, j, k
      real(rp) :: x, y, z, p
      real(rp) :: xx, yy, zz, f

      xx(i) = dx * i
      yy(j) = dy * j
      zz(k) = dz * k

      f(x, y, z) = -1 / ((3 * pi)**2 / lx**2 + (5 * pi)**2 / ly**2 + (7 * pi)**2 / lz**2) * &
         sin(3 * pi * x / lx) * sin(5 * pi * y / ly) * sin(7 * pi * z / lz)

      r = 0
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               x = xx(i)
               y = yy(j)
               z = zz(k)
               p = abs(phi(i, j, k) - f(x, y, z))
               r = r + p**2
            end do
         end do
      end do
      r = sqrt(r) / nx
   end subroutine

   subroutine resexact3d_dirstag(phi, r)
      real(rp), intent(in) :: phi(0:, 0:, 0:)
      real(rp), intent(out) :: r
      integer :: i, j, k
      real(rp) :: x, y, z, p
      real(rp) :: xx, yy, zz, f

      xx(i) = dx / 2 + dx * (i - 1)
      yy(j) = dy / 2 + dy * (j - 1)
      zz(k) = dz / 2 + dz * (k - 1)

      f(x, y, z) = -1 / ((3 * pi)**2 / lx**2 + (5 * pi)**2 / ly**2 + (7 * pi)**2 / lz**2) * &
         sin(3 * pi * x / lx) * sin(5 * pi * y / ly) * sin(7 * pi * z / lz)

      r = 0
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               x = xx(i)
               y = yy(j)
               z = zz(k)
               p = abs(phi(i, j, k) - f(x, y, z))
               r = r + p**2
            end do
         end do
      end do
      r = sqrt(r) / nx
   end subroutine

   subroutine res1d(phi, rhs, aw, ae, btype, r)
      !2nd order finite difference residuum

      integer, parameter :: ea = 1, we = 2

      real(rp), intent(inout) :: phi(0:)
      real(rp), intent(in) :: rhs(:)
      real(rp), intent(in) :: aw, ae
      integer, intent(in) :: btype(2)
      real(rp), intent(out) :: r
      integer i
      real(rp) :: p, ap

      r = 0

      do i = 1, nx
         p = 0
         ap = 0
         if (i > 1) then
            p = p + phi(i - 1) * aw
            ap = ap + aw
         elseif(btype(we) == poisfft_periodic) then
            p = p + phi(nx) * aw
            ap = ap + aw
         elseif(btype(we) == poisfft_dirichletstag) then
            p = p - phi(1) * aw
            ap = ap + aw
         elseif(btype(we) == poisfft_dirichlet) then
            ap = ap + aw
         elseif(btype(we) == poisfft_neumann) then
            p = p + phi(2) * aw
            ap = ap + aw
         end if
         if (i < nx) then
            p = p + phi(i + 1) * ae
            ap = ap + ae
         elseif(btype(ea) == poisfft_periodic) then
            p = p + phi(1) * ae
            ap = ap + ae
         elseif(btype(ea) == poisfft_dirichletstag) then
            p = p - phi(nx) * ae
            ap = ap + ae
         elseif(btype(ea) == poisfft_dirichlet) then
            ap = ap + ae
         elseif(btype(ea) == poisfft_neumann) then
            p = p + phi(nx - 1) * ae
            ap = ap + ae
         end if

         p = p - rhs(i)
         p = abs(-p + ap * phi(i))
         r = max(r, abs(p))
      end do
   end subroutine

   subroutine res2d(phi, rhs, aw, ae, as, an, btype, r)
      !2nd order finite difference residuum
      integer, parameter :: ea = 1, we = 2, so = 3, no = 4

      real(rp), intent(inout) :: phi(0:, 0:)
      real(rp), intent(in) :: rhs(:, :)
      real(rp), intent(in) :: aw, ae
      real(rp), intent(in) :: as, an
      integer, intent(in) :: btype(4)
      real(rp), intent(out) :: r
      integer i, j
      real(rp) :: p, ap

      r = 0

      do j = 1, ny
         do i = 1, nx
            p = 0
            ap = 0
            if (i > 1) then
               p = p + phi(i - 1, j) * aw
               ap = ap + aw
            elseif(btype(we) == poisfft_periodic) then
               p = p + phi(nx, j) * aw
               ap = ap + aw
            elseif(btype(we) == poisfft_dirichletstag) then
               p = p - phi(1, j) * aw
               ap = ap + aw
            elseif(btype(we) == poisfft_dirichlet) then
               ap = ap + aw
            elseif(btype(we) == poisfft_neumann) then
               p = p + phi(2, j) * aw
               ap = ap + aw
            end if
            if (i < nx) then
               p = p + phi(i + 1, j) * ae
               ap = ap + ae
            elseif(btype(ea) == poisfft_periodic) then
               p = p + phi(1, j) * ae
               ap = ap + ae
            elseif(btype(ea) == poisfft_dirichletstag) then
               p = p - phi(nx, j) * ae
               ap = ap + ae
            elseif(btype(ea) == poisfft_dirichlet) then
               ap = ap + ae
            elseif(btype(ea) == poisfft_neumann) then
               p = p + phi(nx - 1, j) * ae
               ap = ap + ae
            end if
            if (j > 1) then
               p = p + phi(i, j - 1) * as
               ap = ap + as
            elseif(btype(so) == poisfft_periodic) then
               p = p + phi(i, ny) * as
               ap = ap + as
            elseif(btype(so) == poisfft_dirichletstag) then
               p = p - phi(i, 1) * as
               ap = ap + as
            elseif(btype(so) == poisfft_dirichlet) then
               ap = ap + as
            elseif(btype(so) == poisfft_neumann) then
               p = p + phi(i, 2) * as
               ap = ap + as
            end if
            if (j < ny) then
               p = p + phi(i, j + 1) * an
               ap = ap + an
            elseif(btype(no) == poisfft_periodic) then
               p = p + phi(i, 1) * an
               ap = ap + an
            elseif(btype(no) == poisfft_dirichletstag) then
               p = p - phi(i, ny) * an
               ap = ap + an
            elseif(btype(so) == poisfft_dirichlet) then
               ap = ap + an
            elseif(btype(so) == poisfft_neumann) then
               p = p + phi(i, ny - 1) * an
               ap = ap + an
            end if

            p = p - rhs(i, j)
            p = abs(-p + ap * phi(i, j))
            r = max(r, abs(p))
         end do
      end do

   end subroutine

   subroutine res3d(phi, rhs, aw, ae, as, an, ab, at, btype, r)
      !2nd order finite difference residuum

      integer, parameter :: ea = 1, we = 2, so = 3, no = 4, bo = 5, to = 6

      real(rp), intent(inout) :: phi(0:, 0:, 0:)
      real(rp), intent(in) :: rhs(:, :, :)
      real(rp), intent(in) :: aw, ae
      real(rp), intent(in) :: as, an
      real(rp), intent(in) :: ab, at
      integer, intent(in) :: btype(6)
      real(rp), intent(out) :: r
      integer i, j, k
      real(rp) :: p, ap

      r = 0
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               p = 0
               ap = 0
               if (i > 1) then
                  p = p + phi(i - 1, j, k) * aw
                  ap = ap + aw
               elseif(btype(we) == poisfft_periodic) then
                  p = p + phi(nx, j, k) * aw
                  ap = ap + aw
               elseif(btype(we) == poisfft_dirichletstag) then
                  p = p - phi(1, j, k) * aw
                  ap = ap + aw
               elseif(btype(we) == poisfft_dirichlet) then
                  ap = ap + aw
               elseif(btype(we) == poisfft_neumann) then
                  p = p + phi(2, j, k) * aw
                  ap = ap + aw
               end if
               if (i < nx) then
                  p = p + phi(i + 1, j, k) * ae
                  ap = ap + ae
               elseif(btype(ea) == poisfft_periodic) then
                  p = p + phi(1, j, k) * ae
                  ap = ap + ae
               elseif(btype(ea) == poisfft_dirichletstag) then
                  p = p - phi(nx, j, k) * ae
                  ap = ap + ae
               elseif(btype(ea) == poisfft_dirichlet) then
                  ap = ap + ae
               elseif(btype(ea) == poisfft_neumann) then
                  p = p + phi(nx - 1, j, k) * ae
                  ap = ap + ae
               end if
               if (j > 1) then
                  p = p + phi(i, j - 1, k) * as
                  ap = ap + as
               elseif(btype(so) == poisfft_periodic) then
                  p = p + phi(i, ny, k) * as
                  ap = ap + as
               elseif(btype(so) == poisfft_dirichletstag) then
                  p = p - phi(i, 1, k) * as
                  ap = ap + as
               elseif(btype(so) == poisfft_dirichlet) then
                  ap = ap + as
               elseif(btype(so) == poisfft_neumann) then
                  p = p + phi(i, 2, k) * as
                  ap = ap + as
               end if
               if (j < ny) then
                  p = p + phi(i, j + 1, k) * an
                  ap = ap + an
               elseif(btype(no) == poisfft_periodic) then
                  p = p + phi(i, 1, k) * an
                  ap = ap + an
               elseif(btype(no) == poisfft_dirichletstag) then
                  p = p - phi(i, ny, k) * an
                  ap = ap + an
               elseif(btype(no) == poisfft_dirichlet) then
                  ap = ap + an
               elseif(btype(no) == poisfft_neumann) then
                  p = p + phi(i, ny - 1, k) * an
                  ap = ap + an
               end if
               if (k > 1) then
                  p = p + phi(i, j, k - 1) * ab
                  ap = ap + ab
               elseif(btype(bo) == poisfft_periodic) then
                  p = p + phi(i, j, nz) * ab
                  ap = ap + ab
               elseif(btype(bo) == poisfft_dirichletstag) then
                  p = p - phi(i, j, 1) * ab
                  ap = ap + ab
               elseif(btype(bo) == poisfft_dirichlet) then
                  ap = ap + ab
               elseif(btype(bo) == poisfft_neumann) then
                  p = p + phi(i, j, 2) * ab
                  ap = ap + ab
               end if
               if (k < nz) then
                  p = p + phi(i, j, k + 1) * at
                  ap = ap + at
               elseif(btype(to) == poisfft_periodic) then
                  p = p + phi(i, j, 1) * at
                  ap = ap + at
               elseif(btype(to) == poisfft_dirichletstag) then
                  p = p - phi(i, j, nz) * at
                  ap = ap + at
               elseif(btype(to) == poisfft_dirichlet) then
                  ap = ap + at
               elseif(btype(to) == poisfft_neumann) then
                  p = p + phi(i, j, nz - 1) * at
                  ap = ap + at
               end if

               p = p - rhs(i, j, k)
               p = abs(-p + ap * phi(i, j, k))
               r = max(r, abs(p))
            end do
         end do
      end do

   end subroutine
end module

module tests
   use kinds
   use globals
   use residues
   use poisfft, &
      poisfft_solver1d => poisfft_solver1d_dp, &
      poisfft_solver2d => poisfft_solver2d_dp, &
      poisfft_solver3d => poisfft_solver3d_dp
   implicit none
   character(20), parameter :: char_bcs(0:4) = &
      ['periodic            ', 'dirichlet regular   ', 'neumann regular     ', 'dirichlet staggered ', 'neumann staggered   ']
contains

   subroutine test1d(bcs, test_proc)
      integer, intent(in) :: bcs(2)
      integer :: i
      interface
         subroutine test_proc(phi, r)
            import
            real(rp), intent(in) :: phi(:)
            real(rp), intent(out) :: r
         end subroutine
      end interface

      write(*, *) '----'
      write(*, '(1x,"1d |",*(a,"|"))') (trim(char_bcs(bcs(i))), i=1, size(bcs))

      call testspectral1d(bcs, test_proc)
      call testfd2_1d(bcs)
      write(*, *)
   end subroutine

   subroutine testspectral1d(bcs, test_proc)
      integer, intent(in) :: bcs(2)
      real(rp) :: r

      interface
         subroutine test_proc(phi, r)
            import
            real(rp), intent(in) :: phi(:)
            real(rp), intent(out) :: r
         end subroutine
      end interface

      call runspectral1d(bcs)
      call test_proc(phi1d, r)

      if (r < nx * 10 * epsilon(1._rp)) then
         write(*, *) 'spectral ok'
      else
         write(*, *) 'spectral fail'
         write(*, *) 'spectral residuum:', r
      end if
   end subroutine

   subroutine runspectral1d(bcs)
      type(poisfft_solver1d) :: solver
      integer, intent(in) :: bcs(2)

      solver = poisfft_solver1d([nx], [lx], bcs)
      call execute(solver, phi1d, rhs1d)
      call finalize(solver)
   end subroutine

   subroutine testfd2_1d(bcs)
      integer, intent(in) :: bcs(2)
      real(rp) :: r

      call runfd2_1d(bcs)
      call res1d(phi1d, rhs1d, dx**(-2), dx**(-2), bcs, r)

      if (r < nx * 10 * epsilon(1._rp)) then
         write(*, *) 'fd2 ok', r
      else
         write(*, *) 'fd2 fail'
         write(*, *) 'fd2 residuum:', r
      end if
   end subroutine

   subroutine runfd2_1d(bcs)
      type(poisfft_solver1d) :: solver
      integer, intent(in) :: bcs(2)

      solver = poisfft_solver1d([nx], [lx], bcs, approximation=2)
      call execute(solver, phi1d, rhs1d)
      call finalize(solver)
   end subroutine

   subroutine test2d(bcs)!, test_pro
      integer, intent(in) :: bcs(4)
      integer :: i
      !     interface
      !       subroutine  test_proc(phi, r)
      !         import
      !         real(rp),intent(in) :: phi(:,:)
      !         real(rp), intent(out) :: r
      !       end subroutine
      !     end interface

      write(*, *) '----'
      write(*, '(1x,"2d |",*(a,"|"))') (trim(char_bcs(bcs(i))), i=1, size(bcs))

      !     call testspectral2d(bcs, test_proc)
      call testfd2_2d(bcs)
      write(*, *)
   end subroutine

   subroutine testspectral2d(bcs, test_proc)
      integer, intent(in) :: bcs(4)
      real(rp) :: r

      interface
         subroutine test_proc(phi, r)
            import
            real(rp), intent(in) :: phi(:, :)
            real(rp), intent(out) :: r
         end subroutine
      end interface

      call runspectral2d(bcs)
      call test_proc(phi2d, r)

      if (r < int(nx, int64) * int(ny, int64) * 10 * epsilon(1._rp)) then
         write(*, *) 'spectral ok'
      else
         write(*, *) 'spectral fail'
         write(*, *) 'spectral residuum:', r
      end if
   end subroutine

   subroutine runspectral2d(bcs)
      type(poisfft_solver2d) :: solver
      integer, intent(in) :: bcs(4)

      solver = poisfft_solver2d([nx, ny], [lx, ly], bcs)
      call execute(solver, phi2d, rhs2d)
      call finalize(solver)
   end subroutine


   subroutine testfd2_2d(bcs)
      integer, intent(in) :: bcs(4)
      real(rp) :: r

      call runfd2_2d(bcs)

      call res2d(phi2d, rhs2d, dx**(-2), dx**(-2), dy**(-2), dy**(-2), bcs, r)
      if (r < int(nx, int64) * int(ny, int64) * 10 * epsilon(1._rp)) then
         write(*, *) 'fd2 ok', r
      else
         write(*, *) 'fd2 fail'
         write(*, *) 'fd2 residuum:', r
      end if
   end subroutine

   subroutine runfd2_2d(bcs)
      type(poisfft_solver2d) :: solver
      integer, intent(in) :: bcs(4)

      solver = poisfft_solver2d([nx, ny], [lx, ly], bcs, approximation=2)
      call execute(solver, phi2d, rhs2d)
      call finalize(solver)
   end subroutine

   subroutine test3d(bcs, test_proc)
      integer, intent(in) :: bcs(6)
      integer :: i
      procedure(resexact3d_dir), optional :: test_proc

      write(*, *) '----'
      write(*, '(1x,"3d |",*(a,"|"))') (trim(char_bcs(bcs(i))), i=1, size(bcs))

      if (present(test_proc)) call testspectral3d(bcs, test_proc)
      call testfd2_3d(bcs)
      write(*, *)
   end subroutine

   subroutine testspectral3d(bcs, test_proc)
      integer, intent(in) :: bcs(6)
      real(rp) :: r

      interface
         subroutine test_proc(phi, r)
            import
            real(rp), intent(in) :: phi(:, :, :)
            real(rp), intent(out) :: r
         end subroutine
      end interface

      call runspectral3d(bcs)
      call test_proc(phi3d, r)

      if (r < int(nx, int64) * int(ny, int64) * int(nz, int64) * 10 * epsilon(1._rp)) then
         write(*, *) 'spectral ok'
      else
         write(*, *) 'spectral fail'
         write(*, *) 'spectral residuum:', r
      end if
   end subroutine

   subroutine runspectral3d(bcs)
      type(poisfft_solver3d) :: solver
      integer, intent(in) :: bcs(6)

      solver = poisfft_solver3d([nx, ny, nz], [lx, ly, lz], bcs)
      call execute(solver, phi3d, rhs3d)
      call finalize(solver)
   end subroutine


   subroutine testfd2_3d(bcs)
      integer, intent(in) :: bcs(6)
      real(rp) :: r

      call runfd2_3d(bcs)

      call res3d(phi3d, rhs3d, dx**(-2), dx**(-2), dy**(-2), dy**(-2), dz**(-2), dz**(-2), bcs, r)
      if (r < int(nx, int64) * int(ny, int64) * int(nz, int64) * 10 * epsilon(1._rp)) then
         write(*, *) 'fd2 ok', r
      else
         write(*, *) 'fd2 fail'
         write(*, *) 'fd2 residuum:', r
      end if
   end subroutine

   subroutine runfd2_3d(bcs)
      type(poisfft_solver3d) :: solver
      integer, intent(in) :: bcs(6)

      solver = poisfft_solver3d([nx, ny, nz], [lx, ly, lz], bcs, approximation=2)
      call execute(solver, phi3d, rhs3d)
      call finalize(solver)
   end subroutine
end module

program testpoisson
   use kinds
   use globals
   use residues
   use tests
   use poisfft
   implicit none

   integer i, j, k, niters
   real(rp) :: x, y, z
   character(len=12) :: arg
   real(rp) :: avg, p

   if (command_argument_count() >= 3) then
      call get_command_argument(1, value=arg)
      read(arg, '(i12)') nx
      call get_command_argument(2, value=arg)
      read(arg, '(i12)') ny
      call get_command_argument(3, value=arg)
      read(arg, '(i12)') nz
   end if

   write(*, *) 'nx,ny,nz', [nx, ny, nz]
   dx = lx / nx; dy = lx / ny; dz = lx / nz

   allocate(rhs3d(nx, ny, nz), phi3d(0:nx + 1, 0:ny + 1, 0:nz + 1))
   allocate(rhs2d(nx, ny), phi2d(0:nx + 1, 0:ny + 1))
   allocate(rhs1d(nx), phi1d(0:nx + 1))

   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            x = (i - 1._rp / 2) * dx
            y = (j - 1._rp / 2) * dy
            z = (k - 1._rp / 2) * dz
            call random_number(rhs3d(i, j, k))
            call random_number(phi3d(i, j, k))
         end do
      end do
   end do
   rhs3d = rhs3d - sum(rhs3d) / size(rhs3d)

   do j = 1, ny
      do i = 1, nx
         x = (i - 1._rp / 2) * dx
         y = (j - 1._rp / 2) * dy
         z = (k - 1._rp / 2) * dz
         call random_number(rhs2d(i, j))
         call random_number(phi2d(i, j))
      end do
   end do
   rhs2d = rhs2d - sum(rhs2d) / size(rhs2d)

   do i = 1, nx
      x = (i - 1._rp / 2) * dx
      call random_number(rhs1d(i))
      call random_number(phi1d(i))
   end do
   rhs1d = rhs1d - sum(rhs1d) / size(rhs1d)

   dx = lx / nx
   do i = 1, nx
      x = dx * (i - 1)
      rhs1d(i) = cos(3 * 2 * pi * (x / lx)) !+ cos(13*2*pi*(x/lx))
   end do
   call test1d([(poisfft_periodic, i=1, 2)], resexact1d_per)

   dx = lx / (nx + 1)
   !   rhs1d = [(dx/2 + dx*(i-1), i = 1,nx)]
   rhs1d = [(sin(5 * pi * (dx * (i)) / lx), i=1, nx)]
   !   rhs1d = [(  -(dx*i-lx/2)**2/(lx/2)**2+1, i = 1,nx)]
   call test1d([(poisfft_dirichlet, i=1, 2)], resexact1d_dir)

   dx = lx / nx
   rhs1d = [(sin(5 * pi * (dx * (i - 0.5_rp)) / lx), i=1, nx)]
   !   rhs1d = [(  -(dx/2 + dx*(i-1)-lx/2)**2/(lx/2)**2+1, i = 1,nx)]
   call test1d([(poisfft_dirichletstag, i=1, 2)], resexact1d_dirstag)

   dx = lx / (nx - 1)
   do i = 1, nx
      x = dx * (i - 1)
      rhs1d(i) = cos(3 * 2 * pi * (x / lx)) !+ cos(13*2*pi*(x/lx))
   end do
   rhs1d = rhs1d - (sum(rhs1d(2:nx - 1)) + rhs1d(1) / 2 + rhs1d(nx) / 2) / (size(rhs1d, kind=size_kind) - 1)
   call test1d([(poisfft_neumann, i=1, 2)], resexact1d_neum)

   dx = lx / nx
   do i = 1, nx
      x = dx / 2 + dx * (i - 1)
      rhs1d(i) = cos(3 * 2 * pi * (x / lx)) !+ cos(13*2*pi*(x/lx))
   end do
   rhs1d = rhs1d - sum(rhs1d) / size(rhs1d, kind=size_kind)
   call test1d([(poisfft_neumannstag, i=1, 2)], resexact1d_neumstag)

   dx = lx / nx
   dy = ly / ny
   call test2d([(poisfft_periodic, i=1, 4)])

   dx = lx / (nx + 1)
   dy = ly / (ny + 1)
   call test2d([(poisfft_dirichlet, i=1, 4)])

   dx = lx / nx
   dy = ly / ny
   call test2d([(poisfft_dirichletstag, i=1, 4)])

   dx = lx / (nx - 1)
   dy = ly / (ny - 1)
   avg = 0
   do j = 1, ny
      do i = 1, nx
         p = rhs2d(i, j)
         if (i == 1 .or. i == nx) p = p / 2
         if (j == 1 .or. j == ny) p = p / 2
         avg = avg + p
      end do
   end do
   rhs2d = rhs2d - avg / (real(nx - 1, rp) * real(ny - 1, rp))
   call test2d([(poisfft_neumann, i=1, 4)])

   dx = lx / nx
   dy = ly / ny
   rhs2d = rhs2d - sum(rhs2d) / (size(rhs2d, kind=size_kind))
   call test2d([(poisfft_neumannstag, i=1, 4)])

   dx = lx / nx
   dy = ly / ny
   dz = lz / nz
   call test3d([(poisfft_periodic, i=1, 6)])

   dx = lx / (nx + 1)
   dy = ly / (ny + 1)
   dz = lz / (nz + 1)
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            x = dx * i; y = dy * j; z = dz * k
            rhs3d(i, j, k) = sin(3 * pi * x / lx) * sin(5 * pi * y / ly) * sin(7 * pi * z / lz)
         end do
      end do
   end do
   call test3d([(poisfft_dirichlet, i=1, 6)], resexact3d_dir)

   dx = lx / nx
   dy = ly / ny
   dz = lz / nz
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            x = dx * (i - 0.5_rp); y = dy * (j - 0.5_rp); z = dz * (k - 0.5_rp)
            rhs3d(i, j, k) = sin(3 * pi * x / lx) * sin(5 * pi * y / ly) * sin(7 * pi * z / lz)
         end do
      end do
   end do
   call test3d([(poisfft_dirichletstag, i=1, 6)], resexact3d_dirstag)

   dx = lx / (nx - 1)
   dy = ly / (ny - 1)
   dz = lz / (nz - 1)
   avg = 0
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            p = rhs3d(i, j, k)
            if (i == 1 .or. i == nx) p = p / 2
            if (j == 1 .or. j == ny) p = p / 2
            if (k == 1 .or. k == nz) p = p / 2
            avg = avg + p
         end do
      end do
   end do
   rhs3d = rhs3d - avg / (real(nx - 1, rp) * real(ny - 1, rp) * real(nz - 1, rp))
   call test3d([(poisfft_neumann, i=1, 6)])

   dx = lx / nx; dy = ly / ny; dz = lz / nz
   rhs3d = rhs3d - sum(rhs3d) / (size(rhs3d, kind=size_kind))
   call test3d([(poisfft_neumannstag, i=1, 6)])

   dx = lx / nx; dy = ly / ny; dz = lz / nz
   call test3d([(poisfft_periodic, i=1, 4), (poisfft_neumannstag, i=5, 6)])

   dx = lx / nx; dy = ly / ny; dz = lz / nz
   call test3d([(poisfft_periodic, i=1, 2), (poisfft_neumannstag, i=3, 4), (poisfft_periodic, i=5, 6)])

   dx = lx / nx; dy = ly / ny; dz = lz / nz
   call test3d([(poisfft_periodic, i=1, 2), (poisfft_neumannstag, i=3, 6)])

end program
