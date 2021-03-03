module kinds
   use poisfft_constants
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



module Residues
  use Kinds
  use Globals
  use PoisFFT
    
  implicit none
  
  integer, parameter :: We = 1,Ea = 2,So = 3,No = 4,Bo = 5,To = 6

contains

  subroutine ResExact1D_Dir(Phi, R)
    real(rp),intent(in) :: Phi(0:)
    real(rp), intent(out) :: R
    integer :: i
    real(rp) :: x, p
    real(rp) :: xx,f
    xx(i) = dx*i
!       f(x) = (x**3 - x*(nx*dx)**2)/6
    f(x) = -Lx**2/(5*pi)**2 * sin(5*pi*x/Lx)
!     f(x) = -(x/(3*Lx**2))*(Lx**3-2*Lx*x**2+x**3)
    R = 0
    do i = 1, nx
      x = xx(i)
      p = abs( Phi(i) - f(x) )
!         print *, x, f(x), Phi(i)
      R = R + p**2
    end do
    R = sqrt(R/nx)
  end subroutine

  subroutine ResExact1D_DirStag(Phi, R)
    real(rp), intent(in) :: Phi(0:)
    real(rp), intent(out) :: R
    integer :: i
    real(rp) :: x, p
    real(rp) :: xx,f
    xx(i) = dx/2 + dx*(i-1)
!       f(x) = (x**3 - x*(nx*dx)**2)/6
    f(x) = -Lx**2/(5*pi)**2 * sin(5*pi*x/Lx)
!     f(x) = -(x/(3*Lx**2))*(Lx**3-2*Lx*x**2+x**3)
    R = 0
    do i = 1, nx
      x = xx(i)
      p = abs( Phi(i) - f(x) )
      R = R + p**2
    end do
    R = sqrt(R/nx)
  end subroutine

  subroutine ResExact1D_Neum(Phi, R)
    real(rp), intent(in) :: Phi(0:)
    real(rp), intent(out) :: R
    integer :: i
    real(rp) :: x, p
    real(rp) :: xx,f,g
!       xx(i) = dx/2 + dx*(i-1) - 0.5*Lx
!       g(x) = (x**3/6 - x*(nx*dx)**2/8)
!       g(x) = (5._rp/16._rp)*Lx**2*x+x**5/(5*Lx**2)-x**3/2
    xx(i) = dx*(i-1)
    g(x) = -(Lx/(6*pi))**2 * cos(6*pi*(x/Lx))
    f(x) = g(x) + (Phi(1)-g(xx(1)))

    R = 0
    do i = 1, nx
      x = xx(i)
      p = abs( Phi(i) - f(x) )
!         print *, x, f(x), Phi(i)
      R = R + p**2
    end do
    R = sqrt(R/nx)
  end subroutine


  subroutine ResExact1D_NeumStag(Phi, R)
    real(rp), intent(in) :: Phi(0:)
    real(rp), intent(out) :: R
    integer :: i
    real(rp) :: x, p
    real(rp) :: xx,f,g
!       xx(i) = dx/2 + dx*(i-1) - 0.5*Lx
!       g(x) = (x**3/6 - x*(nx*dx)**2/8)
!       g(x) = (5._rp/16._rp)*Lx**2*x+x**5/(5*Lx**2)-x**3/2
    xx(i) = dx/2 + dx*(i-1)
    g(x) = -(Lx/(6*pi))**2 * cos(6*pi*(x/Lx))
    f(x) = g(x) + (Phi(1)-g(xx(1)))

    R = 0
    do i = 1, nx
      x = xx(i)
      p = abs( Phi(i) - f(x) )
      R = R + p**2
    end do
    R = sqrt(R/nx)
  end subroutine


  subroutine ResExact1D_NeumStagDirStag(Phi, R)
    real(rp), intent(in) :: Phi(0:)
    real(rp), intent(out) :: R
    integer :: i
    real(rp) :: x, p
    real(rp) :: xx,f
!       xx(i) = dx/2 + dx*(i-1) - 0.5*Lx
!       g(x) = (x**3/6 - x*(nx*dx)**2/8)
!       g(x) = (5._rp/16._rp)*Lx**2*x+x**5/(5*Lx**2)-x**3/2
    xx(i) = dx/2 + dx*(i-1)
    f(x) = -(Lx/(1.5*pi))**2 * cos(1.5*pi*(x/Lx))
 
    R = 0
    do i = 1, nx
      x = xx(i)
      p = abs( Phi(i) - f(x) )  
      R = R + p**2
    end do
    R = sqrt(R/nx)
  end subroutine


  subroutine ResExact1D_DirStagNeumStag(Phi, R)
    real(rp), intent(in) :: Phi(0:)
    real(rp), intent(out) :: R
    integer :: i
    real(rp) :: x, p
    real(rp) :: xx,f
!       xx(i) = dx/2 + dx*(i-1) - 0.5*Lx
!       g(x) = (x**3/6 - x*(nx*dx)**2/8)
!       g(x) = (5._rp/16._rp)*Lx**2*x+x**5/(5*Lx**2)-x**3/2
    xx(i) = dx/2 + dx*(i-1)
    f(x) = -(Lx/(1.5*pi))**2 * sin(1.5*pi*(x/Lx))
 
    R = 0
    do i = 1, nx
      x = xx(i)
      p = abs( Phi(i) - f(x) )  
      R = R + p**2
    end do
    R = sqrt(R/nx)
  end subroutine


  subroutine ResExact1D_per(Phi, R)
    real(rp), intent(in) :: Phi(0:)
    real(rp), intent(out) :: R
    integer :: i
    real(rp) :: x, p
    real(rp) :: xx,f,g
    xx(i) = dx*(i-1)
    g(x) = -(Lx/(6*pi))**2 * cos(6*pi*(x/Lx))
    f(x) = g(x) + (Phi(1)-g(xx(1)))

    R = 0
    do i = 1, nx
      x = xx(i)
      p = abs( Phi(i) - f(x) )
!         print *, x, f(x), Phi(i)
      R = R + p**2
    end do
    R = sqrt(R/nx)
  end subroutine
  
  
  
  
  
  
  
  
  
  subroutine ResExact3D_Dir(Phi, R)
    real(rp),intent(in) :: Phi(0:,0:,0:)
    real(rp), intent(out) :: R
    integer :: i, j, k
    real(rp) :: x, y, z, p
    real(rp) :: xx, yy, zz, f
    
    xx(i) = dx*i
    yy(j) = dy*j
    zz(k) = dz*k

    f(x,y,z) = - 1 / ( (3*pi)**2/Lx**2 + (5*pi)**2/Ly**2 + (7*pi)**2/Lz**2 ) * &
                  sin(3*pi*x/Lx) * sin(5*pi*y/Ly) *sin(7*pi*z/Lz)

    R = 0
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          x = xx(i)
          y = yy(j)
          z = zz(k)
          p = abs( Phi(i, j, k) - f(x,y,z) )
          R = R + p**2
        end do
      end do
    end do
    R = sqrt(R/(nx*ny*nz))
  end subroutine

  subroutine ResExact3D_DirStag(Phi, R)
    real(rp), intent(in) :: Phi(0:,0:,0:)
    real(rp), intent(out) :: R
    integer :: i, j, k
    real(rp) :: x, y, z, p
    real(rp) :: xx, yy, zz, f
    
    xx(i) = dx/2 + dx*(i-1)
    yy(j) = dy/2 + dy*(j-1)
    zz(k) = dz/2 + dz*(k-1)
    
    f(x,y,z) = - 1 / ( (3*pi)**2/Lx**2 + (5*pi)**2/Ly**2 + (7*pi)**2/Lz**2 ) * &
                  sin(3*pi*x/Lx) * sin(5*pi*y/Ly) *sin(7*pi*z/Lz)
                  
    R = 0
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          x = xx(i)
          y = yy(j)
          z = zz(k)
          p = abs( Phi(i, j, k) - f(x,y,z) )
          R = R + p**2
        end do
      end do
    end do
    R = sqrt(R/(nx*ny*nz))
  end subroutine


  subroutine ResExact3D_NeumStag(Phi, R)
    real(rp), intent(in) :: Phi(0:,0:,0:)
    real(rp), intent(out) :: R
    integer :: i, j, k
    real(rp) :: x, y, z, p
    real(rp) :: xx, yy, zz, f
    
    xx(i) = dx/2 + dx*(i-1)
    yy(j) = dy/2 + dy*(j-1)
    zz(k) = dz/2 + dz*(k-1)
    
    f(x,y,z) = - 1 / ( (3*pi)**2/Lx**2 + (5*pi)**2/Ly**2 + (7*pi)**2/Lz**2 ) * &
                  cos(3*pi*x/Lx) * cos(5*pi*y/Ly) *cos(7*pi*z/Lz)
                  
    R = 0
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          x = xx(i)
          y = yy(j)
          z = zz(k)
          p = abs( Phi(i, j, k) - f(x,y,z) )
          R = R + p**2
        end do
      end do
    end do
    R = sqrt(R/(nx*ny*nz))
  end subroutine


  
  subroutine ResExact3D_DsDsP(Phi, R)
    real(rp), intent(in) :: Phi(0:,0:,0:)
    real(rp), intent(out) :: R
    integer :: i, j, k
    real(rp) :: x, y, z, p
    real(rp) :: xx, yy, zz, f
    
    xx(i) = dx/2 + dx*(i-1)
    yy(j) = dy/2 + dy*(j-1)
    zz(k) = dz/2 + dz*(k-1)
    
    f(x,y,z) = - 1 / ( (3*pi)**2/Lx**2 + (5*pi)**2/Ly**2 + (6*pi)**2/Lz**2 ) * &
                  sin(3*pi*x/Lx) * sin(5*pi*y/Ly) *sin(6*pi*z/Lz)
                  
    R = 0
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          x = xx(i)
          y = yy(j)
          z = zz(k)
          p = abs( Phi(i, j, k) - f(x,y,z) )
          R = R + p**2
        end do
      end do
    end do
    R = sqrt(R/(nx*ny*nz))
  end subroutine


  subroutine ResExact3D_DsDsNs(Phi, R)
    real(rp), intent(in) :: Phi(0:,0:,0:)
    real(rp), intent(out) :: R
    integer :: i, j, k
    real(rp) :: x, y, z, p
    real(rp) :: xx, yy, zz, f
    
    xx(i) = dx/2 + dx*(i-1)
    yy(j) = dy/2 + dy*(j-1)
    zz(k) = dz/2 + dz*(k-1)
    
    f(x,y,z) = - 1 / ( (3*pi)**2/Lx**2 + (5*pi)**2/Ly**2 + (6*pi)**2/Lz**2 ) * &
                  sin(3*pi*x/Lx) * sin(5*pi*y/Ly) *cos(6*pi*z/Lz)
                  
    R = 0
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          x = xx(i)
          y = yy(j)
          z = zz(k)
          p = abs( Phi(i, j, k) - f(x,y,z) )
          R = R + p**2
        end do
      end do
    end do
    R = sqrt(R/(nx*ny*nz))
  end subroutine


  
  
  
  
  
  
  
  
  
  
  
  


  subroutine Res1D(Phi,RHS,&
                   Aw,Ae,&
                   Btype,R)
    !2nd order finite difference residuum

    real(rp), intent(inout) :: Phi(0:)
    real(rp), intent(in) :: RHS(:)
    real(rp),intent(in) :: Aw,Ae
    integer,intent(in) :: Btype(2)
    real(rp),intent(out) :: R
    integer i
    real(rp) :: p,Ap

    R = 0

    do i = 1, nx
      p = 0
      Ap = 0
      if (i>1) then
                p = p + Phi(i-1)*Aw
                Ap = Ap + Aw
      elseif (Btype(We)==PoisFFT_Periodic) then
                p = p + Phi(nx)*Aw
                Ap = Ap + Aw
      elseif (Btype(We)==PoisFFT_DirichletStag) then
                p = p - Phi(1)*Aw
                Ap = Ap + Aw
      elseif (Btype(We)==PoisFFT_Dirichlet) then
                Ap = Ap + Aw
      elseif (Btype(We)==PoisFFT_Neumann) then
                p = p + Phi(2)*Aw
                Ap = Ap + Aw
      end if
      if (i<nx) then
                p = p + Phi(i+1)*Ae
                Ap = Ap + Ae
      elseif (Btype(Ea)==PoisFFT_Periodic) then
                p = p + Phi(1)*Ae
                Ap = Ap + Ae
      elseif (Btype(Ea)==PoisFFT_DirichletStag) then
                p = p - Phi(nx)*Ae
                Ap = Ap + Ae
      elseif (Btype(Ea)==PoisFFT_Dirichlet) then
                Ap = Ap + Ae
      elseif (Btype(Ea)==PoisFFT_Neumann) then
                p = p + Phi(nx-1)*Ae
                Ap = Ap + Ae
      end if

      p = p - RHS(i)   
      p = abs(-p +Ap*Phi(i))      
      R = max(R,abs(p))
    end do

  end subroutine Res1D










  subroutine Res2D(Phi,RHS,&
                   Aw,Ae,As,An,&
                   Btype,R)
    !2nd order finite difference residuum

    real(rp), intent(inout) :: Phi(0:,0:)
    real(rp), intent(in) :: RHS(:,:)
    real(rp),intent(in) :: Aw,Ae
    real(rp),intent(in) :: As,An
    integer,intent(in) :: Btype(4)
    real(rp),intent(out) :: R
    integer i,j
    real(rp) :: p,Ap

    R = 0

    do j = 1, ny
      do i = 1, nx
         x = xx(i)
         p = abs(phi(i) - f(x))
         ! print *, x, f(x), phi(i)
         r = r + p**2
      end do
      r = sqrt(r) / nx
   end subroutine

   subroutine resexact1d_neumstag(phi, r)
      real(rp), intent(in) :: phi(0:)
      real(rp), intent(out) :: r
      integer :: i
      real(rp) :: x, p, xx, f, g
      ! xx(i) = dx/2 + dx*(i-1) - 0.5*lx
      ! g(x) = (x**3/6 - x*(nx*dx)**2/8)
      ! g(x) = (5._rp/16._rp)*lx**2*x+x**5/(5*lx**2)-x**3/2
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
         ! print *, x, f(x), phi(i)
         r = r + p**2
      end do
      r = sqrt(r) / nx
   end subroutine

   subroutine resexact3d_dir(phi, r)
      real(rp), intent(in) :: phi(0:, 0:, 0:)
      real(rp), intent(out) :: r
      integer :: i, j, k
      real(rp) :: x, y, z, p, xx, yy, zz, f

      xx(i) = dx * i
      yy(j) = dy * j
      zz(k) = dz * k
      f(x, y, z) = -1 / ((3 * pi)**2 / lx**2 + (5 * pi)**2 / ly**2 + (7 * pi)**2 / lz**2) * &
         sin(3 * pi * x / lx) * sin(5 * pi * y / ly) * sin(7 * pi * z / lz)

      r = 0
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               x = xx(i); y = yy(j); z = zz(k)
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
      real(rp) :: x, y, z, p, xx, yy, zz, f

      xx(i) = dx / 2 + dx * (i - 1)
      yy(j) = dy / 2 + dy * (j - 1)
      zz(k) = dz / 2 + dz * (k - 1)
      f(x, y, z) = -1 / ((3 * pi)**2 / lx**2 + (5 * pi)**2 / ly**2 + (7 * pi)**2 / lz**2) * &
         sin(3 * pi * x / lx) * sin(5 * pi * y / ly) * sin(7 * pi * z / lz)

      r = 0
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               x = xx(i); y = yy(j); z = zz(k)
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

    real(rp), intent(inout) :: Phi(0:,0:,0:)
    real(rp), intent(in) :: RHS(:,:,:)
    real(rp),intent(in) :: Aw,Ae
    real(rp),intent(in) :: As,An
    real(rp),intent(in) :: Ab,At
    integer,intent(in) :: Btype(6)
    real(rp),intent(out) :: R
    integer i,j,k
    real(rp) :: p,Ap

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
      ! 2nd order finite difference residuum
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

module Tests
  use Kinds
  use Globals
  use Residues
  use PoisFFT, PoisFFT_Solver1D => PoisFFT_Solver1D_DP, &
               PoisFFT_Solver2D => PoisFFT_Solver2D_DP, &
               PoisFFT_Solver3D => PoisFFT_Solver3D_DP  
  
  implicit none
  
  character(20), parameter :: char_BCs(0:4) = ["Periodic            ", &
                                             "Dirichlet regular   ", &
                                             "Neumann regular     ", &
                                             "Dirichlet staggered ", &
                                             "Neumann staggered   "]
 
contains
  
  subroutine Test1D(BCs, test_proc)
    integer, intent(in) :: BCs(2)
    integer :: i
    interface
      subroutine  test_proc(Phi, R)
        import
        real(rp),intent(in) :: Phi(:)
        real(rp), intent(out) :: R
      end subroutine
    end interface
    
    write(*,*) "----"
    write(*,'(1x,"1D  |",*(a,"|"))') (trim(char_BCs(BCs(i))), i=1,size(BCs))
    
    call TestSpectral1D(BCs, test_proc)
    call TestFD2_1D(BCs)
    write(*,*)
  end subroutine
  
  subroutine TestSpectral1D(BCs, test_proc)
    integer, intent(in) :: BCs(2)
    real(rp) :: R
    
    interface
      subroutine  test_proc(Phi, R)
        import
        real(rp),intent(in) :: Phi(:)
        real(rp), intent(out) :: R
      end subroutine
    end interface
    
    call RunSpectral1D(BCs)
    
    call test_proc(Phi1D, R)
    
    if (R < nx * 100 * epsilon(1._rp)) then
      write(*,*) "Spectral OK"
    else
      write(*,*) "Spectral FAIL"
      write(*,*) "Spectral residuum:",R
    end if
  end subroutine
  
  subroutine RunSpectral1D(BCs)
    type(PoisFFT_Solver1D) :: Solver
    integer, intent(in) :: BCs(2)
    
    Solver = PoisFFT_Solver1D([nx],[Lx],BCs)

    call Execute(Solver, Phi1D, RHS1D)

    call Finalize(Solver)
  end subroutine

  
  subroutine TestFD2_1D(BCs)
    integer, intent(in) :: BCs(2)
    real(rp) :: R
    
    call RunFD2_1D(BCs)
    
    call Res1D(Phi1D, RHS1D, &
               dx**(-2), dx**(-2), &
               BCs, R)
    
    if (R < nx * 100 * epsilon(1._rp)) then
      write(*,*) "FD2 OK", R
    else
      write(*,*) "FD2 FAIL"
      write(*,*) "FD2 residuum:",R
    end if
  end subroutine
  
  subroutine RunFD2_1D(BCs)
    type(PoisFFT_Solver1D) :: Solver
    integer, intent(in) :: BCs(2)
    
    Solver = PoisFFT_Solver1D([nx],[Lx],BCs, approximation=2)

    call Execute(Solver, Phi1D, RHS1D)

    call Finalize(Solver)
  end subroutine
  
  
  
  
  
  
  
  
  
  subroutine Test2D(BCs)!, test_pro
    integer, intent(in) :: BCs(4)
    integer :: i
!     interface
!       subroutine  test_proc(Phi, R)
!         import
!         real(rp),intent(in) :: Phi(:,:)
!         real(rp), intent(out) :: R
!       end subroutine
!     end interface
    
    write(*,*) "----"
    write(*,'(1x,"2D  |",*(a,"|"))') (trim(char_BCs(BCs(i))), i=1,size(BCs))
    
!     call TestSpectral2D(BCs, test_proc)
    call TestFD2_2D(BCs)
    write(*,*)
  end subroutine
  
  subroutine TestSpectral2D(BCs, test_proc)
    integer, intent(in) :: BCs(4)
    real(rp) :: R
    
    interface
      subroutine  test_proc(Phi, R)
        import
        real(rp),intent(in) :: Phi(:,:)
        real(rp), intent(out) :: R
      end subroutine
    end interface
    
    call RunSpectral2D(BCs)
    
    call test_proc(Phi2D, R)
    
    if (R < int(nx, int64) * int(ny, int64) * 100 * epsilon(1._rp)) then
      write(*,*) "Spectral OK"
    else
      write(*,*) "Spectral FAIL"
      write(*,*) "Spectral residuum:",R
    end if
  end subroutine
  
  subroutine RunSpectral2D(BCs)
    type(PoisFFT_Solver2D) :: Solver
    integer, intent(in) :: BCs(4)
    
    Solver = PoisFFT_Solver2D([nx,ny],[Lx,Ly],BCs)

    call Execute(Solver, Phi2D, RHS2D)

    call Finalize(Solver)
  end subroutine

  
  subroutine TestFD2_2D(BCs)
    integer, intent(in) :: BCs(4)
    real(rp) :: R
    
    call RunFD2_2D(BCs)
    
    call Res2D(Phi2D, RHS2D, &
               dx**(-2), dx**(-2), dy**(-2), dy**(-2), &
               BCs, R)
    if (R < int(nx, int64) * int(ny, int64) * 100 * epsilon(1._rp)) then
      write(*,*) "FD2 OK", R
    else
      write(*,*) "FD2 FAIL"
      write(*,*) "FD2 residuum:",R
    end if
  end subroutine
  
  subroutine RunFD2_2D(BCs)
    type(PoisFFT_Solver2D) :: Solver
    integer, intent(in) :: BCs(4)
    
    Solver = PoisFFT_Solver2D([nx,ny],[Lx,Ly],BCs, approximation=2)

    call Execute(Solver, Phi2D, RHS2D)

    call Finalize(Solver)
  end subroutine
  
  
  
  
  
  
  
  
  
  
  
  
  subroutine Test3D(BCs, test_proc)
    integer, intent(in) :: BCs(6)
    integer :: i
    procedure(ResExact3D_Dir), optional :: test_proc
    
    write(*,*) "----"
    write(*,'(1x,"3D  |",*(a,"|"))') (trim(char_BCs(BCs(i))), i=1,size(BCs))
    
    if (present(test_proc)) call TestSpectral3D(BCs, test_proc)
    call TestFD2_3D(BCs)
    write(*,*)
  end subroutine
  
  subroutine TestSpectral3D(BCs, test_proc)
    integer, intent(in) :: BCs(6)
    real(rp) :: R
    
    interface
      subroutine  test_proc(Phi, R)
        import
        real(rp),intent(in) :: Phi(:,:,:)
        real(rp), intent(out) :: R
      end subroutine
    end interface
    
    call RunSpectral3D(BCs)
    
    call test_proc(Phi3D, R)
    
    if (R < int(nx, int64) * int(ny, int64) * int(nz, int64) * 100 * epsilon(1._rp)) then
      write(*,*) "Spectral OK"
    else
      write(*,*) "Spectral FAIL"
      write(*,*) "Spectral residuum:",R
    end if
  end subroutine
  
  subroutine RunSpectral3D(BCs)
    type(PoisFFT_Solver3D) :: Solver
    integer, intent(in) :: BCs(6)
    
    Solver = PoisFFT_Solver3D([nx,ny,nz], [Lx,Ly,Lz], BCs)

    call Execute(Solver, Phi3D, RHS3D)

    call Finalize(Solver)
  end subroutine

  
  subroutine TestFD2_3D(BCs)
    integer, intent(in) :: BCs(6)
    real(rp) :: R

    call RunFD2_3D(BCs)
   
    call Res3D(Phi3D, RHS3D, &
               dx**(-2), dx**(-2), dy**(-2), dy**(-2), dz**(-2), dz**(-2), &
               BCs, R)
    if (R < int(nx, int64) * int(ny, int64) * int(nz, int64) * 100 * epsilon(1._rp)) then
      write(*,*) "FD2 OK", R
    else
      write(*,*) "FD2 FAIL"
      write(*,*) "FD2 residuum:",R
    end if
  end subroutine
  
  subroutine RunFD2_3D(BCs)
    type(PoisFFT_Solver3D) :: Solver
    integer, intent(in) :: BCs(6)
    
    Solver = PoisFFT_Solver3D([nx,ny,nz], [Lx,Ly,Lz], BCs, &
                              approximation=2)

    call Execute(Solver, Phi3D, RHS3D)

    call Finalize(Solver)
  end subroutine
  
end module

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

   integer :: i, j, k
   real(rp) :: x, y, z
   character(12) :: arg
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
  end do
  RHS3D = RHS3D-sum(RHS3D)/size(RHS3D)




  do j = 1,ny
    do i = 1,nx
     x=(i-1._rp/2)*dx
     y=(j-1._rp/2)*dy
     z=(k-1._rp/2)*dz
     call random_number(RHS2D(i,j))
     call random_number(Phi2D(i,j))
    end do
  end do
  RHS2D = RHS2D-sum(RHS2D)/size(RHS2D)



  do i = 1,nx
    x=(i-1._rp/2)*dx
    call random_number(RHS1D(i))
    call random_number(Phi1D(i))
  end do
  RHS1D = RHS1D-sum(RHS1D)/size(RHS1D)


  dx = Lx / nx
  do i = 1,nx
    x = dx*(i-1)
    RHS1D(i) = cos(3*2*pi*(x/Lx)) !+ cos(13*2*pi*(x/Lx))
  end do
  call Test1D([(PoisFFT_PERIODIC, i = 1,2)], ResExact1D_Per)

  
  dx = Lx / (nx+1)
!   RHS1D = [(dx/2 + dx*(i-1), i = 1,nx)]
  RHS1D = [(sin(5*pi*(dx*(i))/Lx), i = 1,nx)]
!   RHS1D = [(  -(dx*i-Lx/2)**2/(Lx/2)**2+1, i = 1,nx)]
  call Test1D([(PoisFFT_Dirichlet, i = 1,2)], ResExact1D_Dir)


  dx = Lx / nx
  RHS1D = [(sin(5*pi*(dx*(i-0.5_rp))/Lx), i = 1,nx)]
!   RHS1D = [(  -(dx/2 + dx*(i-1)-Lx/2)**2/(Lx/2)**2+1, i = 1,nx)]
  call Test1D([(PoisFFT_DirichletStag, i = 1,2)], ResExact1D_DirStag)


  dx = Lx / (nx-1)
  do i = 1,nx
    x = dx*(i-1)
    RHS1D(i) = cos(3*2*pi*(x/Lx)) !+ cos(13*2*pi*(x/Lx))
  end do
  RHS1D = RHS1D - (sum(RHS1D(2:nx-1))+RHS1D(1)/2+RHS1D(nx)/2)/(size(RHS1D,kind=size_kind)-1)
  call Test1D([(PoisFFT_Neumann, i = 1,2)], ResExact1D_Neum)


  dx = Lx / nx
  do i = 1,nx
    x = dx/2 + dx*(i-1)
    RHS1D(i) = cos(3*2*pi*(x/Lx)) !+ cos(13*2*pi*(x/Lx))
  end do
  RHS1D = RHS1D - sum(RHS1D)/size(RHS1D,kind=size_kind)
  call Test1D([(PoisFFT_NeumannStag, i = 1,2)], ResExact1D_NeumStag)

  
  dx = Lx / nx
  do i = 1,nx
    x = dx/2 + dx*(i-1)
    RHS1D(i) = cos(0.75_rp*2*pi*(x/Lx))
  end do
  call Test1D([PoisFFT_NeumannStag, PoisFFT_DirichletStag], ResExact1D_NeumStagDirStag)


  dx = Lx / nx
  do i = 1,nx
    x = dx/2 + dx*(i-1)
    RHS1D(i) = sin(0.75_rp*2*pi*(x/Lx))
  end do
  call Test1D([PoisFFT_DirichletStag, PoisFFT_NeumannStag], ResExact1D_DirStagNeumStag)

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

  
  

  dx = Lx / nx
  dy = Ly / ny
  call Test2D([(PoisFFT_Periodic, i = 1,4)])

  dx = Lx / (nx+1)
  dy = Ly / (ny+1)
  call Test2D([(PoisFFT_Dirichlet, i = 1,4)])

  dx = Lx / nx
  dy = Ly / ny
  call Test2D([(PoisFFT_DirichletStag, i = 1,4)])

  dx = Lx / (nx-1)
  dy = Ly / (ny-1)
  avg = 0
  do j = 1,ny
    do i = 1,nx
      p = RHS2D(i,j)
      if (i==1.or.i==nx) p = p / 2
      if (j==1.or.j==ny) p = p / 2
      avg = avg + p
    end do
  end do
  RHS2D = RHS2D - avg/(real(nx-1,rp)*real(ny-1,rp))
  call Test2D([(PoisFFT_Neumann, i = 1,4)])

  dx = Lx / nx
  dy = Ly / ny
  RHS2D = RHS2D - sum(RHS2D)/(size(RHS2D,kind=size_kind))
  call Test2D([(PoisFFT_NeumannStag, i = 1,4)])


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

   dx = lx / nx; dy = ly / ny
   call test2d([(poisfft_periodic, i=1, 4)])

   dx = lx / (nx + 1); dy = ly / (ny + 1)
   call test2d([(poisfft_dirichlet, i=1, 4)])

   dx = lx / nx; dy = ly / ny
   call test2d([(poisfft_dirichletstag, i=1, 4)])

   dx = lx / (nx - 1); dy = ly / (ny - 1)
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

   dx = lx / nx; dy = ly / ny
   rhs2d = rhs2d - sum(rhs2d) / (size(rhs2d, kind=size_kind))
   call test2d([(poisfft_neumannstag, i=1, 4)])

   dx = lx / nx; dy = ly / ny; dz = lz / nz
   call test3d([(poisfft_periodic, i=1, 6)])

  dx = Lx / (nx+1)
  dy = Ly / (ny+1)
  dz = Lz / (nz+1)
  do k = 1,nz
    do j = 1,ny
      do i = 1,nx
        x = dx*i; y = dy*j; z = dz*k
        RHS3D(i,j,k) = sin(3*pi*x/Lx) * sin(5*pi*y/Ly) *sin(7*pi*z/Lz)
      end do
    end do
  end do
  call Test3D([(PoisFFT_Dirichlet, i = 1,6)], ResExact3D_Dir)


  dx = Lx / nx
  dy = Ly / ny
  dz = Lz / nz
  do k = 1,nz
    do j = 1,ny
      do i = 1,nx
        x = dx*(i-0.5_rp); y = dy*(j-0.5_rp); z = dz*(k-0.5_rp)
        RHS3D(i,j,k) = sin(3*pi*x/Lx) * sin(5*pi*y/Ly) *sin(6*pi*z/Lz)
      end do
    end do
  end do
  call Test3D([(PoisFFT_DirichletStag, i = 1,4),(PoisFFT_Periodic, i = 5,6)], ResExact3D_DsDsP)


  dx = Lx / nx
  dy = Ly / ny
  dz = Lz / nz
  do k = 1,nz
    do j = 1,ny
      do i = 1,nx
        x = dx*(i-0.5_rp); y = dy*(j-0.5_rp); z = dz*(k-0.5_rp)
        RHS3D(i,j,k) = sin(3*pi*x/Lx) * sin(5*pi*y/Ly) *cos(6*pi*z/Lz)
      end do
    end do
  end do
  call Test3D([(PoisFFT_DirichletStag, i = 1,4),(PoisFFT_NeumannStag, i = 5,6)], ResExact3D_DsDsNs)


  dx = Lx / nx
  dy = Ly / ny
  dz = Lz / nz
  do k = 1,nz
    do j = 1,ny
      do i = 1,nx
        x = dx*(i-0.5_rp); y = dy*(j-0.5_rp); z = dz*(k-0.5_rp)
        RHS3D(i,j,k) = sin(3*pi*x/Lx) * sin(5*pi*y/Ly) *cos(1.5*pi*z/Lz)
      end do
    end do
  end do
  call Test3D([(PoisFFT_DirichletStag, i = 1,4),PoisFFT_NeumannStag, PoisFFT_DirichletStag])


  dx = Lx / nx
  dy = Ly / ny
  dz = Lz / nz
  do k = 1,nz
    do j = 1,ny
      do i = 1,nx
        x = dx*(i-0.5_rp); y = dy*(j-0.5_rp); z = dz*(k-0.5_rp)
        RHS3D(i,j,k) = sin(3*pi*x/Lx) * sin(5*pi*y/Ly) *sin(7*pi*z/Lz)
      end do
    end do
  end do
  call Test3D([(PoisFFT_DirichletStag, i = 1,6)], ResExact3D_DirStag)


  dx = Lx / (nx-1)
  dy = Ly / (ny-1)
  dz = Lz / (nz-1)
  avg = 0
  do k = 1,nz
    do j = 1,ny
      do i = 1,nx
        p = RHS3D(i,j,k)
        if (i==1.or.i==nx) p = p / 2
        if (j==1.or.j==ny) p = p / 2
        if (k==1.or.k==nz) p = p / 2
        avg = avg + p
      end do
   end do
   call test3d([(poisfft_dirichlet, i=1, 6)], resexact3d_dir)

   dx = lx / nx; dy = ly / ny; dz = lz / nz
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            x = dx * (i - 0.5_rp); y = dy * (j - 0.5_rp); z = dz * (k - 0.5_rp)
            rhs3d(i, j, k) = sin(3 * pi * x / lx) * sin(5 * pi * y / ly) * sin(7 * pi * z / lz)
         end do
      end do
   end do
   call test3d([(poisfft_dirichletstag, i=1, 6)], resexact3d_dirstag)

  dx = Lx / nx
  dy = Ly / ny
  dz = Lz / nz
  do k = 1,nz
    do j = 1,ny
      do i = 1,nx
        x = dx*(i-0.5_rp); y = dy*(j-0.5_rp); z = dz*(k-0.5_rp)
        RHS3D(i,j,k) = cos(3*pi*x/Lx) * cos(5*pi*y/Ly) * cos(7*pi*z/Lz)
      end do
    end do
  end do
  RHS3D = RHS3D - sum(RHS3D)/(size(RHS3D,kind=size_kind))
  call Test3D([(PoisFFT_NeumannStag, i = 1,6)], ResExact3D_NeumStag)

   dx = lx / nx; dy = ly / ny; dz = lz / nz
   call test3d([(poisfft_periodic, i=1, 4), (poisfft_neumannstag, i=5, 6)])

  dx = Lx / nx
  dy = Ly / ny
  dz = Lz / nz
  call Test3D([(PoisFFT_NeumannStag, i = 1,2),(PoisFFT_PERIODIC, i = 3,6)])

  dx = Lx / nx
  dy = Ly / ny
  dz = Lz / nz
  call Test3D([(PoisFFT_PERIODIC, i = 1,2),(PoisFFT_NeumannStag, i = 3,4),(PoisFFT_PERIODIC, i = 5,6)])

   dx = lx / nx; dy = ly / ny; dz = lz / nz
   call test3d([(poisfft_periodic, i=1, 2), (poisfft_neumannstag, i=3, 6)])

end program
