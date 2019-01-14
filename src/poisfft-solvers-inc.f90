subroutine poisfft_solver1d_fullperiodic(self, phi, rhs)
   type(poisfft_solver1d), intent(inout) :: self
   real(RP), intent(out) :: phi(:)
   real(RP), intent(in) :: rhs(:)
   integer i
#ifdef MPI
   interface
      subroutine mpi_gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, ierror)
         integer :: sendbuf, recvbuf(*)
         integer :: sendcount, sendtype, recvcount, recvtype, root
         integer :: comm, ierror
      end subroutine

      subroutine mpi_gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm, ierror)
         import
         real(RP) :: sendbuf(*), recvbuf(*)
         integer :: sendcount, sendtype, recvcounts(*), displs(*)
         integer :: recvtype, root, comm, ierror
      end subroutine

      subroutine mpi_scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, ierror)
         import
         real(RP) :: sendbuf(*), recvbuf(*)
         integer :: sendcounts(*), displs(*), sendtype
         integer :: recvcount, recvtype, root, comm, ierror
      end subroutine
   end interface

   integer, allocatable :: displs(:), counts(:)
   real(RP), allocatable :: tmp(:)
   integer :: ie

   if (self % mpi_transpose_needed) then
      call mpi_comm_rank(self % mpi % comm, self % mpi % rank, ie)
      call mpi_comm_size(self % mpi % comm, self % mpi % np, ie)
      allocate(displs(self % mpi % np))
      allocate(counts(self % mpi % np))
      if (self % mpi % rank == 0) then; allocate(tmp(1:self % gnx)); else; allocate(tmp(1)); end if

      call mpi_gather(self % offx, 1, MPI_INTEGER, displs, 1, MPI_INTEGER, 0, self % mpi % comm, ie)
      call mpi_gather(self % nx, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, self % mpi % comm, ie)
      call mpi_gatherv(rhs, self % nx, MPI_RP, tmp, counts, displs, MPI_RP, 0, self % mpi % comm, ie)
      if (ie /= 0) stop "Error in MPI_Gatherv!"

      if (self % mpi % rank == 0) self % cwork = cmplx(tmp, real(0., RP), CP)

      if (self % mpi % rank == 0) then
         call execute(self % forward, self % cwork)
         if (self % offx == 0) self % cwork(1) = 0
         forall(i=2:self % gnx) self % cwork(i) = self % cwork(i) / self % denomx(i)
         call execute(self % backward, self % cwork)
      end if

      tmp = real(self % cwork, RP) / self % norm_factor
      call mpi_scatterv(tmp, counts, displs, MPI_RP, phi, self % nx, MPI_RP, 0, self % mpi % comm, ie)
   end if
#else
   ! Forward FFT of RHS
   self % cwork = cmplx(rhs, real(0., RP), CP)
   call execute(self % forward, self % cwork)
   forall(i=2:self % nx) self % cwork(i) = self % cwork(i) / self % denomx(i)
   call execute(self % backward, self % cwork)
   phi = real(self % cwork, RP) / self % norm_factor
#endif
end subroutine

subroutine poisfft_solver1d_fulldirichlet(self, phi, rhs)
   type(poisfft_solver1d), intent(inout) :: self
   real(RP), intent(out) :: phi(:)
   real(RP), intent(in) :: rhs(:)
   integer :: i

   ! forward fft of rhs
   self % rwork = rhs
   call execute(self % forward, self % rwork)
   forall(i=1:self % nx) self % rwork(i) = self % rwork(i) / self % denomx(i)
   call execute(self % backward, self % rwork)
   phi = self % rwork / self % norm_factor
end subroutine

subroutine poisfft_solver1d_fullneumann(self, phi, rhs)
   type(poisfft_solver1d), intent(inout) :: self
   real(RP), intent(out) :: phi(:)
   real(RP), intent(in) :: rhs(:)
   integer :: i

   ! forward fft of rhs
   self % rwork = rhs
   call execute(self % forward, self % rwork)
   self % rwork(1) = 0
   forall(i=2:self % nx) self % rwork(i) = self % rwork(i) / self % denomx(i)
   call execute(self % backward, self % rwork)
   phi = self % rwork / self % norm_factor
end subroutine

subroutine poisfft_solver2d_fullperiodic(self, phi, rhs)
   type(poisfft_solver2d), intent(inout) :: self
   real(RP), intent(out) :: phi(:, :)
   real(RP), intent(in) :: rhs(:, :)
   integer :: i, j

   ! Forward FFT of RHS
   self % cwork = cmplx(RHS, real(0., RP), CP)

#ifdef MPI
   call execute_mpi(self % forward)
#else
   call execute(self % forward, self % cwork)
#endif

   if (self % offx == 0 .and. self % offy == 0) self % cwork(1, 1) = 0
   do j = 1, self % ny
      do i = max(3 - j - self % offx - self % offy, 1), self % nx
         self % cwork(i, j) = self % cwork(i, j) / (self % denomx(i) + self % denomy(j))
      end do
   end do

#ifdef MPI
   call execute_mpi(self % backward)
#else
   call execute(self % backward, self % cwork)
#endif
   phi = real(self % cwork, RP) / self % norm_factor
end subroutine

subroutine poisfft_solver2d_fulldirichlet(self, phi, rhs)
   type(poisfft_solver2d), intent(inout) :: self
   real(RP), intent(out) :: phi(:, :)
   real(RP), intent(in) :: rhs(:, :)
   integer :: i, j

   ! forward fft of rhs
   self % rwork = rhs
   call execute(self % forward, self % rwork)
   do j = 1, self % ny
      do i = 1, self % nx
         self % rwork(i, j) = self % rwork(i, j) / (self % denomx(i) + self % denomy(j))
      end do
   end do

   call execute(self % backward, self % rwork)
   phi = self % rwork / self % norm_factor
end subroutine

subroutine poisfft_solver2d_fullneumann(self, phi, rhs)
   type(poisfft_solver2d), intent(inout) :: self
   real(RP), intent(out) :: phi(:, :)
   real(RP), intent(in) :: rhs(:, :)
   integer :: i, j

   ! forward fft of rhs
   self % rwork = rhs
   call execute(self % forward, self % rwork)
   self % rwork(1, 1) = 0
   do j = 1, self % ny
      do i = max(3 - j, 1), self % nx
         self % rwork(i, j) = self % rwork(i, j) / (self % denomx(i) + self % denomy(j))
      end do
   end do
   call execute(self % backward, self % rwork)
   phi = self % rwork / self % norm_factor
end subroutine

subroutine poisfft_solver3d_fullperiodic(self, phi, rhs)
   type(poisfft_solver3d), intent(inout) :: self
   real(RP), intent(out) :: phi(:, :, :)
   real(RP), intent(in) :: rhs(:, :, :)
   integer :: i, j, k

   !$omp parallel private(i,j,k)
   !$omp workshare
   self % cwork = cmplx(RHS, real(0., RP), CP)
   !$omp end workshare
   !$omp end parallel
#ifdef MPI
   call execute_mpi(self % forward)
#else
   call execute(self % forward, self % cwork)
#endif
   !$omp parallel private(i,j,k)
   if (self % offx == 0 .and. self % offy == 0 .and. self % offz == 0) then
#define xwork cwork
#include "loop_nest_3d.f90"
#undef xwork
      !NOTE: if IEEE FPE exceptions are disabled all this is not necessary and
      ! the loop can be over all indexes because the infinity or NaN is changed to 0 below

      !$omp single
      self % cwork(1, 1, 1) = 0
      !$omp end single
   else
      !$omp do
      do k = 1, self % nz; do j = 1, self % ny; do i = 1, self % nx
         self % cwork(i, j, k) = self % cwork(i, j, k) / (self % denomx(i) + self % denomy(j) + self % denomz(k))
      end do; end do; end do
      !$omp end do
   end if
   !$omp end parallel
#ifdef MPI
   call execute_mpi(self % backward)
#else
   call execute(self % backward, self % cwork)
#endif
   !$omp parallel
   !$omp workshare
   phi = real(self % cwork, RP) / self % norm_factor
   !$omp end workshare
   !$omp end parallel
end subroutine

subroutine poisfft_solver3d_fulldirichlet(self, phi, rhs)
   type(poisfft_solver3d), intent(inout) :: self
   real(RP), intent(out) :: phi(:, :, :)
   real(RP), intent(in) :: rhs(:, :, :)
   integer :: i, j, k

   !$omp parallel private(i,j,k)
   !$omp workshare
   self % rwork = RHS
   !$omp end workshare
   !$omp end parallel
#ifdef MPI
   call execute_mpi(self % forward)
#else
   call execute(self % forward, self % rwork)
#endif
   !$omp parallel private(i,j,k)
   !$omp do
   do k = 1, self % nz
      do j = 1, self % ny
         do i = 1, self % nx
            self % rwork(i, j, k) = self % rwork(i, j, k) / (self % denomx(i) + self % denomy(j) + self % denomz(k))
         end do
      end do
   end do
   !$omp end do
   !$omp end parallel
#ifdef MPI
   call execute_mpi(self % backward)
#else
   call execute(self % backward, self % rwork)
#endif
   !$omp parallel private(i,j,k)
   !$omp workshare
   phi = self % rwork / self % norm_factor
   !$omp end workshare
   !$omp end parallel
end subroutine

subroutine poisfft_solver3d_fullneumann(self, phi, rhs)
   type(poisfft_solver3d), intent(inout) :: self
   real(RP), intent(out) :: phi(:, :, :)
   real(RP), intent(in) :: rhs(:, :, :)
   integer :: i, j, k

   !$omp parallel private(i,j,k)

   !$omp workshare
   self % rwork = RHS
   !$omp end workshare

   !$omp end parallel
#ifdef MPI
   call execute_mpi(self % forward)
#else
   call execute(self % forward, self % rwork)
#endif
   !$omp parallel private(i,j,k)
   !$omp single
   if (self % offx == 0 .and. self % offy == 0 .and. self % offz == 0) self % rwork(1, 1, 1) = 0
   !$omp end single

   if (self % offx == 0 .and. self % offy == 0 .and. self % offz == 0) then
#define xwork rwork
#include "loop_nest_3d.f90"
#undef xwork
      !NOTE: if IEEE FPE exceptions are disabled all this is not necessary and
      ! the loop can be over all indexes because the infinity or NaN is changed to 0 below

      !$omp single
      self % rwork(1, 1, 1) = 0
      !$omp end single
   else
      !$omp do
      do k = 1, self % nz
         do j = 1, self % ny
            do i = 1, self % nx
               self % rwork(i, j, k) = self % rwork(i, j, k) / (self % denomx(i) + self % denomy(j) + self % denomz(k))
            end do
         end do
      end do
      !$omp end do
   end if
   !$omp end parallel
#ifdef MPI
   call execute_mpi(self % backward)
#else
   call execute(self % backward, self % rwork)
#endif
   !$omp parallel private(i,j,k)
   !$omp workshare
   phi = self % rwork / self % norm_factor
   !$omp end workshare
   !$omp end parallel
end subroutine

#ifdef MPI
subroutine poisfft_solver3d_ppns(self, phi, rhs)
   type(poisfft_solver3d), intent(inout) :: self
   real(RP), intent(out) :: phi(:, :, :)
   real(RP), intent(in) :: rhs(:, :, :)
   integer :: i, j, k, tid

   if (self % solvers1d(3) % mpi_transpose_needed) then
      call transform_1d_real_z(self % solvers1d(3:), phi, rhs, forward=.true., use_rhs=.true.)
   else
      !$omp parallel private(tid,i,j,k)
      tid = 0
      !$ tid = omp_get_thread_num()
      !$omp do collapse(2)
      do j = 1, self % ny
         do i = 1, self % nx
            self % solvers1d(3 + tid) % rwork = rhs(i, j, :)
            call execute(self % solvers1d(3 + tid) % forward, self % solvers1d(3 + tid) % rwork)
            phi(i, j, :) = self % solvers1d(3 + tid) % rwork
         end do
      end do
      !$omp end parallel
   end if

   if (self % ny < self % gny) then
      do k = 1, self % nz
         !$omp parallel workshare
         self % solvers2d(1) % cwork(1:self % nx, 1:self % ny) = cmplx(phi(1:self % nx, 1:self % ny, k), real(0., RP), CP)
         !$omp end parallel workshare

         call execute_mpi(self % solvers2d(1) % forward, self % solvers2d(1) % cwork)

         if (k == 1 .and. self % offz == 0) then
            !$omp parallel do
            do j = 1, self % ny
               do i = max(3 - j - self % offx - self % offy, 1), self % nx
                  self % solvers2d(1) % cwork(i, j) = self % solvers2d(1) % cwork(i, j) / (self % denomx(i) + self % denomy(j))
               end do
            end do
            !$omp end parallel do
            !NOTE: if IEEE FPE exceptions are disabled all this is not necessary and
            ! the loop can be over all indexes because the infinity or NaN is changed to 0 below
            if (self % offx == 0 .and. self % offy == 0) self % solvers2d(1) % cwork(1, 1) = 0
         else
            !$omp parallel do collapse(2)
            do j = 1, self % ny
               do i = 1, self % nx
                  self % solvers2d(1) % cwork(i, j) = self % solvers2d(1) % cwork(i, j) / (self % denomx(i) + self % denomy(j) + self % denomz(k))
               end do
            end do
            !$omp end parallel do
         endif
         call execute_mpi(self % solvers2d(1) % backward, self % solvers2d(1) % cwork)
         !$omp parallel workshare
         phi(:, :, k) = real(self % solvers2d(1) % cwork, RP) / self % norm_factor
         !$omp end parallel workshare
      end do
   else

      !$omp parallel private(tid,i,j,k)
      tid = 1
      !$ tid = omp_get_thread_num()+1

      !$omp do
      do k = 1, self % nz
         self % solvers2d(tid) % cwork(1:self % nx, 1:self % ny) = cmplx(phi(1:self % nx, 1:self % ny, k), real(0., RP), CP)
         call execute(self % solvers2d(tid) % forward, self % solvers2d(tid) % cwork)

         if (k == 1 .and. self % offz == 0) then
            do j = 2, self % ny
               do i = 2, self % nx
                  self % solvers2d(tid) % cwork(i, j) = self % solvers2d(tid) % cwork(i, j) / (self % denomx(i) + self % denomy(j))
               end do
            end do
            do i = 2, self % nx
               self % solvers2d(tid) % cwork(i, 1) = self % solvers2d(tid) % cwork(i, 1) / self % denomx(i)
            end do
            do j = 2, self % ny
               self % solvers2d(tid) % cwork(1, j) = self % solvers2d(tid) % cwork(1, j) / self % denomy(j)
            end do
            !NOTE: if IEEE FPE exceptions are disabled all this is not necessary and
            ! the loop can be over all indexes because the infinity or NaN is changed to 0 below
            self % solvers2d(tid) % cwork(1, 1) = 0
         else
            do j = 1, self % ny
               do i = 1, self % nx
                  self % solvers2d(tid) % cwork(i, j) = self % solvers2d(tid) % cwork(i, j) / (self % denomx(i) + self % denomy(j) + self % denomz(k))
               end do
            end do
         endif
         call execute(self % solvers2d(tid) % backward, self % solvers2d(tid) % cwork)
         phi(:, :, k) = real(self % solvers2d(tid) % cwork, RP) / self % norm_factor
      end do
      !$omp end parallel
   end if

   if (self % solvers1d(3) % mpi_transpose_needed) then
      call transform_1d_real_z(self % solvers1d(3:), phi, rhs, forward=.false., use_rhs=.false.)
   else
      !$omp parallel private(tid,i,j,k)
      tid = 0
      !$ tid = omp_get_thread_num()
      !$omp do collapse(2)
      do j = 1, self % ny
         do i = 1, self % nx
            self % solvers1d(3 + tid) % rwork = phi(i, j, :)
            call execute(self % solvers1d(3 + tid) % backward, self % solvers1d(3 + tid) % rwork)
            phi(i, j, :) = self % solvers1d(3 + tid) % rwork
         end do
      end do
      !$omp end parallel
   end if
end subroutine

#else

subroutine poisfft_solver3d_ppns(self, phi, rhs)
   type(poisfft_solver3d), intent(inout) :: self
   real(RP), intent(out) :: phi(:, :, :)
   real(RP), intent(in) :: rhs(:, :, :)
   integer :: i, j, k, tid !thread id

   !$omp parallel private(tid,i,j,k)
   tid = 1
   !$ tid = omp_get_thread_num()+1

   ! Forward FFT of RHS in z dimension according to Wilhelmson, Ericksen, JCP 1977
   !$omp do
   do j = 1, self % ny
      do i = 1, self % nx
         self % solvers1d(tid) % rwork = rhs(i, j, :)
         call execute(self % solvers1d(tid) % forward, self % solvers1d(tid) % rwork)
         phi(i, j, :) = self % solvers1d(tid) % rwork
      end do
   end do
   !$omp end do

   !$omp do
   do k = 1, self % nz
      self % solvers2d(tid) % cwork(1:self % nx, 1:self % ny) = cmplx(phi(1:self % nx, 1:self % ny, k), real(0., RP), CP)

      call execute(self % solvers2d(tid) % forward, self % solvers2d(tid) % cwork)
      if (k == 1) then
         do j = 2, self % ny
            do i = 2, self % nx
               self % solvers2d(tid) % cwork(i, j) = self % solvers2d(tid) % cwork(i, j) / (self % denomx(i) + self % denomy(j))
            end do
         end do
         do i = 2, self % nx
            self % solvers2d(tid) % cwork(i, 1) = self % solvers2d(tid) % cwork(i, 1) / self % denomx(i)
         end do
         do j = 2, self % ny
            self % solvers2d(tid) % cwork(1, j) = self % solvers2d(tid) % cwork(1, j) / self % denomy(j)
         end do
         !NOTE: if IEEE FPE exceptions are disabled all this is not necessary and
         ! the loop can be over all indexes because the infinity or NaN is changed to 0 below
         self % solvers2d(tid) % cwork(1, 1) = 0
      else
         do j = 1, self % ny
            do i = 1, self % nx
               self % solvers2d(tid) % cwork(i, j) = self % solvers2d(tid) % cwork(i, j) / (self % denomx(i) + self % denomy(j) + self % denomz(k))
            end do
         end do
      endif
      call execute(self % solvers2d(tid) % backward, self % solvers2d(tid) % cwork)
      phi(:, :, k) = real(self % solvers2d(tid) % cwork, RP) / self % norm_factor
   end do
   !$omp end do

   !$omp do
   do j = 1, self % ny
      do i = 1, self % nx
         self % solvers1d(tid) % rwork = phi(i, j, :)
         call execute(self % solvers1d(tid) % backward, self % solvers1d(tid) % rwork)
         phi(i, j, :) = self % solvers1d(tid) % rwork
      end do
   end do
   !$omp end do
   !$omp end parallel
end subroutine

#endif

subroutine poisfft_solver3d_pnsp(self, phi, rhs)
   type(poisfft_solver3d), intent(inout) :: self
   real(RP), intent(out) :: phi(:, :, :)
   real(RP), intent(in) :: rhs(:, :, :)
   integer :: i, j, k, tid !thread id

   !$omp parallel private(tid,i,j,k)
   tid = 1
   !$ tid = omp_get_thread_num()+1

   ! Forward FFT of RHS in y dimension
   !$omp do
   do k = 1, self % nz
      do i = 1, self % nx
         self % solvers1d(tid) % rwork = rhs(i, :, k)
         call execute(self % solvers1d(tid) % forward, self % solvers1d(tid) % rwork)
         phi(i, :, k) = self % solvers1d(tid) % rwork
      end do
   end do
   !$omp end do

   !$omp do
   do j = 1, self % ny
      self % solvers2d(tid) % cwork(1:self % nx, 1:self % nz) = cmplx(phi(1:self % nx, 1:self % nz, k), real(0., RP), CP)
      call execute(self % solvers2d(tid) % forward, self % solvers2d(tid) % cwork)
      if (j == 1) then
         do k = 2, self % nz
            do i = 2, self % nx
               self % solvers2d(tid) % cwork(i, k) = self % solvers2d(tid) % cwork(i, k) / (self % denomx(i) + self % denomz(k))
            end do
         end do
         do i = 2, self % nx
            self % solvers2d(tid) % cwork(i, 1) = self % solvers2d(tid) % cwork(i, 1) / self % denomx(i)
         end do
         do k = 2, self % nz
            self % solvers2d(tid) % cwork(1, k) = self % solvers2d(tid) % cwork(1, k) / self % denomz(k)
         end do
         !NOTE: if IEEE FPE exceptions are disabled all this is not necessary and
         ! the loop can be over all indexes because the infinity or NaN is changed to 0 below
         self % solvers2d(tid) % cwork(1, 1) = 0
      else
         do k = 1, self % nz
            do i = 1, self % nx
               self % solvers2d(tid) % cwork(i, k) = self % solvers2d(tid) % cwork(i, k) / (self % denomx(i) + self % denomy(j) + self % denomz(k))
            end do
         end do
      endif
      call execute(self % solvers2d(tid) % backward, self % solvers2d(tid) % cwork)
      phi(:, :, k) = real(self % solvers2d(tid) % cwork, RP) / self % norm_factor
   end do
   !$omp end do

   !$omp do
   do k = 1, self % nz
      do i = 1, self % nx
         self % solvers1d(tid) % rwork = phi(i, :, k)
         call execute(self % solvers1d(tid) % backward, self % solvers1d(tid) % rwork)
         phi(i, :, k) = self % solvers1d(tid) % rwork
      end do
   end do
   !$omp end do
   !$omp end parallel
end subroutine

#ifdef MPI
subroutine poisfft_solver3d_pnsns(self, phi, rhs)
   type(poisfft_solver3d), intent(inout), target :: self
   real(RP), intent(out) :: phi(:, :, :)
   real(RP), intent(in) :: rhs(:, :, :)
   integer :: i, j, k, tid

   if (self % solvers1d(3) % mpi_transpose_needed) then
      call transform_1d_real_z(self % solvers1d(3 :: 3), phi, rhs, forward=.true., use_rhs=.true.)
   else
      !$omp parallel private(tid,i,j,k)
      tid = 0
      !$ tid = omp_get_thread_num()
      !$omp do
      do j = 1, self % ny
         do i = 1, self % nx
            self % solvers1d(3 + 3 * tid) % rwork = rhs(i, j, :)
            call execute(self % solvers1d(3 + 3 * tid) % forward, self % solvers1d(3 + 3 * tid) % rwork)
            phi(i, j, :) = self % solvers1d(3 + 3 * tid) % rwork
         end do
      end do
      !$omp end parallel
   end if

   if (self % solvers1d(2) % mpi_transpose_needed) then
      call transform_1d_real_y(self % solvers1d(2 :: 3), phi, rhs, forward=.true., use_rhs=.false.)
   else
      !$omp parallel private(tid,i,j,k)
      tid = 0
      !$ tid = omp_get_thread_num()
      !$omp do
      do k = 1, self % nz
         do i = 1, self % nx
            self % solvers1d(2 + 3 * tid) % rwork = phi(i, :, k)
            call execute(self % solvers1d(2 + 3 * tid) % forward, self % solvers1d(2 + 3 * tid) % rwork)
            phi(i, :, k) = self % solvers1d(2 + 3 * tid) % rwork
         end do
      end do
      !$omp end parallel
   end if

   !$omp parallel private(tid,i,j,k)
   tid = 0
   !$ tid = omp_get_thread_num()
   !$omp do
   do k = 1, self % nz
      do j = 1, self % ny
         self % solvers1d(1 + 3 * tid) % cwork = cmplx(phi(:, j, k), real(0., RP), CP)
         call execute(self % solvers1d(1 + 3 * tid) % forward, self % solvers1d(1 + 3 * tid) % cwork)

         do i = max(4 - j - k - self % offy - self % offz, 1), self % nx
            self % solvers1d(1 + 3 * tid) % cwork(i) = self % solvers1d(1 + 3 * tid) % cwork(i) / (self % denomx(i) + self % denomy(j) + self % denomz(k))
         end do
         !note: if ieee fpe exceptions are disabled all this is not necessary and
         ! the loop can be over all indexes because the infinity or nan is changed to 0 below
         if (self % offy == 0 .and. self % offz == 0 .and. j == 1 .and. k == 1) self % solvers1d(1 + 3 * tid) % cwork(1) = 0

         call execute(self % solvers1d(1 + 3 * tid) % backward, self % solvers1d(1 + 3 * tid) % cwork)
         phi(:, j, k) = real(self % solvers1d(1 + 3 * tid) % cwork, RP) / self % norm_factor
      end do
   end do
   !$omp end parallel

   if (self % Solvers1D(2) % mpi_transpose_needed) then
      call transform_1d_real_y(self % Solvers1D(2 :: 3), Phi, RHS, forward=.false., use_rhs=.false.)
   else
      !$omp parallel private(tid,i,j,k)
      tid = 0
      !$ tid = omp_get_thread_num()
      !$omp do
      do k = 1, self % nz
         do i = 1, self % nx
            self % Solvers1D(2 + 3 * tid) % rwork = Phi(i, :, k)
            call Execute(self % Solvers1D(2 + 3 * tid) % backward, self % Solvers1D(2 + 3 * tid) % rwork)
            Phi(i, :, k) = self % Solvers1D(2 + 3 * tid) % rwork
         end do
      end do
      !$omp end parallel
   end if

   if (self % solvers1d(3) % mpi_transpose_needed) then
      call transform_1d_real_z(self % solvers1d(3 :: 3), phi, rhs, forward=.false., use_rhs=.false.)
   else
      !$omp parallel private(tid,i,j,k)
      tid = 0
      !$ tid = omp_get_thread_num()
      !$omp do
      do j = 1, self % ny
         do i = 1, self % nx
            self % solvers1d(3 + 3 * tid) % rwork = phi(i, j, :)
            call execute(self % solvers1d(3 + 3 * tid) % backward, self % solvers1d(3 + 3 * tid) % rwork)
            phi(i, j, :) = self % solvers1d(3 + 3 * tid) % rwork
         end do
      end do
      !$omp end parallel
   end if
end subroutine
 !mpi
#else
 !threads
subroutine poisfft_solver3d_pnsns(self, phi, rhs)
   type(poisfft_solver3d), intent(inout) :: self
   real(RP), intent(out) :: phi(:, :, :)
   real(RP), intent(in) :: rhs(:, :, :)
   integer :: i, j, k, tid !thread id

   !$omp parallel private(tid,i,j,k)
   tid = 1
   !$ tid = omp_get_thread_num()+1

   ! forward fft of rhs in x dimension
   !$omp do
   do i = 1, self % nx
      self % solvers2d(tid) % rwork = rhs(i, :, :)
      call execute(self % solvers2d(tid) % forward, self % solvers2d(tid) % rwork)
      phi(i, :, :) = self % solvers2d(tid) % rwork
   end do
   !$omp end do

   !$omp do
   do k = 1, self % nz
      do j = 1, self % ny
         self % solvers1d(tid) % cwork = cmplx(phi(:, j, k), real(0., RP), CP)
         call execute(self % solvers1d(tid) % forward, self % solvers1d(tid) % cwork)
         do i = max(4 - j - k - self % offy - self % offz, 1), self % nx
            self % solvers1d(tid) % cwork(i) = self % solvers1d(tid) % cwork(i) / (self % denomx(i) + self % denomy(j) + self % denomz(k))
         end do
         !note: if ieee fpe exceptions are disabled all this is not necessary and
         ! the loop can be over all indexes because the infinity or nan is changed to 0 below
         if (self % offy == 0 .and. self % offz == 0 .and. j == 1 .and. k == 1) self % solvers1d(tid) % cwork(1) = 0
         call execute(self % solvers1d(tid) % backward, self % solvers1d(tid) % cwork)
         phi(:, j, k) = real(self % solvers1d(tid) % cwork, RP) / self % norm_factor
      end do
   end do
   !$omp end do

   !$omp do
   do i = 1, self % nx
      self % solvers2d(tid) % rwork = phi(i, :, :)
      call execute(self % solvers2d(tid) % backward, self % solvers2d(tid) % rwork)
      phi(i, :, :) = self % solvers2d(tid) % rwork
   end do
   !$omp end do
   !$omp end parallel

end subroutine

#endif

#ifdef MPI
subroutine transform_1d_real_y(D1D, Phi, RHS, forward, use_rhs)
   type(PoisFFT_Solver1D), intent(inout), target :: D1D(:)
   real(RP), intent(inout) :: Phi(:, :, :)
   real(RP), intent(in) :: RHS(:, :, :)
   logical, intent(in) :: forward, use_rhs
   integer :: nx, ny, nz, i, j, k, l, ie, tid
   interface
      subroutine mpi_alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm, ierror)
         import
         real(RP) :: sendbuf(*), recvbuf(*)
         integer :: sendcounts(*), sdispls(*), sendtype
         integer :: recvcounts(*), rdispls(*), recvtype
         integer :: comm, ierror
      end subroutine
   end interface

   nx = size(Phi, 1)
   ny = size(Phi, 2)
   nz = size(Phi, 3)

#define m D1D(1)%mpi

   !$omp parallel private(i,j,k,tid)
   tid = 1
   !$ tid = omp_get_thread_num()+1

   !step1 local transpose
   if (use_rhs) then
      !$omp do
      do k = 1, nz; do i = 1, nx; do j = 1, ny
         m % tmp1(j, k, i) = RHS(i, j, k)
      end do; end do; end do
   else
      !$omp do
      do k = 1, nz; do i = 1, nx; do j = 1, ny
         m % tmp1(j, k, i) = Phi(i, j, k)
      end do; end do; end do
   end if

   !step2 exchange
   !$omp single
   call mpi_alltoallv(m % tmp1, m % scounts, m % sdispls, MPI_RP, m % tmp2, m % rcounts, m % rdispls, MPI_RP, m % comm, ie)
   !$omp end single

   !step3 local reordering of blocks
   !$omp do collapse(3)
   do l = 1, m % np
      do k = 0, m % rnxs(1) - 1
         do j = 0, nz - 1
            do i = 0, m % rnzs(l) - 1
               m % rwork(i + m % sumrnzs(l) + 1, j + 1, k + 1) = &
                  m % tmp2(i + j * m % rnzs(l) + k * (nz * m % rnzs(l)) + m % rdispls(l))
            end do
         end do
      end do
   end do

   ! Forward FFT of RHS in z dimension according to Wilhelmson, Ericksen, JCP 1977
   if (forward) then
      !$omp do
      do k = 1, size(m % rwork, 3)
         do j = 1, size(m % rwork, 2)
            !todo: is this really alligned?
            call execute(d1d(tid) % forward, m % rwork(:, j, k))
         end do
      end do
   else
      !$omp do
      do k = 1, size(m % rwork, 3)
         do j = 1, size(m % rwork, 2)
            call execute(d1d(tid) % backward, m % rwork(:, j, k))
         end do
      end do
   end if

   !step3' local reordering of blocks
   !$omp do collapse(3)
   do l = 1, m % np
      do k = 0, m % rnxs(1) - 1
         do j = 0, nz - 1
            do i = 0, m % rnzs(l) - 1
               m % tmp2(i + j * m % rnzs(l) + k * (nz * m % rnzs(l)) + m % rdispls(l)) = &
                  m % rwork(i + m % sumrnzs(l) + 1, j + 1, k + 1)
            end do
         end do
      end do
   end do

   !step2' exchange
   !$omp single
   call mpi_alltoallv(m % tmp2, m % rcounts, m % rdispls, MPI_RP, m % tmp1, m % scounts, m % sdispls, MPI_RP, m % comm, ie)
   !$omp end single

   !$omp do
   do k = 1, nz; do i = 1, nx; do j = 1, ny
      Phi(i, j, k) = m % tmp1(j, k, i)
   end do; end do; end do
   !$omp end parallel

#undef m
end subroutine

subroutine transform_1d_real_z(d1d, phi, rhs, forward, use_rhs)
   type(poisfft_solver1d), intent(inout), target :: d1d(:)
   real(RP), intent(inout) :: phi(:, :, :)
   real(RP), intent(in) :: rhs(:, :, :)
   logical, intent(in) :: forward, use_rhs
   integer :: nx, ny, nz, i, j, k, l, ie, tid
   interface
      subroutine mpi_alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm, ierror)
         import
         real(RP) :: sendbuf(*), recvbuf(*)
         integer :: sendcounts(*), sdispls(*), sendtype
         integer :: recvcounts(*), rdispls(*), recvtype
         integer :: comm, ierror
      end subroutine
   end interface

   nx = size(phi, 1)
   ny = size(phi, 2)
   nz = size(phi, 3)

#define m d1d(1) % mpi

   !$omp parallel private(i,j,k,tid)
   tid = 1
   !$ tid = omp_get_thread_num()+1

   !step1 local transpose
   if (use_rhs) then
      !$omp do
      do j = 1, ny; do i = 1, nx; do k = 1, nz
         m % tmp1(k, j, i) = rhs(i, j, k)
      end do; end do; end do
   else
      !$omp do
      do j = 1, ny; do i = 1, nx; do k = 1, nz
         m % tmp1(k, j, i) = phi(i, j, k)
      end do; end do; end do

   end if

   !step2 exchange
   !$omp single
   call mpi_alltoallv(m % tmp1, m % scounts, m % sdispls, MPI_RP, m % tmp2, m % rcounts, m % rdispls, MPI_RP, m % comm, ie)
   !$omp end single

   !step3 local reordering of blocks
   !$omp do collapse(3)
   do l = 1, m % np
      do k = 0, m % rnxs(1) - 1
         do j = 0, ny - 1
            do i = 0, m % rnzs(l) - 1
               m % rwork(i + m % sumrnzs(l) + 1, j + 1, k + 1) = &
                  m % tmp2(i + j * m % rnzs(l) + k * (ny * m % rnzs(l)) + m % rdispls(l))
            end do
         end do
      end do
   end do

   ! forward fft of rhs in z dimension according to wilhelmson, ericksen, jcp 1977
   if (forward) then
      !$omp do
      do k = 1, size(m % rwork, 3)
         do j = 1, size(m % rwork, 2)
            call execute(d1d(tid) % forward, m % rwork(:, j, k))
         end do
      end do
   else
      !$omp do
      do k = 1, size(m % rwork, 3)
         do j = 1, size(m % rwork, 2)
            call execute(d1d(tid) % backward, m % rwork(:, j, k))
         end do
      end do
   end if

   !step3' local reordering of blocks
   !$omp do collapse(3)
   do l = 1, m % np
      do k = 0, m % rnxs(1) - 1
         do j = 0, ny - 1
            do i = 0, m % rnzs(l) - 1
               m % tmp2(i + j * m % rnzs(l) + k * (ny * m % rnzs(l)) + m % rdispls(l)) = m % rwork(i + m % sumrnzs(l) + 1, j + 1, k + 1)
            end do
         end do
      end do
   end do

   !step2' exchange
   !$omp single
   call mpi_alltoallv(m % tmp2, m % rcounts, m % rdispls, MPI_RP, m % tmp1, m % scounts, m % sdispls, MPI_RP, m % comm, ie)
   !$omp end single

   !$omp do
   do j = 1, ny; do k = 1, nz; do i = 1, nx
      phi(i, j, k) = m % tmp1(k, j, i)
   end do; end do; end do
   !$omp end parallel
#undef m
end subroutine

#endif

