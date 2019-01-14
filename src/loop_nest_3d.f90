!$omp do
do k = 2, self % nz
   do j = 2, self % ny
      do i = 2, self % nx
         self % XWORK(i, j, k) = self % XWORK(i, j, k) / (self % denomx(i) + self % denomy(j) + self % denomz(k))
      end do
   end do
end do
!$omp end do

!$omp do
do j = 2, self % ny
   do i = 2, self % nx
      self % XWORK(i, j, 1) = self % XWORK(i, j, 1) / (self % denomx(i) + self % denomy(j))
   end do
end do
!$omp end do

!$omp do
do k = 2, self % nz
   do i = 2, self % nx
      self % XWORK(i, 1, k) = self % XWORK(i, 1, k) / (self % denomx(i) + self % denomz(k))
   end do
end do
!$omp end do

!$omp do
do k = 2, self % nz
   do j = 2, self % ny
      self % XWORK(1, j, k) = self % XWORK(1, j, k) / (self % denomy(j) + self % denomz(k))
   end do
end do
!$omp end do

!$omp do
do i = 2, self % nx
   self % XWORK(i, 1, 1) = self % XWORK(i, 1, 1) / (self % denomx(i))
end do
!$omp end do

!$omp do
do j = 2, self % ny
   self % XWORK(1, j, 1) = self % XWORK(1, j, 1) / (self % denomy(j))
end do
!$omp end do

!$omp do
do k = 2, self % nz
   self % XWORK(1, 1, k) = self % XWORK(1, 1, k) / (self % denomz(k))
end do
!$omp end do


