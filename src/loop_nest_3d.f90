!$omp do
do k = 2, self % nz
   do j = 2, self % ny
      do i = 2, self % nx
         self % xwork(i, j, k) = self % xwork(i, j, k) / (self % denomx(i) + self % denomy(j) + self % denomz(k))
      end do
   end do
end do
!$omp end do

!$omp do
do j = 2, self % ny
   do i = 2, self % nx
      self % xwork(i, j, 1) = self % xwork(i, j, 1) / (self % denomx(i) + self % denomy(j))
   end do
end do
!$omp end do

!$omp do
do k = 2, self % nz
   do i = 2, self % nx
      self % xwork(i, 1, k) = self % xwork(i, 1, k) / (self % denomx(i) + self % denomz(k))
   end do
end do
!$omp end do

!$omp do
do k = 2, self % nz
   do j = 2, self % ny
      self % xwork(1, j, k) = self % xwork(1, j, k) / (self % denomy(j) + self % denomz(k))
   end do
end do
!$omp end do

!$omp do
do i = 2, self % nx
   self % xwork(i, 1, 1) = self % xwork(i, 1, 1) / (self % denomx(i))
end do
!$omp end do

!$omp do
do j = 2, self % ny
   self % xwork(1, j, 1) = self % xwork(1, j, 1) / (self % denomy(j))
end do
!$omp end do

!$omp do
do k = 2, self % nz
   self % xwork(1, 1, k) = self % xwork(1, 1, k) / (self % denomz(k))
end do
!$omp end do


