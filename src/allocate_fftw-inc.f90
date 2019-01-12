#if dimensions == 1
#define poisfft_solverxd poisfft_solver1d
#define gdims [self%gnx]
#elif dimensions == 2
#define poisfft_solverxd poisfft_solver2d
#define gdims [self%gny,self%gnx]
#else
#define poisfft_solverxd poisfft_solver3d
#define gdims [self%gnz,self%gny,self%gnx]
#endif

#if realcomplex == 1
#define datatype real(RP)
#define dataf    real(1,RP)
#define xwork    rwork
#else
#define datatype complex(CP)
#define dataf    cmplx(1,1,CP)
#define xwork    cwork
#endif

type(poisfft_solverxd), intent(inout) :: self
type(c_ptr) :: p
#if defined(MPI) && dimensions>1
integer(c_size_t) :: cnt
integer(c_intptr_t) :: a1(dimensions), a2(dimensions), a3(dimensions), a4(dimensions)
!                      local_ni        local_i_start,  local_no,       local_o_start

! alloc_local = pfft_local_size_dft_3d(n, comm_cart_3d, pfft_transposed_none, local_ni, local_i_start, local_no, local_o_start)
cnt = pfft_local_size_dft(dimensions, int(gdims, c_intptr_t), self % mpi % comm, pfft_transposed_none, a1, a2, a3, a4)
if (.not. all(a1(dimensions:1:-1) == self % nxyz)) then
   print *, a1(dimensions:1:-1), '!=', self%nxyz
   stop 'Error. Inconsistent size of local arrays!'
end if
p = fftw_malloc(storage_size(dataf) / storage_size('a') * cnt)
#elif defined(MPI) && dimensions==1
p = fftw_malloc(storage_size(dataf) / storage_size('a') * int(self % gnx, c_size_t))
#else
p = fftw_malloc(storage_size(dataf) / storage_size('a') * self % cnt)
#endif

if (c_associated(p)) then
#if defined(MPI) && dimensions==1
   call c_f_pointer(p, self % xwork, gdims)
#else
   call c_f_pointer(p, self % xwork, self % nxyz)
#endif
else
   stop 'Data allocate error, fftw_malloc returned NULL.'
endif

self % xwork = 0
#undef xwork
#undef datatype
#undef dataf
#undef gdims
#undef poisfft_solverxd
