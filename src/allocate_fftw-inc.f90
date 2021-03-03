#if (DIM==1)
#define POISFFT_SOLVERXD poisfft_solver1d
#define GDIMS [self%gnx]
#elif (DIM==2)
#define POISFFT_SOLVERXD poisfft_solver2d
#define GDIMS [self%gny,self%gnx]
#else
#define POISFFT_SOLVERXD poisfft_solver3d
#define GDIMS [self%gnz,self%gny,self%gnx]
#endif

#if (ISREAL==1)
#define DATATYPE real(RP)
#define DATAF    real(1,RP)
#define XWORK    rwork
#else
#define DATATYPE complex(CP)
#define DATAF    cmplx(1,1,CP)
#define XWORK    cwork
#endif

type(POISFFT_SOLVERXD), intent(inout) :: self
type(c_ptr) :: p
#if defined(MPI) && (DIM>1)
integer(c_size_t) :: cnt
integer(c_intptr_t) :: a1(DIM), a2(DIM), a3(DIM), a4(DIM)
!                   local_ni local_i_start local_no local_o_start

      class( PoisFFT_SolverXD ), intent(inout)        :: D
      type(c_ptr) :: p
#if defined(MPI) && dimensions>1
      integer(c_size_t) :: cnt
      integer(c_intptr_t) :: a1(dimensions),a2(dimensions),a3(dimensions),a4(dimensions)

      cnt =  pfft_local_size_dft(dimensions,int(gdims,c_intptr_t), &
                  D%mpi%comm, PFFT_TRANSPOSED_NONE, &
                  a1,a2,a3,a4)
      if (.not.all(a1(dimensions:1:-1)==D%nxyz)) then
        write(*,*) "PoisFFT Error. Inconsistent size of local arrays from pfft_local_size_dft!"
        write(*,*) a1(dimensions:1:-1)
        write(*,*) D%nxyz
        stop
      end if  
      p = fftw_malloc( storage_size( dataf )/storage_size('a') * cnt )
#elif defined(MPI) && dimensions==1
      p = fftw_malloc( storage_size( dataf )/storage_size('a') * int(D%gnx, c_size_t ))
#else
p = fftw_malloc(storage_size(DATAF) / storage_size('_') * self % cnt)
#endif

if (c_associated(p)) then
#if defined(MPI) && (DIM==1)
   call c_f_pointer(p, self % XWORK, GDIMS)
#else
   call c_f_pointer(p, self % XWORK, self % nxyz)
#endif
else
   stop 'Data allocate error, fftw_malloc returned NULL.'
endif

self % XWORK = 0
#undef XWORK
#undef DATATYPE
#undef DATAF
#undef GDIMS
#undef POISFFT_SOLVERXD
