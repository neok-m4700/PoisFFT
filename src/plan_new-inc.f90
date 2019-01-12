#if (dimensions == 1)

#define POISFFT_SOLVERXD PoisFFT_Solver1D
#define POISFFT_PLANXD   PoisFFT_Plan1D
#define REALPLANTYPES    plantypes(1)
#define NXYZS            self%gnx
#define COLONS           :

#elif (dimensions == 2)

#define POISFFT_SOLVERXD PoisFFT_Solver2D
#define POISFFT_PLANXD   PoisFFT_Plan2D
#define REALPLANTYPES    plantypes(1),plantypes(2)
#define NXYZS            self%gny,self%gnx
#define COLONS           :,:

#else

#define POISFFT_SOLVERXD PoisFFT_Solver3D
#define POISFFT_PLANXD   PoisFFT_Plan3D
#define REALPLANTYPES    plantypes(1),plantypes(2),plantypes(3)
#define NXYZS            self%gnz,self%gny,self%gnx
#define COLONS           :,:,:

#endif

#if (PREC == 2)
#define PFFT_CMPLX pfft_plan_dft
#define PFFT_REAL pfft_plan_r2r
#define FFTW_CMPLX_MPI_2D fftw_mpi_plan_dft_2d
#else
#define PFFT_CMPLX pfftf_plan_dft
#define PFFT_REAL pfftf_plan_r2r
#define FFTW_CMPLX_MPI_2D fftwf_mpi_plan_dft_2d
#endif

type(POISFFT_PLANXD) :: plan

type(POISFFT_SOLVERXD), intent(inout) :: self
integer, intent(in), dimension(:) :: plantypes
logical, intent(in), optional :: distributed
logical :: distr

distr = .false.

if (plantypes(1) == fft_complex) then

   if (size(plantypes) < 2) then
      write(*, *) "Error: not enough flags when creating POISFFT_PLANXD"
      STOP
   endif

   plan % dir = plantypes(2)

#if defined(MPI) && dimensions > 1
   distr = merge(distributed, .true., present(distributed))

#if dimensions == 2
   if (distr) then
      if (self % mpi % comm_dim == 1) then
         plan % planptr = FFTW_CMPLX_MPI_2D(int(self % gny, c_intptr_t), int(self % gnx, c_intptr_t), &
            self % cwork, self % cwork, self % mpi % comm, plan % dir, fftw_measure)
         plan % method = fft_distributed_fftw
      else
         plan % planptr = PFFT_CMPLX(dimensions, int([NXYZS], c_intptr_t), self % cwork, self % cwork, self % mpi % comm, &
            plan % dir, pfft_transposed_none + pfft_measure + pfft_preserve_input)
      end if
#elif dimensions == 3
   if (distr) then
      plan % planptr = PFFT_CMPLX(dimensions, int([NXYZS], c_intptr_t), self % cwork, self % cwork, self % mpi % comm, &
         plan % dir, pfft_transposed_none + pfft_measure + pfft_preserve_input)
#endif
   else
      plan % planptr = fftw_plan_gen(NXYZS, self % cwork, self % cwork, plan % dir, fftw_measure)
   end if
#else
   plan % planptr = fftw_plan_gen(NXYZS, self % cwork, self % cwork, plan % dir, fftw_measure)
#endif

else
   if (size(plantypes) < dimensions) then
      write(*, *) "Error: not enough flags when creating POISFFT_PLANXD, there must be one per dimension."
      STOP
   endif

#if defined(MPI) && dimensions > 1
   distr = merge(distributed, .true., present(distributed))

   if (distr) then
      plan % planptr = PFFT_REAL(dimensions, int([NXYZS], c_intptr_t), self % rwork, self % rwork, self % mpi % comm, &
         plantypes, pfft_transposed_none + pfft_measure + pfft_preserve_input)
   else
      plan % planptr = fftw_plan_gen(NXYZS, self % rwork, self % rwork, REALPLANTYPES, fftw_measure)
   end if
#else
   plan % planptr = fftw_plan_gen(NXYZS, self % rwork, self % rwork, REALPLANTYPES, fftw_measure)
#endif

endif

plan % planowner = .true.
plan % distributed = distr

if (.not. c_associated(plan % planptr)) stop "Error, FFT plan not created!"

#undef COLONS
#undef REALPLANTYPES
#undef NXYZS
#undef POISFFT_SOLVERXD
#undef POISFFT_PLANXD
#undef SLICEXD
#undef PFFT_CMPLX
#undef PFFT_REAL
#undef FFTW_CMPLX_MPI_2D
