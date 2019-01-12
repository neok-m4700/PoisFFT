module pfft
#ifdef MPI
   use fftw3

#include "pfft.f03"

!NOTE: only for older versions of PFFT
#ifdef MISSING_PFFT_PLAN_WITH_NTHREADS
   interface
      subroutine pfft_plan_with_nthreads(nthreads) bind(c, name="pfft_plan_with_nthreads")
         import
         integer(c_int), value :: nthreads
      end subroutine
   end interface

#endif

!NOTE: only for older versions of PFFT
#ifdef MISSING_PFFT_R2R
   interface
      type(c_ptr) function pfft_plan_r2r(rnk, nos, in, out, comm_cart, kinds, pfft_flags) bind(c, name='pfft_plan_r2r_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: nos
         real(c_double), dimension(*), intent(out) :: in
         real(c_double), dimension(*), intent(out) :: out
         integer(c_int32_t), value :: comm_cart
         integer(c_fftw_r2r_kind), dimension(*), intent(in) :: kinds
         integer(c_int), value :: pfft_flags
      end function

      type(c_ptr) function pfftf_plan_r2r(rnk, nos, in, out, comm_cart, kinds, pfft_flags) bind(c, name='pfftf_plan_r2r_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: nos
         real(c_float), dimension(*), intent(out) :: in
         real(c_float), dimension(*), intent(out) :: out
         integer(c_int32_t), value :: comm_cart
         integer(c_fftw_r2r_kind), dimension(*), intent(in) :: kinds
         integer(c_int), value :: pfft_flags
      end function
   end interface

#endif

#endif
end module