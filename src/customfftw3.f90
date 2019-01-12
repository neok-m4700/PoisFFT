module fftw3
   use iso_c_binding
   integer, parameter :: c_fftw_r2r_kind = c_int32_t

   integer(c_int), parameter :: fftw_r2hc = 0
   integer(c_int), parameter :: fftw_hc2r = 1
   integer(c_int), parameter :: fftw_dht = 2
   integer(c_int), parameter :: fftw_redft00 = 3
   integer(c_int), parameter :: fftw_redft01 = 4
   integer(c_int), parameter :: fftw_redft10 = 5
   integer(c_int), parameter :: fftw_redft11 = 6
   integer(c_int), parameter :: fftw_rodft00 = 7
   integer(c_int), parameter :: fftw_rodft01 = 8
   integer(c_int), parameter :: fftw_rodft10 = 9
   integer(c_int), parameter :: fftw_rodft11 = 10
   integer(c_int), parameter :: fftw_forward = -1
   integer(c_int), parameter :: fftw_backward = +1
   integer(c_int), parameter :: fftw_measure = 0
   integer(c_int), parameter :: fftw_destroy_input = 1
   integer(c_int), parameter :: fftw_unaligned = 2
   integer(c_int), parameter :: fftw_conserve_memory = 4
   integer(c_int), parameter :: fftw_exhaustive = 8
   integer(c_int), parameter :: fftw_preserve_input = 16
   integer(c_int), parameter :: fftw_patient = 32
   integer(c_int), parameter :: fftw_estimate = 64
   integer(c_int), parameter :: fftw_estimate_patient = 128
   integer(c_int), parameter :: fftw_believe_pcost = 256
   integer(c_int), parameter :: fftw_no_dft_r2hc = 512
   integer(c_int), parameter :: fftw_no_nonthreaded = 1024
   integer(c_int), parameter :: fftw_no_buffering = 2048
   integer(c_int), parameter :: fftw_no_indirect_op = 4096
   integer(c_int), parameter :: fftw_allow_large_generic = 8192
   integer(c_int), parameter :: fftw_no_rank_splits = 16384
   integer(c_int), parameter :: fftw_no_vrank_splits = 32768
   integer(c_int), parameter :: fftw_no_vrecurse = 65536
   integer(c_int), parameter :: fftw_no_simd = 131072
   integer(c_int), parameter :: fftw_no_slow = 262144
   integer(c_int), parameter :: fftw_no_fixed_radix_large_n = 524288
   integer(c_int), parameter :: fftw_allow_pruning = 1048576
   integer(c_int), parameter :: fftw_wisdom_only = 2097152

   interface fftw_execute_gen
      subroutine fftw_execute_r2r(p, in, out) bind(c, name='fftw_execute_r2r')
         import
         type(c_ptr), value :: p
         real(c_double), dimension(*), intent(inout) :: in
         real(c_double), dimension(*), intent(inout) :: out
      end subroutine

      subroutine fftwf_execute_r2r(p, in, out) bind(c, name='fftwf_execute_r2r')
         import
         type(c_ptr), value :: p
         real(c_float), dimension(*), intent(inout) :: in
         real(c_float), dimension(*), intent(inout) :: out
      end subroutine

      subroutine fftw_execute_dft(p, in, out) bind(c, name='fftw_execute_dft')
         import
         type(c_ptr), value :: p
         complex(c_double_complex), dimension(*), intent(inout) :: in
         complex(c_double_complex), dimension(*), intent(inout) :: out
      end subroutine

      subroutine fftwf_execute_dft(p, in, out) bind(c, name='fftwf_execute_dft')
         import
         type(c_ptr), value :: p
         complex(c_float_complex), dimension(*), intent(inout) :: in
         complex(c_float_complex), dimension(*), intent(inout) :: out
      end subroutine

      subroutine fftw_execute_dft_r2c(p, in, out) bind(c, name='fftw_execute_dft_r2c')
         import
         type(c_ptr), value :: p
         real(c_double), dimension(*), intent(inout) :: in
         complex(c_double_complex), dimension(*), intent(out) :: out
      end subroutine

      subroutine fftw_execute_dft_c2r(p, in, out) bind(c, name='fftw_execute_dft_c2r')
         import
         type(c_ptr), value :: p
         complex(c_double_complex), dimension(*), intent(inout) :: in
         real(c_double), dimension(*), intent(out) :: out
      end subroutine

      subroutine fftwf_execute_dft_r2c(p, in, out) bind(c, name='fftwf_execute_dft_r2c')
         import
         type(c_ptr), value :: p
         real(c_float), dimension(*), intent(inout) :: in
         complex(c_float_complex), dimension(*), intent(out) :: out
      end subroutine

      subroutine fftwf_execute_dft_c2r(p, in, out) bind(c, name='fftwf_execute_dft_c2r')
         import
         type(c_ptr), value :: p
         complex(c_float_complex), dimension(*), intent(inout) :: in
         real(c_float), dimension(*), intent(out) :: out
      end subroutine
   end interface

   interface fftw_plan_gen
      type(c_ptr) function fftw_plan_dft_1d(n, in, out, sign, flags) bind(c, name='fftw_plan_dft_1d')
         import
         integer(c_int), value :: n
         complex(c_double_complex), dimension(*), intent(inout) :: in
         complex(c_double_complex), dimension(*), intent(inout) :: out
         integer(c_int), value :: sign
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_plan_dft_2d(n0, n1, in, out, sign, flags) bind(c, name='fftw_plan_dft_2d')
         import
         integer(c_int), value :: n0
         integer(c_int), value :: n1
         complex(c_double_complex), dimension(n1, *), intent(inout) :: in
         complex(c_double_complex), dimension(n1, *), intent(inout) :: out
         integer(c_int), value :: sign
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_plan_dft_3d(n0, n1, n2, in, out, sign, flags) bind(c, name='fftw_plan_dft_3d')
         import
         integer(c_int), value :: n0
         integer(c_int), value :: n1
         integer(c_int), value :: n2
         complex(c_double_complex), dimension(n2, n1, *), intent(inout) :: in
         complex(c_double_complex), dimension(n2, n1, *), intent(inout) :: out
         integer(c_int), value :: sign
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_plan_dft_1d(n, in, out, sign, flags) bind(c, name='fftwf_plan_dft_1d')
         import
         integer(c_int), value :: n
         complex(c_float_complex), dimension(*), intent(inout) :: in
         complex(c_float_complex), dimension(*), intent(inout) :: out
         integer(c_int), value :: sign
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_plan_dft_2d(n0, n1, in, out, sign, flags) bind(c, name='fftwf_plan_dft_2d')
         import
         integer(c_int), value :: n0
         integer(c_int), value :: n1
         complex(c_float_complex), dimension(n1, *), intent(inout) :: in
         complex(c_float_complex), dimension(n1, *), intent(inout) :: out
         integer(c_int), value :: sign
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_plan_dft_3d(n0, n1, n2, in, out, sign, flags) bind(c, name='fftwf_plan_dft_3d')
         import
         integer(c_int), value :: n0
         integer(c_int), value :: n1
         integer(c_int), value :: n2
         complex(c_float_complex), dimension(n2, n1, *), intent(inout) :: in
         complex(c_float_complex), dimension(n2, n1, *), intent(inout) :: out
         integer(c_int), value :: sign
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_plan_r2r_1d(n, in, out, kind, flags) bind(c, name='fftw_plan_r2r_1d')
         import
         integer(c_int), value :: n
         real(c_double), dimension(*), intent(inout) :: in
         real(c_double), dimension(*), intent(inout) :: out
         integer(c_fftw_r2r_kind), value :: kind
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_plan_r2r_2d(n0, n1, in, out, kind0, kind1, flags) bind(c, name='fftw_plan_r2r_2d')
         import
         integer(c_int), value :: n0
         integer(c_int), value :: n1
         real(c_double), dimension(n1, *), intent(inout) :: in
         real(c_double), dimension(n1, *), intent(inout) :: out
         integer(c_fftw_r2r_kind), value :: kind0
         integer(c_fftw_r2r_kind), value :: kind1
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_plan_r2r_3d(n0, n1, n2, in, out, kind0, kind1, kind2, flags) bind(c, name='fftw_plan_r2r_3d')
         import
         integer(c_int), value :: n0
         integer(c_int), value :: n1
         integer(c_int), value :: n2
         real(c_double), dimension(n2, n1, *), intent(inout) :: in
         real(c_double), dimension(n2, n1, *), intent(inout) :: out
         integer(c_fftw_r2r_kind), value :: kind0
         integer(c_fftw_r2r_kind), value :: kind1
         integer(c_fftw_r2r_kind), value :: kind2
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_plan_r2r_1d(n, in, out, kind, flags) bind(c, name='fftwf_plan_r2r_1d')
         import
         integer(c_int), value :: n
         real(c_float), dimension(*), intent(inout) :: in
         real(c_float), dimension(*), intent(inout) :: out
         integer(c_fftw_r2r_kind), value :: kind
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_plan_r2r_2d(n0, n1, in, out, kind0, kind1, flags) bind(c, name='fftwf_plan_r2r_2d')
         import
         integer(c_int), value :: n0
         integer(c_int), value :: n1
         real(c_float), dimension(n1, *), intent(inout) :: in
         real(c_float), dimension(n1, *), intent(inout) :: out
         integer(c_fftw_r2r_kind), value :: kind0
         integer(c_fftw_r2r_kind), value :: kind1
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_plan_r2r_3d(n0, n1, n2, in, out, kind0, kind1, kind2, flags) bind(c, name='fftwf_plan_r2r_3d')
         import
         integer(c_int), value :: n0
         integer(c_int), value :: n1
         integer(c_int), value :: n2
         real(c_float), dimension(n2, n1, *), intent(inout) :: in
         real(c_float), dimension(n2, n1, *), intent(inout) :: out
         integer(c_fftw_r2r_kind), value :: kind0
         integer(c_fftw_r2r_kind), value :: kind1
         integer(c_fftw_r2r_kind), value :: kind2
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_plan_many_dft_r2c(rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, flags) bind(c, name='fftw_plan_many_dft_r2c')
         import
         integer(c_int), value :: rank
         integer(c_int), dimension(*), intent(in) :: n
         integer(c_int), value :: howmany
         real(c_double), dimension(*), intent(out) :: in
         integer(c_int), dimension(*), intent(in) :: inembed
         integer(c_int), value :: istride
         integer(c_int), value :: idist
         complex(c_double_complex), dimension(*), intent(out) :: out
         integer(c_int), dimension(*), intent(in) :: onembed
         integer(c_int), value :: ostride
         integer(c_int), value :: odist
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_plan_dft_r2c(rank, n, in, out, flags) bind(c, name='fftw_plan_dft_r2c')
         import
         integer(c_int), value :: rank
         integer(c_int), dimension(*), intent(in) :: n
         real(c_double), dimension(*), intent(out) :: in
         complex(c_double_complex), dimension(*), intent(out) :: out
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_plan_dft_r2c_1d(n, in, out, flags) bind(c, name='fftw_plan_dft_r2c_1d')
         import
         integer(c_int), value :: n
         real(c_double), dimension(*), intent(out) :: in
         complex(c_double_complex), dimension(*), intent(out) :: out
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_plan_dft_r2c_2d(n0, n1, in, out, flags) bind(c, name='fftw_plan_dft_r2c_2d')
         import
         integer(c_int), value :: n0
         integer(c_int), value :: n1
         real(c_double), dimension(n1, *), intent(out) :: in
         complex(c_double_complex), dimension(n1, *), intent(out) :: out
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_plan_dft_r2c_3d(n0, n1, n2, in, out, flags) bind(c, name='fftw_plan_dft_r2c_3d')
         import
         integer(c_int), value :: n0
         integer(c_int), value :: n1
         integer(c_int), value :: n2
         real(c_double), dimension(n2, n1, *), intent(out) :: in
         complex(c_double_complex), dimension(n2, n1, *), intent(out) :: out
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_plan_many_dft_c2r(rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, flags) bind(c, name='fftw_plan_many_dft_c2r')
         import
         integer(c_int), value :: rank
         integer(c_int), dimension(*), intent(in) :: n
         integer(c_int), value :: howmany
         complex(c_double_complex), dimension(*), intent(out) :: in
         integer(c_int), dimension(*), intent(in) :: inembed
         integer(c_int), value :: istride
         integer(c_int), value :: idist
         real(c_double), dimension(*), intent(out) :: out
         integer(c_int), dimension(*), intent(in) :: onembed
         integer(c_int), value :: ostride
         integer(c_int), value :: odist
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_plan_dft_c2r(rank, n, in, out, flags) bind(c, name='fftw_plan_dft_c2r')
         import
         integer(c_int), value :: rank
         integer(c_int), dimension(*), intent(in) :: n
         complex(c_double_complex), dimension(*), intent(out) :: in
         real(c_double), dimension(*), intent(out) :: out
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_plan_dft_c2r_1d(n, in, out, flags) bind(c, name='fftw_plan_dft_c2r_1d')
         import
         integer(c_int), value :: n
         complex(c_double_complex), dimension(*), intent(out) :: in
         real(c_double), dimension(*), intent(out) :: out
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_plan_dft_c2r_2d(n0, n1, in, out, flags) bind(c, name='fftw_plan_dft_c2r_2d')
         import
         integer(c_int), value :: n0
         integer(c_int), value :: n1
         complex(c_double_complex), dimension(n1, *), intent(out) :: in
         real(c_double), dimension(n1, *), intent(out) :: out
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_plan_dft_c2r_3d(n0, n1, n2, in, out, flags) bind(c, name='fftw_plan_dft_c2r_3d')
         import
         integer(c_int), value :: n0
         integer(c_int), value :: n1
         integer(c_int), value :: n2
         complex(c_double_complex), dimension(n2, n1, *), intent(out) :: in
         real(c_double), dimension(n2, n1, *), intent(out) :: out
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_plan_many_dft_r2c(rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, flags) bind(c, name='fftwf_plan_many_dft_r2c')
         import
         integer(c_int), value :: rank
         integer(c_int), dimension(*), intent(in) :: n
         integer(c_int), value :: howmany
         real(c_float), dimension(*), intent(out) :: in
         integer(c_int), dimension(*), intent(in) :: inembed
         integer(c_int), value :: istride
         integer(c_int), value :: idist
         complex(c_float_complex), dimension(*), intent(out) :: out
         integer(c_int), dimension(*), intent(in) :: onembed
         integer(c_int), value :: ostride
         integer(c_int), value :: odist
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_plan_dft_r2c(rank, n, in, out, flags) bind(c, name='fftwf_plan_dft_r2c')
         import
         integer(c_int), value :: rank
         integer(c_int), dimension(*), intent(in) :: n
         real(c_float), dimension(*), intent(out) :: in
         complex(c_float_complex), dimension(*), intent(out) :: out
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_plan_dft_r2c_1d(n, in, out, flags) bind(c, name='fftwf_plan_dft_r2c_1d')
         import
         integer(c_int), value :: n
         real(c_float), dimension(*), intent(out) :: in
         complex(c_float_complex), dimension(*), intent(out) :: out
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_plan_dft_r2c_2d(n0, n1, in, out, flags) bind(c, name='fftwf_plan_dft_r2c_2d')
         import
         integer(c_int), value :: n0
         integer(c_int), value :: n1
         real(c_float), dimension(n1, *), intent(out) :: in
         complex(c_float_complex), dimension(n1, *), intent(out) :: out
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_plan_dft_r2c_3d(n0, n1, n2, in, out, flags) bind(c, name='fftwf_plan_dft_r2c_3d')
         import
         integer(c_int), value :: n0
         integer(c_int), value :: n1
         integer(c_int), value :: n2
         real(c_float), dimension(n2, n1, *), intent(out) :: in
         complex(c_float_complex), dimension(n2, n1, *), intent(out) :: out
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_plan_many_dft_c2r(rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, flags) bind(c, name='fftwf_plan_many_dft_c2r')
         import
         integer(c_int), value :: rank
         integer(c_int), dimension(*), intent(in) :: n
         integer(c_int), value :: howmany
         complex(c_float_complex), dimension(*), intent(out) :: in
         integer(c_int), dimension(*), intent(in) :: inembed
         integer(c_int), value :: istride
         integer(c_int), value :: idist
         real(c_float), dimension(*), intent(out) :: out
         integer(c_int), dimension(*), intent(in) :: onembed
         integer(c_int), value :: ostride
         integer(c_int), value :: odist
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_plan_dft_c2r(rank, n, in, out, flags) bind(c, name='fftwf_plan_dft_c2r')
         import
         integer(c_int), value :: rank
         integer(c_int), dimension(*), intent(in) :: n
         complex(c_float_complex), dimension(*), intent(out) :: in
         real(c_float), dimension(*), intent(out) :: out
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_plan_dft_c2r_1d(n, in, out, flags) bind(c, name='fftwf_plan_dft_c2r_1d')
         import
         integer(c_int), value :: n
         complex(c_float_complex), dimension(*), intent(out) :: in
         real(c_float), dimension(*), intent(out) :: out
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_plan_dft_c2r_2d(n0, n1, in, out, flags) bind(c, name='fftwf_plan_dft_c2r_2d')
         import
         integer(c_int), value :: n0
         integer(c_int), value :: n1
         complex(c_float_complex), dimension(n1, *), intent(out) :: in
         real(c_float), dimension(n1, *), intent(out) :: out
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_plan_dft_c2r_3d(n0, n1, n2, in, out, flags) bind(c, name='fftwf_plan_dft_c2r_3d')
         import
         integer(c_int), value :: n0
         integer(c_int), value :: n1
         integer(c_int), value :: n2
         complex(c_float_complex), dimension(n2, n1, *), intent(out) :: in
         real(c_float), dimension(n2, n1, *), intent(out) :: out
         integer(c_int), value :: flags
      end function

   end interface fftw_plan_gen

   interface
      type(c_ptr) function fftw_plan_many_dft(rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, sign, flags) bind(c, name='fftw_plan_many_dft')
         import
         integer(c_int), value :: rank
         integer(c_int), dimension(*), intent(in) :: n
         integer(c_int), value :: howmany
         complex(c_double_complex), dimension(*), intent(inout) :: in
         integer(c_int), dimension(*), intent(in) :: inembed
         integer(c_int), value :: istride
         integer(c_int), value :: idist
         complex(c_double_complex), dimension(*), intent(inout) :: out
         integer(c_int), dimension(*), intent(in) :: onembed
         integer(c_int), value :: ostride
         integer(c_int), value :: odist
         integer(c_int), value :: sign
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_plan_many_r2r(rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, kind, flags) bind(c, name='fftw_plan_many_r2r')
         import
         integer(c_int), value :: rank
         integer(c_int), dimension(*), intent(in) :: n
         integer(c_int), value :: howmany
         real(c_double), dimension(*), intent(inout) :: in
         integer(c_int), dimension(*), intent(in) :: inembed
         integer(c_int), value :: istride
         integer(c_int), value :: idist
         real(c_double), dimension(*), intent(inout) :: out
         integer(c_int), dimension(*), intent(in) :: onembed
         integer(c_int), value :: ostride
         integer(c_int), value :: odist
         integer(c_fftw_r2r_kind), dimension(*), intent(in) :: kind
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_plan_many_dft(rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, sign, flags) bind(c, name='fftwf_plan_many_dft')
         import
         integer(c_int), value :: rank
         integer(c_int), dimension(*), intent(in) :: n
         integer(c_int), value :: howmany
         complex(c_float_complex), dimension(*), intent(inout) :: in
         integer(c_int), dimension(*), intent(in) :: inembed
         integer(c_int), value :: istride
         integer(c_int), value :: idist
         complex(c_float_complex), dimension(*), intent(inout) :: out
         integer(c_int), dimension(*), intent(in) :: onembed
         integer(c_int), value :: ostride
         integer(c_int), value :: odist
         integer(c_int), value :: sign
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_plan_many_r2r(rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, kind, flags) bind(c, name='fftwf_plan_many_r2r')
         import
         integer(c_int), value :: rank
         integer(c_int), dimension(*), intent(in) :: n
         integer(c_int), value :: howmany
         real(c_float), dimension(*), intent(inout) :: in
         integer(c_int), dimension(*), intent(in) :: inembed
         integer(c_int), value :: istride
         integer(c_int), value :: idist
         real(c_float), dimension(*), intent(inout) :: out
         integer(c_int), dimension(*), intent(in) :: onembed
         integer(c_int), value :: ostride
         integer(c_int), value :: odist
         integer(c_fftw_r2r_kind), dimension(*), intent(in) :: kind
         integer(c_int), value :: flags
      end function
   end interface

   interface
      type(c_ptr) function fftw_malloc(n) bind(c, name='fftw_malloc')
         import
         integer(c_size_t), value :: n
      end function

      subroutine fftw_free(p) bind(c, name='fftw_free')
         import
         type(c_ptr), value :: p
      end subroutine

      subroutine fftwf_destroy_plan(p) bind(c, name='fftwf_destroy_plan')
         import
         type(c_ptr), value :: p
      end subroutine

      subroutine fftw_destroy_plan(p) bind(c, name='fftw_destroy_plan')
         import
         type(c_ptr), value :: p
      end subroutine

      subroutine fftw_plan_with_nthreads(nthreads) bind(c, name='fftw_plan_with_nthreads')
         import
         integer(c_int), value :: nthreads
      end subroutine

      integer(c_int) function fftw_init_threads() bind(c, name='fftw_init_threads')
         import
      end function

      subroutine fftw_cleanup_threads() bind(c, name='fftw_cleanup_threads')
         import
      end subroutine
   end interface

#ifdef MPI
   integer(c_intptr_t), parameter :: fftw_mpi_default_block = 0
   integer(c_int), parameter :: fftw_mpi_scrambled_in = 134217728
   integer(c_int), parameter :: fftw_mpi_scrambled_out = 268435456
   integer(c_int), parameter :: fftw_mpi_transposed_in = 536870912
   integer(c_int), parameter :: fftw_mpi_transposed_out = 1073741824

   type, bind(c) :: fftw_mpi_ddim
      integer(c_intptr_t) n, ib, ob
   end type fftw_mpi_ddim

   interface
      subroutine fftw_mpi_init() bind(c, name='fftw_mpi_init')
         import
      end subroutine

      subroutine fftw_mpi_cleanup() bind(c, name='fftw_mpi_cleanup')
         import
      end subroutine

      integer(c_intptr_t) function fftw_mpi_local_size_many_transposed(rnk, n, howmany, block0, block1, comm, local_n0, local_0_start, local_n1, local_1_start) bind(c, name='fftw_mpi_local_size_many_transposed_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         integer(c_intptr_t), value :: howmany
         integer(c_intptr_t), value :: block0
         integer(c_intptr_t), value :: block1
         integer(c_int32_t), value :: comm
         integer(c_intptr_t), intent(out) :: local_n0
         integer(c_intptr_t), intent(out) :: local_0_start
         integer(c_intptr_t), intent(out) :: local_n1
         integer(c_intptr_t), intent(out) :: local_1_start
      end function

      integer(c_intptr_t) function fftw_mpi_local_size_many(rnk, n, howmany, block0, comm, local_n0, local_0_start) bind(c, name='fftw_mpi_local_size_many_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         integer(c_intptr_t), value :: howmany
         integer(c_intptr_t), value :: block0
         integer(c_int32_t), value :: comm
         integer(c_intptr_t), intent(out) :: local_n0
         integer(c_intptr_t), intent(out) :: local_0_start
      end function

      integer(c_intptr_t) function fftw_mpi_local_size_transposed(rnk, n, comm, local_n0, local_0_start, local_n1, local_1_start) bind(c, name='fftw_mpi_local_size_transposed_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         integer(c_int32_t), value :: comm
         integer(c_intptr_t), intent(out) :: local_n0
         integer(c_intptr_t), intent(out) :: local_0_start
         integer(c_intptr_t), intent(out) :: local_n1
         integer(c_intptr_t), intent(out) :: local_1_start
      end function

      integer(c_intptr_t) function fftw_mpi_local_size(rnk, n, comm, local_n0, local_0_start) bind(c, name='fftw_mpi_local_size_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         integer(c_int32_t), value :: comm
         integer(c_intptr_t), intent(out) :: local_n0
         integer(c_intptr_t), intent(out) :: local_0_start
      end function

      integer(c_intptr_t) function fftw_mpi_local_size_many_1d(n0, howmany, comm, sign, flags, local_ni, local_i_start, local_no, local_o_start) bind(c, name='fftw_mpi_local_size_many_1d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: howmany
         integer(c_int32_t), value :: comm
         integer(c_int), value :: sign
         integer(c_int), value :: flags
         integer(c_intptr_t), intent(out) :: local_ni
         integer(c_intptr_t), intent(out) :: local_i_start
         integer(c_intptr_t), intent(out) :: local_no
         integer(c_intptr_t), intent(out) :: local_o_start
      end function

      integer(c_intptr_t) function fftw_mpi_local_size_1d(n0, comm, sign, flags, local_ni, local_i_start, local_no, local_o_start) bind(c, name='fftw_mpi_local_size_1d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_int32_t), value :: comm
         integer(c_int), value :: sign
         integer(c_int), value :: flags
         integer(c_intptr_t), intent(out) :: local_ni
         integer(c_intptr_t), intent(out) :: local_i_start
         integer(c_intptr_t), intent(out) :: local_no
         integer(c_intptr_t), intent(out) :: local_o_start
      end function

      integer(c_intptr_t) function fftw_mpi_local_size_2d(n0, n1, comm, local_n0, local_0_start) bind(c, name='fftw_mpi_local_size_2d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         integer(c_int32_t), value :: comm
         integer(c_intptr_t), intent(out) :: local_n0
         integer(c_intptr_t), intent(out) :: local_0_start
      end function

      integer(c_intptr_t) function fftw_mpi_local_size_2d_transposed(n0, n1, comm, local_n0, local_0_start, local_n1, local_1_start) bind(c, name='fftw_mpi_local_size_2d_transposed_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         integer(c_int32_t), value :: comm
         integer(c_intptr_t), intent(out) :: local_n0
         integer(c_intptr_t), intent(out) :: local_0_start
         integer(c_intptr_t), intent(out) :: local_n1
         integer(c_intptr_t), intent(out) :: local_1_start
      end function

      integer(c_intptr_t) function fftw_mpi_local_size_3d(n0, n1, n2, comm, local_n0, local_0_start) bind(c, name='fftw_mpi_local_size_3d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         integer(c_intptr_t), value :: n2
         integer(c_int32_t), value :: comm
         integer(c_intptr_t), intent(out) :: local_n0
         integer(c_intptr_t), intent(out) :: local_0_start
      end function

      integer(c_intptr_t) function fftw_mpi_local_size_3d_transposed(n0, n1, n2, comm, local_n0, local_0_start, local_n1, local_1_start) bind(c, name='fftw_mpi_local_size_3d_transposed_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         integer(c_intptr_t), value :: n2
         integer(c_int32_t), value :: comm
         integer(c_intptr_t), intent(out) :: local_n0
         integer(c_intptr_t), intent(out) :: local_0_start
         integer(c_intptr_t), intent(out) :: local_n1
         integer(c_intptr_t), intent(out) :: local_1_start
      end function

      type(c_ptr) function fftw_mpi_plan_many_transpose(n0, n1, howmany, block0, block1, in, out, comm, flags) bind(c, name='fftw_mpi_plan_many_transpose_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         integer(c_intptr_t), value :: howmany
         integer(c_intptr_t), value :: block0
         integer(c_intptr_t), value :: block1
         real(c_double), dimension(*), intent(inout) :: in
         real(c_double), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_mpi_plan_transpose(n0, n1, in, out, comm, flags) bind(c, name='fftw_mpi_plan_transpose_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         real(c_double), dimension(*), intent(inout) :: in
         real(c_double), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function

      subroutine fftw_mpi_gather_wisdom(comm_) bind(c, name='fftw_mpi_gather_wisdom_f03')
         import
         integer(c_int32_t), value :: comm_
      end subroutine

      subroutine fftw_mpi_broadcast_wisdom(comm_) bind(c, name='fftw_mpi_broadcast_wisdom_f03')
         import
         integer(c_int32_t), value :: comm_
      end subroutine

      subroutine fftw_mpi_execute_dft(p, in, out) bind(c, name='fftw_mpi_execute_dft')
         import
         type(c_ptr), value :: p
         complex(c_double_complex), dimension(*), intent(inout) :: in
         complex(c_double_complex), dimension(*), intent(inout) :: out
      end subroutine

      subroutine fftw_mpi_execute_dft_r2c(p, in, out) bind(c, name='fftw_mpi_execute_dft_r2c')
         import
         type(c_ptr), value :: p
         real(c_double), dimension(*), intent(inout) :: in
         complex(c_double_complex), dimension(*), intent(inout) :: out
      end subroutine

      subroutine fftw_mpi_execute_dft_c2r(p, in, out) bind(c, name='fftw_mpi_execute_dft_c2r')
         import
         type(c_ptr), value :: p
         complex(c_double_complex), dimension(*), intent(inout) :: in
         real(c_double), dimension(*), intent(inout) :: out
      end subroutine

      subroutine fftw_mpi_execute_r2r(p, in, out) bind(c, name='fftw_mpi_execute_r2r')
         import
         type(c_ptr), value :: p
         real(c_double), dimension(*), intent(inout) :: in
         real(c_double), dimension(*), intent(inout) :: out
      end subroutine

   end interface

   type, bind(c) :: fftwf_mpi_ddim
      integer(c_intptr_t) n, ib, ob
   end type fftwf_mpi_ddim

   interface
      subroutine fftwf_mpi_init() bind(c, name='fftwf_mpi_init')
         import
      end subroutine

      subroutine fftwf_mpi_cleanup() bind(c, name='fftwf_mpi_cleanup')
         import
      end subroutine

      integer(c_intptr_t) function fftwf_mpi_local_size_many_transposed(rnk, n, howmany, block0, block1, comm, local_n0, local_0_start, local_n1, local_1_start) bind(c, name='fftwf_mpi_local_size_many_transposed_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         integer(c_intptr_t), value :: howmany
         integer(c_intptr_t), value :: block0
         integer(c_intptr_t), value :: block1
         integer(c_int32_t), value :: comm
         integer(c_intptr_t), intent(out) :: local_n0
         integer(c_intptr_t), intent(out) :: local_0_start
         integer(c_intptr_t), intent(out) :: local_n1
         integer(c_intptr_t), intent(out) :: local_1_start
      end function

      integer(c_intptr_t) function fftwf_mpi_local_size_many(rnk, n, howmany, block0, comm, local_n0, local_0_start) bind(c, name='fftwf_mpi_local_size_many_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         integer(c_intptr_t), value :: howmany
         integer(c_intptr_t), value :: block0
         integer(c_int32_t), value :: comm
         integer(c_intptr_t), intent(out) :: local_n0
         integer(c_intptr_t), intent(out) :: local_0_start
      end function

      integer(c_intptr_t) function fftwf_mpi_local_size_transposed(rnk, n, comm, local_n0, local_0_start, local_n1, local_1_start) bind(c, name='fftwf_mpi_local_size_transposed_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         integer(c_int32_t), value :: comm
         integer(c_intptr_t), intent(out) :: local_n0
         integer(c_intptr_t), intent(out) :: local_0_start
         integer(c_intptr_t), intent(out) :: local_n1
         integer(c_intptr_t), intent(out) :: local_1_start
      end function

      integer(c_intptr_t) function fftwf_mpi_local_size(rnk, n, comm, local_n0, local_0_start) bind(c, name='fftwf_mpi_local_size_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         integer(c_int32_t), value :: comm
         integer(c_intptr_t), intent(out) :: local_n0
         integer(c_intptr_t), intent(out) :: local_0_start
      end function

      integer(c_intptr_t) function fftwf_mpi_local_size_many_1d(n0, howmany, comm, sign, flags, local_ni, local_i_start, local_no, local_o_start) bind(c, name='fftwf_mpi_local_size_many_1d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: howmany
         integer(c_int32_t), value :: comm
         integer(c_int), value :: sign
         integer(c_int), value :: flags
         integer(c_intptr_t), intent(out) :: local_ni
         integer(c_intptr_t), intent(out) :: local_i_start
         integer(c_intptr_t), intent(out) :: local_no
         integer(c_intptr_t), intent(out) :: local_o_start
      end function

      integer(c_intptr_t) function fftwf_mpi_local_size_1d(n0, comm, sign, flags, local_ni, local_i_start, local_no, local_o_start) bind(c, name='fftwf_mpi_local_size_1d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_int32_t), value :: comm
         integer(c_int), value :: sign
         integer(c_int), value :: flags
         integer(c_intptr_t), intent(out) :: local_ni
         integer(c_intptr_t), intent(out) :: local_i_start
         integer(c_intptr_t), intent(out) :: local_no
         integer(c_intptr_t), intent(out) :: local_o_start
      end function

      integer(c_intptr_t) function fftwf_mpi_local_size_2d(n0, n1, comm, local_n0, local_0_start) bind(c, name='fftwf_mpi_local_size_2d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         integer(c_int32_t), value :: comm
         integer(c_intptr_t), intent(out) :: local_n0
         integer(c_intptr_t), intent(out) :: local_0_start
      end function

      integer(c_intptr_t) function fftwf_mpi_local_size_2d_transposed(n0, n1, comm, local_n0, local_0_start, local_n1, local_1_start) bind(c, name='fftwf_mpi_local_size_2d_transposed_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         integer(c_int32_t), value :: comm
         integer(c_intptr_t), intent(out) :: local_n0
         integer(c_intptr_t), intent(out) :: local_0_start
         integer(c_intptr_t), intent(out) :: local_n1
         integer(c_intptr_t), intent(out) :: local_1_start
      end function

      integer(c_intptr_t) function fftwf_mpi_local_size_3d(n0, n1, n2, comm, local_n0, local_0_start) bind(c, name='fftwf_mpi_local_size_3d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         integer(c_intptr_t), value :: n2
         integer(c_int32_t), value :: comm
         integer(c_intptr_t), intent(out) :: local_n0
         integer(c_intptr_t), intent(out) :: local_0_start
      end function

      integer(c_intptr_t) function fftwf_mpi_local_size_3d_transposed(n0, n1, n2, comm, local_n0, local_0_start, local_n1, local_1_start) bind(c, name='fftwf_mpi_local_size_3d_transposed_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         integer(c_intptr_t), value :: n2
         integer(c_int32_t), value :: comm
         integer(c_intptr_t), intent(out) :: local_n0
         integer(c_intptr_t), intent(out) :: local_0_start
         integer(c_intptr_t), intent(out) :: local_n1
         integer(c_intptr_t), intent(out) :: local_1_start
      end function

      type(c_ptr) function fftwf_mpi_plan_many_transpose(n0, n1, howmany, block0, block1, in, out, comm, flags) bind(c, name='fftwf_mpi_plan_many_transpose_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         integer(c_intptr_t), value :: howmany
         integer(c_intptr_t), value :: block0
         integer(c_intptr_t), value :: block1
         real(c_float), dimension(*), intent(inout) :: in
         real(c_float), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_mpi_plan_transpose(n0, n1, in, out, comm, flags) bind(c, name='fftwf_mpi_plan_transpose_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         real(c_float), dimension(*), intent(inout) :: in
         real(c_float), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function



      subroutine fftwf_mpi_gather_wisdom(comm_) bind(c, name='fftwf_mpi_gather_wisdom_f03')
         import
         integer(c_int32_t), value :: comm_
      end subroutine

      subroutine fftwf_mpi_broadcast_wisdom(comm_) bind(c, name='fftwf_mpi_broadcast_wisdom_f03')
         import
         integer(c_int32_t), value :: comm_
      end subroutine

      subroutine fftwf_mpi_execute_dft(p, in, out) bind(c, name='fftwf_mpi_execute_dft')
         import
         type(c_ptr), value :: p
         complex(c_float_complex), dimension(*), intent(inout) :: in
         complex(c_float_complex), dimension(*), intent(inout) :: out
      end subroutine

      subroutine fftwf_mpi_execute_dft_r2c(p, in, out) bind(c, name='fftwf_mpi_execute_dft_r2c')
         import
         type(c_ptr), value :: p
         real(c_float), dimension(*), intent(inout) :: in
         complex(c_float_complex), dimension(*), intent(inout) :: out
      end subroutine

      subroutine fftwf_mpi_execute_dft_c2r(p, in, out) bind(c, name='fftwf_mpi_execute_dft_c2r')
         import
         type(c_ptr), value :: p
         complex(c_float_complex), dimension(*), intent(inout) :: in
         real(c_float), dimension(*), intent(inout) :: out
      end subroutine

      subroutine fftwf_mpi_execute_r2r(p, in, out) bind(c, name='fftwf_mpi_execute_r2r')
         import
         type(c_ptr), value :: p
         real(c_float), dimension(*), intent(inout) :: in
         real(c_float), dimension(*), intent(inout) :: out
      end subroutine
   end interface

   interface fftw_mpi_plan_gen
      type(c_ptr) function fftw_mpi_plan_many_dft(rnk, n, howmany, block, tblock, in, out, comm, sign, flags) bind(c, name='fftw_mpi_plan_many_dft_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         integer(c_intptr_t), value :: howmany
         integer(c_intptr_t), value :: block
         integer(c_intptr_t), value :: tblock
         complex(c_double_complex), dimension(*), intent(inout) :: in
         complex(c_double_complex), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: sign
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_mpi_plan_dft(rnk, n, in, out, comm, sign, flags) bind(c, name='fftw_mpi_plan_dft_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         complex(c_double_complex), dimension(*), intent(inout) :: in
         complex(c_double_complex), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: sign
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_mpi_plan_dft_1d(n0, in, out, comm, sign, flags) bind(c, name='fftw_mpi_plan_dft_1d_f03')
         import
         integer(c_intptr_t), value :: n0
         complex(c_double_complex), dimension(*), intent(inout) :: in
         complex(c_double_complex), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: sign
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_mpi_plan_dft_2d(n0, n1, in, out, comm, sign, flags) bind(c, name='fftw_mpi_plan_dft_2d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         complex(c_double_complex), dimension(*), intent(inout) :: in
         complex(c_double_complex), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: sign
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_mpi_plan_dft_3d(n0, n1, n2, in, out, comm, sign, flags) bind(c, name='fftw_mpi_plan_dft_3d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         integer(c_intptr_t), value :: n2
         complex(c_double_complex), dimension(*), intent(inout) :: in
         complex(c_double_complex), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: sign
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_mpi_plan_many_r2r(rnk, n, howmany, iblock, oblock, in, out, comm, kind, flags) bind(c, name='fftw_mpi_plan_many_r2r_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         integer(c_intptr_t), value :: howmany
         integer(c_intptr_t), value :: iblock
         integer(c_intptr_t), value :: oblock
         real(c_double), dimension(*), intent(inout) :: in
         real(c_double), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_fftw_r2r_kind), dimension(*), intent(in) :: kind
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_mpi_plan_r2r(rnk, n, in, out, comm, kind, flags) bind(c, name='fftw_mpi_plan_r2r_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         real(c_double), dimension(*), intent(inout) :: in
         real(c_double), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_fftw_r2r_kind), dimension(*), intent(in) :: kind
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_mpi_plan_r2r_2d(n0, n1, in, out, comm, kind0, kind1, flags) bind(c, name='fftw_mpi_plan_r2r_2d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         real(c_double), dimension(*), intent(inout) :: in
         real(c_double), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_fftw_r2r_kind), value :: kind0
         integer(c_fftw_r2r_kind), value :: kind1
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_mpi_plan_r2r_3d(n0, n1, n2, in, out, comm, kind0, kind1, kind2, flags) bind(c, name='fftw_mpi_plan_r2r_3d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         integer(c_intptr_t), value :: n2
         real(c_double), dimension(*), intent(inout) :: in
         real(c_double), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_fftw_r2r_kind), value :: kind0
         integer(c_fftw_r2r_kind), value :: kind1
         integer(c_fftw_r2r_kind), value :: kind2
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_mpi_plan_many_dft_r2c(rnk, n, howmany, iblock, oblock, in, out, comm, flags) bind(c, name='fftw_mpi_plan_many_dft_r2c_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         integer(c_intptr_t), value :: howmany
         integer(c_intptr_t), value :: iblock
         integer(c_intptr_t), value :: oblock
         real(c_double), dimension(*), intent(inout) :: in
         complex(c_double_complex), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_mpi_plan_dft_r2c(rnk, n, in, out, comm, flags) bind(c, name='fftw_mpi_plan_dft_r2c_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         real(c_double), dimension(*), intent(inout) :: in
         complex(c_double_complex), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_mpi_plan_dft_r2c_2d(n0, n1, in, out, comm, flags) bind(c, name='fftw_mpi_plan_dft_r2c_2d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         real(c_double), dimension(*), intent(inout) :: in
         complex(c_double_complex), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_mpi_plan_dft_r2c_3d(n0, n1, n2, in, out, comm, flags) bind(c, name='fftw_mpi_plan_dft_r2c_3d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         integer(c_intptr_t), value :: n2
         real(c_double), dimension(*), intent(inout) :: in
         complex(c_double_complex), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_mpi_plan_many_dft_c2r(rnk, n, howmany, iblock, oblock, in, out, comm, flags) bind(c, name='fftw_mpi_plan_many_dft_c2r_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         integer(c_intptr_t), value :: howmany
         integer(c_intptr_t), value :: iblock
         integer(c_intptr_t), value :: oblock
         complex(c_double_complex), dimension(*), intent(inout) :: in
         real(c_double), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_mpi_plan_dft_c2r(rnk, n, in, out, comm, flags) bind(c, name='fftw_mpi_plan_dft_c2r_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         complex(c_double_complex), dimension(*), intent(inout) :: in
         real(c_double), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_mpi_plan_dft_c2r_2d(n0, n1, in, out, comm, flags) bind(c, name='fftw_mpi_plan_dft_c2r_2d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         complex(c_double_complex), dimension(*), intent(inout) :: in
         real(c_double), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftw_mpi_plan_dft_c2r_3d(n0, n1, n2, in, out, comm, flags) bind(c, name='fftw_mpi_plan_dft_c2r_3d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         integer(c_intptr_t), value :: n2
         complex(c_double_complex), dimension(*), intent(inout) :: in
         real(c_double), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function



      type(c_ptr) function fftwf_mpi_plan_many_dft(rnk, n, howmany, block, tblock, in, out, comm, sign, flags) bind(c, name='fftwf_mpi_plan_many_dft_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         integer(c_intptr_t), value :: howmany
         integer(c_intptr_t), value :: block
         integer(c_intptr_t), value :: tblock
         complex(c_float_complex), dimension(*), intent(inout) :: in
         complex(c_float_complex), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: sign
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_mpi_plan_dft(rnk, n, in, out, comm, sign, flags) bind(c, name='fftwf_mpi_plan_dft_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         complex(c_float_complex), dimension(*), intent(inout) :: in
         complex(c_float_complex), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: sign
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_mpi_plan_dft_1d(n0, in, out, comm, sign, flags) bind(c, name='fftwf_mpi_plan_dft_1d_f03')
         import
         integer(c_intptr_t), value :: n0
         complex(c_float_complex), dimension(*), intent(inout) :: in
         complex(c_float_complex), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: sign
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_mpi_plan_dft_2d(n0, n1, in, out, comm, sign, flags) bind(c, name='fftwf_mpi_plan_dft_2d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         complex(c_float_complex), dimension(*), intent(inout) :: in
         complex(c_float_complex), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: sign
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_mpi_plan_dft_3d(n0, n1, n2, in, out, comm, sign, flags) bind(c, name='fftwf_mpi_plan_dft_3d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         integer(c_intptr_t), value :: n2
         complex(c_float_complex), dimension(*), intent(inout) :: in
         complex(c_float_complex), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: sign
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_mpi_plan_many_r2r(rnk, n, howmany, iblock, oblock, in, out, comm, kind, flags) bind(c, name='fftwf_mpi_plan_many_r2r_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         integer(c_intptr_t), value :: howmany
         integer(c_intptr_t), value :: iblock
         integer(c_intptr_t), value :: oblock
         real(c_float), dimension(*), intent(inout) :: in
         real(c_float), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_fftw_r2r_kind), dimension(*), intent(in) :: kind
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_mpi_plan_r2r(rnk, n, in, out, comm, kind, flags) bind(c, name='fftwf_mpi_plan_r2r_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         real(c_float), dimension(*), intent(inout) :: in
         real(c_float), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_fftw_r2r_kind), dimension(*), intent(in) :: kind
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_mpi_plan_r2r_2d(n0, n1, in, out, comm, kind0, kind1, flags) bind(c, name='fftwf_mpi_plan_r2r_2d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         real(c_float), dimension(*), intent(inout) :: in
         real(c_float), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_fftw_r2r_kind), value :: kind0
         integer(c_fftw_r2r_kind), value :: kind1
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_mpi_plan_r2r_3d(n0, n1, n2, in, out, comm, kind0, kind1, kind2, flags) bind(c, name='fftwf_mpi_plan_r2r_3d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         integer(c_intptr_t), value :: n2
         real(c_float), dimension(*), intent(inout) :: in
         real(c_float), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_fftw_r2r_kind), value :: kind0
         integer(c_fftw_r2r_kind), value :: kind1
         integer(c_fftw_r2r_kind), value :: kind2
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_mpi_plan_many_dft_r2c(rnk, n, howmany, iblock, oblock, in, out, comm, flags) bind(c, name='fftwf_mpi_plan_many_dft_r2c_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         integer(c_intptr_t), value :: howmany
         integer(c_intptr_t), value :: iblock
         integer(c_intptr_t), value :: oblock
         real(c_float), dimension(*), intent(inout) :: in
         complex(c_float_complex), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_mpi_plan_dft_r2c(rnk, n, in, out, comm, flags) bind(c, name='fftwf_mpi_plan_dft_r2c_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         real(c_float), dimension(*), intent(inout) :: in
         complex(c_float_complex), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_mpi_plan_dft_r2c_2d(n0, n1, in, out, comm, flags) bind(c, name='fftwf_mpi_plan_dft_r2c_2d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         real(c_float), dimension(*), intent(inout) :: in
         complex(c_float_complex), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_mpi_plan_dft_r2c_3d(n0, n1, n2, in, out, comm, flags) bind(c, name='fftwf_mpi_plan_dft_r2c_3d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         integer(c_intptr_t), value :: n2
         real(c_float), dimension(*), intent(inout) :: in
         complex(c_float_complex), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_mpi_plan_many_dft_c2r(rnk, n, howmany, iblock, oblock, in, out, comm, flags) bind(c, name='fftwf_mpi_plan_many_dft_c2r_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         integer(c_intptr_t), value :: howmany
         integer(c_intptr_t), value :: iblock
         integer(c_intptr_t), value :: oblock
         complex(c_float_complex), dimension(*), intent(inout) :: in
         real(c_float), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_mpi_plan_dft_c2r(rnk, n, in, out, comm, flags) bind(c, name='fftwf_mpi_plan_dft_c2r_f03')
         import
         integer(c_int), value :: rnk
         integer(c_intptr_t), dimension(*), intent(in) :: n
         complex(c_float_complex), dimension(*), intent(inout) :: in
         real(c_float), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_mpi_plan_dft_c2r_2d(n0, n1, in, out, comm, flags) bind(c, name='fftwf_mpi_plan_dft_c2r_2d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         complex(c_float_complex), dimension(*), intent(inout) :: in
         real(c_float), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function

      type(c_ptr) function fftwf_mpi_plan_dft_c2r_3d(n0, n1, n2, in, out, comm, flags) bind(c, name='fftwf_mpi_plan_dft_c2r_3d_f03')
         import
         integer(c_intptr_t), value :: n0
         integer(c_intptr_t), value :: n1
         integer(c_intptr_t), value :: n2
         complex(c_float_complex), dimension(*), intent(inout) :: in
         real(c_float), dimension(*), intent(inout) :: out
         integer(c_int32_t), value :: comm
         integer(c_int), value :: flags
      end function
   end interface

#endif
end module
