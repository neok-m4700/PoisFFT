module fft_sp
#define PREC 1

#include "fft-inc.f90"

#undef PREC
end module

module fft_dp
#define PREC 2

#include "fft-inc.f90"

#undef PREC
end module
