module poisfft_sp
   use fft_sp
#define PREC 1
#include "poisfft-inc.f90"
#undef PREC
end module

module poisfft_dp
   use fft_dp
#define PREC 2
#include "poisfft-inc.f90"
#undef PREC
end module
