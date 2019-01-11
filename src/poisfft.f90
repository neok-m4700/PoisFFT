
module poisfft_sp
   use FFT_SP
#define PREC 1

#include "poisfft-inc.f90"

#undef PREC
end module

module poisfft_dp
   use FFT_DP
#define PREC 2

#include "poisfft-inc.f90"

#undef PREC
end module
