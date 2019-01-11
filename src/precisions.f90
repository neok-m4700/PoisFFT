module poisfft_precisions
  use iso_c_binding
  implicit none

  integer, parameter :: dcp =  c_double_complex
  integer, parameter :: drp =  c_double

  integer, parameter :: scp =  c_float_complex
  integer, parameter :: srp =  c_float

end module
