module poisfft_constants
  use iso_c_binding
  use iso_fortran_env, only : &
     fin => input_unit, &
     fout => output_unit, &
     ferr => error_unit
  implicit none

  integer, parameter :: dcp =  c_double_complex
  integer, parameter :: drp =  c_double
  integer, parameter :: scp =  c_float_complex
  integer, parameter :: srp =  c_float

  integer, parameter :: poisfft_periodic = 0
  integer, parameter :: poisfft_dirichlet = 1
  integer, parameter :: poisfft_neumann = 2
  integer, parameter :: poisfft_dirichletstag = 3
  integer, parameter :: poisfft_neumannstag = 4

  integer, parameter :: poisfft_spectral = 0
  integer, parameter :: poisfft_finitedifference2 = 2
  integer, parameter :: poisfft_finitedifference4 = 4
end module
