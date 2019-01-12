module poisfft_parameters
   integer, parameter :: poisfft_periodic = 0
   integer, parameter :: poisfft_dirichlet = 1
   integer, parameter :: poisfft_neumann = 2
   integer, parameter :: poisfft_dirichletstag = 3
   integer, parameter :: poisfft_neumannstag = 4

   integer, parameter :: poisfft_spectral = 0
   integer, parameter :: poisfft_finitedifference2 = 2
   integer, parameter :: poisfft_finitedifference4 = 4
end module
