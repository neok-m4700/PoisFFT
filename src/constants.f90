module poisfft_constants
   use iso_c_binding
   use iso_fortran_env, only : &
      fin => input_unit, &
      fout => output_unit, &
      ferr => error_unit
   implicit none

   enum, bind(c)
      enumerator :: dcp = c_double_complex, drp = c_double, scp = c_float_complex, srp = c_float
   end enum

   enum, bind(c)
      enumerator :: poisfft_periodic = 0, &
         poisfft_dirichlet = 1, poisfft_neumann = 2, &
         poisfft_dirichletstag = 3, poisfft_neumannstag = 4
   end enum

   enum, bind(c)
      enumerator :: poisfft_spectral = 0, poisfft_finitedifference2 = 2, poisfft_finitedifference4 = 4
   end enum

end module
