#if dims==1
#define colons :
#elif dims==2
#define colons :,:
#else
#define colons :,:,:
#endif

type(c_ptr), value :: d
type(c_ptr), value :: phi, rhs
real(rp), pointer :: f_phi(colons), f_rhs(colons)
integer(c_int), optional :: ngphi(dims), ngrhs(dims)
integer :: i

call c_f_pointer(d, f_d)

if (present(ngphi)) then
   call c_f_pointer(phi, f_phi, [(f_d % nxyz(i) + 2 * ngphi(dims + 1 - i), i=1, dims)])
else
   call c_f_pointer(phi, f_phi, f_d % nxyz)
end if

if (present(ngrhs)) then
   call c_f_pointer(rhs, f_rhs, [(f_d % nxyz(i) + 2 * ngrhs(dims + 1 - i), i=1, dims)])
else
   call c_f_pointer(rhs, f_rhs, f_d % nxyz)
end if

call execute(f_d, f_phi, f_rhs)
#undef colons
