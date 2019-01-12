#if (DIM==1)
#define COLONS :
#elif (DIM==2)
#define COLONS :,:
#else
#define COLONS :,:,:
#endif

type(c_ptr), value :: self
type(c_ptr), value :: phi, rhs
real(RPC), pointer :: f_phi(COLONS), f_rhs(COLONS)
integer(c_int), optional :: ngphi(DIM), ngrhs(DIM)
integer :: i

call c_f_pointer(self, f_self)

if (present(ngphi)) then
   call c_f_pointer(phi, f_phi, [(f_self % nxyz(i) + 2 * ngphi(DIM + 1 - i), i=1, DIM)])
else
   call c_f_pointer(phi, f_phi, f_self % nxyz)
end if

if (present(ngrhs)) then
   call c_f_pointer(rhs, f_rhs, [(f_self % nxyz(i) + 2 * ngrhs(DIM + 1 - i), i=1, DIM)])
else
   call c_f_pointer(rhs, f_rhs, f_self % nxyz)
end if

call execute(f_self, f_phi, f_rhs)
#undef COLONS
