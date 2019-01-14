type(c_ptr), intent(inout) :: self

call c_f_pointer(self, f_self)
call finalize(f_self)
deallocate(f_self)
self = c_null_ptr
