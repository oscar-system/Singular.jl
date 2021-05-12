function check_error()
   if libSingular.have_error()
      error(libSingular.get_and_clear_error())
   end
end
