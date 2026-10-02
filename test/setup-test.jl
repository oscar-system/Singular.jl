@testset "omalloc page size compatibility" begin
   wrapper_page_size = Singular.libSingular.omalloc_page_size()
   @test wrapper_page_size == Singular.Setup.omalloc_page_size()
   @test isnothing(Singular.libSingular._check_omalloc_page_size())

   mismatched_page_size = wrapper_page_size == 4096 ? 16384 : 4096
   mktemp() do config_path, io
      write(io, "#define SIZEOF_SYSTEM_PAGE $(mismatched_page_size)\n")
      close(io)

      exception = try
         Singular.libSingular._check_omalloc_page_size(config_path)
         nothing
      catch error
         error
      end
      @test exception isa ErrorException
      @test occursin("requires $(wrapper_page_size) bytes", exception.msg)
      @test occursin("built for $(mismatched_page_size) bytes", exception.msg)
   end

   mktemp() do config_path, io
      write(io, "#define SOME_OTHER_SETTING 4096\n")
      close(io)
      @test_throws ErrorException Singular.Setup.omalloc_page_size(config_path)
   end
end
