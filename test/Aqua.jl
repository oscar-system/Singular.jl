using Aqua

@testset "Aqua.jl" begin
   Aqua.test_all(
      Singular;
      ambiguities=false,         # TODO: fix ambiguities
      unbound_args=false,        # TODO: fix unbound args
   )
end
