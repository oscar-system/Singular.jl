@testset "singular_jll_is_overridden" begin
   # The prebuilt libsingular_julia must not be reused against a Singular_jll
   # that has been overridden with a custom Singular build.
   mktempdir() do tmpdir
      depot = joinpath(tmpdir, "depot")
      overrides_toml = joinpath(depot, "artifacts", "Overrides.toml")
      mkpath(dirname(overrides_toml))
      write(
         overrides_toml,
         """
         [$(Base.PkgId(Singular.Singular_jll).uuid)]
         Singular = "$(joinpath(tmpdir, "singular_override"))"
         """,
      )

      artifacts = Singular.Setup.Pkg.Artifacts
      old_depot_path = copy(DEPOT_PATH)
      old_overrides = deepcopy(artifacts.ARTIFACT_OVERRIDES[])

      try
         empty!(DEPOT_PATH)
         push!(DEPOT_PATH, depot)
         artifacts.ARTIFACT_OVERRIDES[] = nothing
         @test Singular.Setup.singular_jll_is_overridden()

         # same depot, but without the override
         rm(overrides_toml)
         artifacts.ARTIFACT_OVERRIDES[] = nothing
         @test !Singular.Setup.singular_jll_is_overridden()
      finally
         empty!(DEPOT_PATH)
         append!(DEPOT_PATH, old_depot_path)
         artifacts.ARTIFACT_OVERRIDES[] = old_overrides
      end
   end
end
