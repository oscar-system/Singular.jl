@testset "freealg" begin
  R, (u11,u12,u13,u14,
  u21,u22,u23,u24,
  u31,u32,u33,u34,
  u41,u42,u43,u44) = FreeAlgebra(QQ, ["u11", "u12", "u13", "u14",
  "u21", "u22", "u23", "u24",
  "u31", "u32", "u33", "u34",
  "u41", "u42", "u43", "u44"], 7)
  #rs1 = u11 + u12 + u13 + u14 - 1
  rs2 = u21 + u22 + u23 + u24 - 1
  rs3 = u31 + u32 + u33 + u34 - 1
  rs4 = u41 + u42 + u43 + u44 - 1
  cs1 = u11 + u21 + u31 + u41 - 1
  cs2 = u12 + u22 + u32 + u42 - 1
  cs3 = u13 + u23 + u33 + u43 - 1
  cs4 = u14 + u24 + u34 + u44 - 1
  
  Singular.libSingular.set_option("OPT_REDTAIL",false,R.ptr)
  
  J1 = Ideal(R,[rs2,rs3,rs4,cs1,cs2,cs3,cs4])
  J1=std(J1)
  
  @test J1[7] == u11 + u21 + u31 + u41 - 1
  
end
