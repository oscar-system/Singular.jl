function test_call_interpreter()
    print("call_interpreter...")
    Singular.libSingular.call_interpreter("ring r;")
    result = Singular.libSingular.call_interpreter("r;")
    @test result[1] == false
    Singular.libSingular.call_interpreter("poly f = x;")
    result = Singular.libSingular.call_interpreter("f;")
    @test result[1] == false
    @test result[2] == "x\n"
    @test result[3] == ""
    @test result[4] == ""
    @test length(result) == 4
    println("PASS")
    println("")
end
