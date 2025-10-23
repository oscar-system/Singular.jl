#=
   This file generates (at precompilation time) wrappers for all Singular library functions
   All wrappers are given by LIBRARYNAME.functionname. We currently use all caps library names
   to avoid confusion between the LIBRARYNAME modules and other types, for example "Ideal" or "Poly".
=#

#=
   In the dictionary `input_manipulator_funcs` for each library a dictionary
   can be added, and for each function in this library a preprocessing function
   can be added. The function gets executed on the input, before the library function
   is called.
=#
input_manipulator_funcs = Dict(
    :dummy => Dict(
        # :dummy => i->[ x + 1 for x in i]
    )
)

# The same as `input_manipulator_funcs`, but the function gets executed on the output of
# a function.
output_manipulator_funcs = Dict(
    :dummy => Dict(
        # :dummy => i -> i + 1
    )
)

#=
   For each library `lib` a module `LIB` is created, which contains wrappers for all globally exposed
   functions of the library. Each library function `foo` can be called either as `foo(ring, args...)`,
   where `ring` is the ring the arguments belong to, or as `foo(args...)` if no ring is needed or the
   ring can be determined from the input arguments. Furthermore, if applicable, input and output manipulator
   functions are added to the call.
=#
for (libname, funcs) in Setup.libraryfunctiondictionary
    libname_caps = Symbol( "Lib" * uppercasefirst(string(libname)))
    func_calls = Any[]
    libname_dot_lib = string(libname) * ".lib"
    for i in funcs
        if i[1] == "g"
            func_name = i[2]
            symb = Symbol(func_name)
            input_manipulator = get(get(input_manipulator_funcs, libname, Dict()), symb, identity)
            output_manipulator = get(get(output_manipulator_funcs, libname, Dict()), symb, identity)
            push!(func_calls, :($symb(args...) = $(output_manipulator)(low_level_caller($(libname_dot_lib), $func_name, $(input_manipulator)(args))) ))
            push!(func_calls, :($symb(ring::PolyRingUnion, args...) = $(output_manipulator)(low_level_caller_ring($(libname_dot_lib), $func_name, ring, $(input_manipulator)(args))) ))
        end
    end
    eval(:(baremodule $libname_caps
        import ..Singular: PolyRingUnion, low_level_caller, low_level_caller_ring
        import Base: *
        $(func_calls...)
    end))
end
