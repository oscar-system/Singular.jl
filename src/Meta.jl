#=
   This file generates (at precompilation time) wrappers for all Singular library functions
   All wrappers are given by LIBRARYNAME.functionname. We currently use all caps library names
   to avoid confusion between the LIBRARYNAME modules and other types, for example "Ideal" or "Poly".
=#

# This script calls singulars parse_libs and creates a dictionary
# containing Singulars libraries and the functions contained in those libraries.
include("libraryfuncdictionary.jl")

#=
   In the dictionary `input_manipulator_funcs` for each library a dictionary
   can be added, and for each function in this library a preprocesing function
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
for (name, funcs) in libraryfunctiondictionary
    name_caps = Symbol( "Lib" * uppercasefirst(string(name)))
    func_calls = Any[]
    name_string = string(name) * ".lib"
    for i in funcs
        if i[1] == "g"
            func_name = i[2]
            symb = Symbol(func_name)
            input_manipulator = get(get(input_manipulator_funcs, name, Dict()), symb, identity)
            output_manipulator = get(get(output_manipulator_funcs, name, Dict()), symb, identity)
            push!(func_calls, :($symb(args...) = $(output_manipulator)(low_level_caller($(name_string), $func_name, $(input_manipulator)(args)...)) ))
            push!(func_calls, :($symb(ring::PolyRing, args...) = $(output_manipulator)(low_level_caller_rng($(name_string), $func_name, ring, $(input_manipulator)(args)...)) ))
        end
    end
    eval(:(baremodule $name_caps
        import ..Singular: PolyRing, low_level_caller, low_level_caller_rng
        import Base: *
        $(func_calls...)
    end))
end
