include("libraryfuncdictionary.jl")

input_manipulator_funcs = Dict(
    :dummy => Dict(
        # :dummy => i->[ x + 1 for x in i]
    )
)

output_manipulator_funcs = Dict(
    :dummy => Dict(
        # :dummy => i -> i + 1
    )
)

for (name,funcs) in libraryfunctiondictionary
    name_caps = Symbol( "Lib" * uppercasefirst(string(name)))
    func_calls = Any[]
    name_string = string(name) * ".lib"
    for i in funcs
        if i[1] == "g"
            func_name = i[2]
            symb = Symbol(func_name)
            input_manipulator = haskey(input_manipulator_funcs,name) && haskey(input_manipulator_funcs[name],symb) ? input_manipulator_funcs[name][symb] : identity
            output_manipulator = haskey(output_manipulator_funcs,name) && haskey(output_manipulator_funcs[name],symb) ? output_manipulator_funcs[name][symb] : identity
            push!(func_calls, :($symb(args...) = $(output_manipulator)(low_level_caller($(name_string),$func_name,$(input_manipulator)(args)...)) ))
            push!(func_calls, :($symb(ring::PolyRing,args...) = $(output_manipulator)(low_level_caller_rng($(name_string),$func_name,ring,$(input_manipulator)(args)...)) ))
        end
    end
    eval(:(baremodule $name_caps
        import ..Singular: PolyRing, low_level_caller, low_level_caller_rng
        import Base: *
        $(func_calls...)
    end))
end
