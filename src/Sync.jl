# This module provides an automatic wrapping mechanism for Singular
# that supports mutual exclusion in multi-threaded code.
#
# In order to use it, simply substitute calls to
#
#       Singular.f(...)
#
# with calls to:
#
#       Singular.Sync.f(...)
#
# Future versions may offer automatic locking for Singular by default
# in multi-threaded programs.

module Sync

    # import the normal modules in use by Singular so that forwarded
    # definitions have access to the necessary types defined by them.

    using Singular
    using AbstractAlgebra
    using Nemo
    using LinearAlgebra
    using CxxWrap
    using SuiteSparse
    using SharedArrays
    using SparseArrays

    # We need to import reduce() because we will later overload it.
    import Base: reduce

    const _mutex = ReentrantLock()
    const _reserved = Set([
        :eval, :include, :__init__,
    ])

    function _lock()
        Base.lock(_mutex)
    end

    function _unlock()
        Base.unlock(_mutex)
    end

    function _wrap_methods(mod::Module, fn::Function)
        function qualified_typename(typename)
            if typename isa Core.TypeName
              reduce((a, b) -> :($a.$b),
                  [ fullname(typename.module)..., typename.name ])
            else
              typename
            end
        end
        function bounded(typevar::TypeVar)
            ub = strip_type_bounds(typevar.ub)
            name = typevar.name
            lb = typevar.lb
            if typevar.lb === Union{}
                if ub === Any
                    :( $(name) )
                else
                    :( $(name) <: $(ub) )
                end
            else
                lb = strip_type_bounds(lb)
                :( $(lb) <: $(name) <: $(ub) )
            end
        end
        function isunbounded(typevar::TypeVar)
            typevar.lb === Union{} && typevar.ub === Any
        end
        function union_subtypes_aux(type)
            if type.b isa Union
              [ type.a, union_subtypes(type.b)... ]
            else
              [ type.a, type.b ]
            end
        end
        function union_subtypes(type)
            map(strip_type_bounds, union_subtypes_aux(type))
        end
        # The function signature is annotated with superfluous type
        # bounds that we need to strip out in order for it to compile
        # again.
        # Example:
        #
        # f(x::X{T}) where {T <: Int}
        #
        # becomes
        #
        # f(x::X{T <: Int}) where {T <: Int}
        #
        # internally. Similarly, if Y has two type parameters, then
        #
        # Y{T} becomes Y{T, S} where S internally.
        function strip_type_bounds(type)
            if type isa TypeVar
                # basic case
                type.name
            elseif type isa UnionAll
                if Base.print_without_params(type)
                    # early return for types like RoundMode or Rational.
                    return type
                end
                # type is basic parameterized type represented in a
                # UnionAll linked list, with type.var denoting the
                # type parameter and type.body the tail of the list.
                # First, we unwrap it analogously to Base.unwrap_unionall(),
                # but preserving type parameters.
                vars = [ type.var ]
                type = type.body
                while type isa UnionAll
                    push!(vars, type.var)
                    type = type.body
                end
                if type === Union{}
                    # bottom type
                    type
                elseif type isa Union
                    # union type
                    :(Union{$(union_subtypes(type)...)}
                        where {$(map(bounded, vars)...)})
                else
                    # X{...} where X is not a union type.
                    parameters = [ type.parameters... ]
                    # strip off trailing parameters in types of the
                    # kind Y{T, S} where S.
                    while !isempty(parameters) && !isempty(vars) &&
                            last(parameters) == last(vars) &&
                            isunbounded(last(vars))
                        pop!(vars)
                        pop!(parameters)
                    end
                    # type.name is not necessarily a symbol, but a
                    # qualified type, e.g. Module.T, and we need to
                    # reconstruct it as an Expr().
                    typename = qualified_typename(type.name)
                    # Return type.
                    :($(typename){$(map(strip_type_bounds, parameters)...)}
                        where {$(map(bounded, vars)...)})
                end
            elseif type isa Union
                    # unparameterized union type
                    :(Union{$(union_subtypes(type)...)})
            elseif type isa DataType
                if length(type.parameters) == 0
                    type
                else
                    # parameterized type, but all parameters are concrete.
                    # we do not need to manage typevars of the type, but
                    # the overall function may still have typevars, so
                    # we still need to strip bounds recursively.
                    if type isa Union
                        :(Union{$(union_subtypes(type)...)})
                    else
                        typename = qualified_typename(type.name)
                        :($(typename){$(map(strip_type_bounds,
                            type.parameters)...)})
                    end
                end
            elseif type isa Symbol
                # This is a special case for Val{:lex} types, as
                # symbols are treated differently in quoted expressions.
                Meta.quot(type)
            else
                type
            end
        end
        for method in methods(fn).ms
            funcname = method.name
            # Functions that have type variables store them in a linked list
            # in their signature, where sig.var is the type variable and
            # sig.body is the link to the next item in the list. Such linked
            # lists are of type UnionAll. See Base.unwrap_unionall() and
            # related functions for similar functionality.
            typevars = []
            sig = method.sig
            while sig isa UnionAll
                typevar = sig.var
                push!(typevars, bounded(typevar))
                sig = sig.body
            end
            # After stripping out all the type variables, the remaining
            # signature contains the function type and the argument types.
            resulttype = sig.parameters[1]
            argtypes = []
            rawargtypes = []
            for argtype in sig.parameters[2:end]
                push!(rawargtypes, argtype)
                push!(argtypes, strip_type_bounds(argtype))
            end
            has_kwargs = length(Base.kwarg_decl(method)) > 0
            is_vararg = method.isva
            argnames = []
            args = []
            argc = 0
            for argtype in argtypes
                argname = Symbol("a$(argc)")
                argc += 1
                # There is no need to deal with optional positional
                # arguments; Julia will create separate methods for all
                # arities and supply the defaults in those methods.
                if is_vararg && argc == length(argtypes)
                    # If we are dealing with a vararg method, the last
                    # argument needs to be annotated with ... to pass on the
                    # varargs tuple properly.
                    push!(argnames, :($argname...))
                else
                    push!(argnames, argname)
                end
                push!(args, :($argname :: $(argtype)))
            end
            # TODO: We don't enumerate keyword args; they are processed
            # in the body of the wrapped method, anyway. However, this
            # impedes tab completion with using the REPL.
            #
            # The main challenge is extracting the name of the arguments,
            # as the names can contain additional type information.
            kwargs = []
            if has_kwargs
                push!(kwargs, :(kw...))
            end
            funcdef = quote
                function $funcname($(args...);$(kwargs...)) where {$(typevars...)}
                    try
                        Sync._lock()
                        $mod.$funcname($(argnames...);$(kwargs...))
                    finally
                        Sync._unlock()
                    end
                end
                $(if Base.isexported(mod, funcname)
                    :(export $funcname)
                end)
            end
            Sync.eval(funcdef)
        end
    end

    function _forward_value(mod :: Module, name :: Symbol)
        if !(Base.isbindingresolved(Sync, name))
            code = :( const $(name) = $(mod).$(name) )
            Sync.eval(code)
        end
    end

    function _wrap_module(mod :: Module)

        function find_def(name)
            if '#' in string(name) || name in _reserved
                return nothing
            end
            try
                result = Base.eval(mod, name)
            catch
                nothing
            end
        end


        errors = []
        unwrapped = Set{Symbol}()
        for name in names(mod; all = true)
            def = find_def(name)
            if def isa Function
                try
                    _wrap_methods(mod, def)
                catch
                    _forward_value(mod, name)
                    push!(errors, name)
                    # debugging code
                    # println("error wrapping $name[$method_num]")
                    # for (exception, backtrace) in Base.catch_stack()
                    #     showerror(stdout, exception, backtrace)
                    #     println()
                    # end
                end
            elseif def !== nothing
              _forward_value(mod, name)
            end
        end
    end

end

macro sync(expr)
    quote
        try
            Sync._lock()
            $(esc(expr))
        finally
            Sync._unlock()
        end
    end
end


@static if false
    eval(:(module WrapExample
        display(x::Int) = println("Int: $x")
        display(x::String) = println("String: $x")
        export display
        display2(;alpha = 1, kw...) = println("kw: alpha = $(alpha); $(kw)")
        display3(x::T) where {T <: Int} = println("varargs: $(x)")
        display4(x...) where Q = println(x)
    end))

    Sync._wrap_module(WrapExample)

    Sync.display(10)
    Sync.display("Hello, world!")
    Sync.display2(; alpha = 2, beta = 99)
    Sync.display3(10)
    Sync.display4(1, 2)
end
