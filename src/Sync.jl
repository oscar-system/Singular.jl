module Sync

    const _mutex = ReentrantLock()

    function _lock()
        Base.lock(_mutex)
    end

    function _unlock()
        Base.unlock(_mutex)
    end

    function _wrap_methods(mod::Module, fn::Function)
        for method in methods(fn).ms
            funcname = method.name
            # Functions that have type variables store them in a linked list
            # in their signature, where sig.var is the type variable and
            # sig.body is the link to the next item in the list. Such linked
            # lists are of type UnionAll.
            typevars = []
            sig = method.sig
            while sig isa UnionAll
                typevar = sig.var
                push!(typevars, :( $(typevar.lb) <: $(typevar.name) <: $(typevar.ub) ))
                sig = sig.body
            end
            # After stripping out all the type variables, the remaining
            # signature contains the function type and the argument types.
            resulttype = sig.parameters[1]
            argtypes = []
            for argtype in sig.parameters[2:end]
                # We need to strip type bounds from types involving type
                # variables, as they are already part of the where clause.
                #
                # TODO: This does not handle the case T{A<:Bound} correctly,
                # only T<:Bound.
                #
                # Parameterized types with type variables are again
                # represented as UnionAll linked lists.
                if argtype isa TypeVar
                    push!(argtypes, argtype.name)
                else
                    push!(argtypes, argtype)
                end
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
                push!(args, :($argname :: $argtype))
            end
            # We don't enumerate keyword args; they are processed in the
            # body of the wrapped method, anyway.
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

    function _wrap_module(mod :: Module)
        reserved = Set([:eval, :include, :__init__])

        function find_def(name)
            if '#' in string(name) || name in reserved
                return nothing
            end
            try
                result = Base.eval(mod, name)
            catch
                nothing
            end
        end


        for name in names(mod; all = true)
            def = find_def(name)
            if def isa Function
                try
                    _wrap_methods(mod, def)
                catch
                    # debugging code
                    # println("error wrapping $name")
                    for (exception, backtrace) in Base.catch_stack()
                        # showerror(stdout, exception, backtrace)
                        # println()
                    end
                end
            end
        end
    end

end

@static if false
    eval(:(module WrapExample
        display(x::Int) = println("Int: $x")
        display(x::String) = println("String: $x")
        export display
        # f(x::T, Y::U) where {T, U} = 0
        display2(;alpha = 1, kw...) = println("kw: alpha = $(alpha); $(kw)")
        display3(x::T) where {T <: Int} = println("varargs: $(x)")
        display4(x...) = println(x)
    end))

    Sync._wrap_module(WrapExample)

    Sync.display(10)
    Sync.display("Hello, world!")
    Sync.display2(; alpha = 2, beta = 99)
    Sync.display3(10)
    Sync.display4(1, 2)
end
