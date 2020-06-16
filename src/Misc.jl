"""
	@mustimplement f(x::T,args...)

Macro helping to defining a public interface. Let `T` be an abstract type then
`@mustimplement f(x::T,...)` will issue a meaningfull error, if a specific
subtype `S<:T` is implemented, but the function `f(x::S,...)` is not.
"""
macro mustimplement(sig)
    fname = sig.args[1]
    arg1 = sig.args[2]
    while isa(fname,Expr)
        arg1 = fname.args[2]
        fname = fname.args[1]
    end
    while isa(arg1,Expr)
        arg1 = arg1.args[1]
    end
    @info "" arg1
    :($(esc(sig)) = error($(Expr(:quote,arg1)),"::",typeof($(esc(arg1)))," must implement ", $(Expr(:quote,sig))))
end
