"""
    getTable(params::Dict,tablename::String;tableOptional::Bool=false)

Returns the table named `tablename` from the contents of a TOML file, which was
previously parsed into the dictionary `params`. If `tableOptional` is `true` 
the output will be either a dictionary with the contents of the table or 
`missing`. If `tableOptional` is true the function will throw an error instead
of returning missing.
"""
function getTable(params::Dict,tablename::String;tableOptional::Bool=false)
    for part in split(tablename,".")
        params = get(params,part,missing)
        ismissing(params) & !tableOptional && throw(MissingException("Non-optional table [$(tablename)] missing in configuration file."))
    end
    return params
end

"""
    getValue(params::Dict,tablename::String,key::Sring;valueOptional::Bool=false)

Returns a the value specified by a `key` from the table named `tablename` from
the contents of a TOML file, which was previously parsed into the dictionary 
`params`. If `valueOptional` is `true` the output will be either the value or 
`missing`. If `valueOptional` is true the function will throw an error instead
of returning missing.
"""
function getValue(params::Dict,tablename::String,key::String;valueOptional::Bool=false)
    params = getTable(params,tablename,tableOptional=valueOptional)
    if ismissing(params)
	    return missing
    else
	    value = get(params,key,missing)
        ismissing(value) & !valueOptional && throw(MissingException("key $(key) expected in table [$(tablename)] but missing."))
	    return value
    end
end

"""
    getType(params::Dict,tablename::String)

Tries to parse the key "type" from the table named `tablename` from the 
contents of a TOML file, which was previously parsed into the dictionary 
`params` into an actual julia type. Returns missing if it fails to do so, e.g.
if the key is missing in the TOML file.
"""
getType(params::Dict,tablename::String) = eval(Symbol(getValue(params,tablename,"type",valueOptional=true)))

"""
    initialize(T::Type,params::Dict,floatType::Type)

Tries to initialize the type `T` from the contents of a TOML file, which was
previously parsed into the dictionary `params`. To this end all values obtained
from the TOML are converted to the type `floatType` or arrays thereof. So `T`
needs to have a suitable constructor.
"""
function initialize(T::Type,tablename::String,params::Dict,floatType::Type)
    fieldNamesT = map(string,fieldnames(T))
    values = map(key->getValue(params,tablename,key),fieldNamesT)
    fieldValues = []
    for (i,value) in enumerate(values)
        TField = fieldtype(T,i)
        if TField <: Integer && isprimitivetype(TField)
            push!(fieldValues,TField(value))
        elseif typeof(value) <: Real
            push!(fieldValues,floatType(value))
        elseif typeof(value) <: Array{T} where {T<:Array}
            push!(fieldValues,floatType.(hcat(value...)))
        elseif typeof(value) <: Array{T} where {T<:Real}
            push!(fieldValues,floatType.(value))
        elseif typeof(value) <: Dict{String,Any}
            if length(keys(value))==1
                typeName = collect(keys(value))[1]
                type = eval(Symbol(typeName))
                push!(fieldValues,initialize(type,typeName,value,floatType))
            else
                throw(ErrorException("Recursive initialization failed. Dictionariy $value expected to have a single key only."))
            end
        else
            throw(ErrorException("I do not know how to handle value $value when trying to initialize type $T from table $tablename"))
        end
    end
    return T(fieldValues...)
end

"""
    getTableName(T::Type)

Outputs the tablename which is expected to hold the data to `initialize` the
type `T` from the contents of a TOML file.
"""
function getTableName(T::Type)
    out = String[]
    while T != Any
        typeString = string(T)
        m = match(r"^(?:\w*\.)*(\w+)", typeString)
	    pushfirst!(out,m.captures[1])
        T = supertype(T)
    end
    return join(out,".")
end 
