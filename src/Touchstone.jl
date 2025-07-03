
## Touchstone V1 only
## MA only

# work in progress


function get_N_ports(filename::String)
        
    ext = splitext(uppercase( filename))[2]
    m = match(r"\.S(\d+)P", ext)

    if m !== nothing
        N = parse(Int, m.captures[1])
        println("Number of ports: ", N)
    else
        error("Could not determine number of ports from extension: $ext")
    end
    return N
end


function parse_touchstone(filename::String)
    N = get_N_ports(filename)

lines = uppercase.(readlines(filename))
# remove coments and empty lines 
filter!(line -> !isempty(strip(line)),lines)
comments = filter(line -> startswith(strip(line),"!"),lines)

lines_no_comments = filter(line -> !startswith(strip(line),"!"),lines)

# Option line
option_line = filter(line -> startswith(strip(line),"#"),lines_no_comments)

tokens = split(option_line...)

# frequency
frequency_multipliers = [1e9,1e6,1e3,1]
frequency_prefixes = ["GHZ","MHZ","KHZ","HZ"]
default_frequency_prefix = "GHZ"
default_frequency_multiplier = 1e9

freq_mult_idx = findfirst([prfx in tokens for prfx in frequency_prefixes])
if isnothing(freq_mult_idx)
    frequency_multiplier = default_frequency_multiplier
    frequency_prefix = default_frequency_prefix
else
    frequency_multiplier = frequency_multipliers[freq_mult_idx]
    frequency_prefix = frequency_prefixes[freq_mult_idx]
end

# parameter
parameter_types = ["S","Y","Z","H","G"]
default_parameter_type = "S"
par_type_idx = findfirst([pt in tokens for pt in parameter_types])

parameter_type = isnothing(par_type_idx) ? default_parameter_type : parameter_types[par_type_idx]

# format
formats = ["DB","MA","RI"]
default_format = "MA"

fmt_idx = findfirst([fmt in tokens for fmt in formats])
format = isnothing(fmt_idx) ? default_format : formats[fmt_idx]

# Ref impedance
default_ref_resistance = 50
R_idx = findfirst(tokens .== "R")
ref_R = isnothing(R_idx) ? default_ref_resistance : parse(Float64,tokens[R_idx+1])

data_lines = filter(line -> !startswith(strip(line),"#"),lines_no_comments)
data = split.(strip.(data_lines))

# if N == 2 then order is F, S11, S21, S12, S22 !!!!!!
# if not, then order is s11 s12 s13 ... S21 S22 S23 ... and so on

block_length = (N*N*2+1)
N_entries = sum(length.(data))

if N_entries%block_length == 0
    N_freqs = N_entries÷block_length
else
    error("simething is wrong with data, expect N_entries%blocklength == 0")
end


freqs = zeros(Float64,N_freqs)
params = []

block = []

i = 1
while !isempty(data)
    append!(block,popfirst!(data))
    if length(block) == block_length
        freqs[i] = parse(Float64,popfirst!(block))
        if format == "MA"
            mag_angle = parse.(Float64,block)
            params_tmp = ComplexF64[]

            for j in 1:2:length(mag_angle)
                append!(params_tmp,mag_angle[j]*cis(mag_angle[j+1]/180*π))
            end
            params_tmp = reshape(params_tmp,N,N)
            if N > 2
                params_tmp=Matrix(transpose(params_tmp))
            end
            push!(params,params_tmp)
        else
            error("$format in not implemented yet")
        end
        i += 1
        empty!(block)
    end
end
return params, freqs
end

