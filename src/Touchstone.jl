
## Touchstone V1 only


# work in progress

struct NetworkParameters
    N_ports::Int64
    R_ref::Vector{Float64}
    frequencies::Vector{Float64}
    # frequency_prefix::String # GHZ, MHZ, KHZ, HZ
    params::Array{ComplexF64, 3}
    parameter_type::String # S or Z or Y or G or H 

    noise_freqiencies::Union{Nothing, Vector{Float64}}
    NF_min::Union{Nothing, Vector{Float64}}
    Gamma_opt::Union{Nothing, Vector{ComplexF64}}
    r_norm_noise::Union{Nothing, Vector{Float64}}
    metadata::Dict{String, Any}

end

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


function parse_touchstone(filename::String)::NetworkParameters

    #Calculate N_ports
    N_ports = get_N_ports(filename)

    ## find options and data
    data = []
    option_line = []
    found_option_line = false
    for line in split.(eachline(filename))

        if line[1] == "!" # do nothing, it is comment
        elseif line[1] == "#" && !found_option_line
            option_line = uppercase.(line)
        else 
            comm_idx = findfirst(isequal("!"),line)
            isnothing(comm_idx) ? push!(data,line) : push!(data,line[1:comm_idx-1])
        end
    end

    # Parse options
    frequency_prefix = intersect(option_line, ["GHZ", "MHZ", "KHZ", "HZ"])
    frequency_prefix = isempty(frequency_prefix) ? "GHZ" : frequency_prefix[1]
    if frequency_prefix == "GHZ"
        frequency_scale = 1e9
    elseif frequency_prefix == "MHZ"
        frequency_scale = 1e6
    elseif frequency_prefix == "KHZ"
        frequency_scale = 1e3
    elseif frequency_prefix == "HZ"
        frequency_scale = 1e0
    end
    parameter_type = intersect(option_line, ["S", "Z", "Y", "H","G"])
    parameter_type = isempty(parameter_type) ? "S" : parameter_type[1] 
    format = intersect(option_line, ["MA", "DB", "RI"])
    format = isempty(format) ? "MA" : format[1]

    R_ref_idx = findfirst(==("R"), option_line)
    R_ref = ones(N_ports) * (isnothing(R_ref_idx) ? 50.0 : parse(Float64, option_line[R_ref_idx + 1]))


    # declare fields
    frequencies = Float64[]
    noise_frequencies = Float64[]
    NF_min = Float64[]
    Gamma_opt = ComplexF64[]
    r_norm_noise = Float64[]
    params = Matrix{ComplexF64}[]
    block = []

    if N_ports == 1
        first_line_length = 3
    elseif N_ports == 3
        first_line_length = 7
    else
        first_line_length = 9
    end
    block_length = 2*N_ports^2

    # parse data 
    for (line_idx,line) in enumerate(data)
        # Frequency line
        if length(line) == first_line_length
            current_frequency = parse(Float64, line[1])
            parsed_data = parse.(Float64, line[2:end])
    
            if isempty(frequencies)
                push!(frequencies, current_frequency*frequency_scale)
            elseif current_frequency*frequency_scale > frequencies[end]
                push!(frequencies, current_frequency*frequency_scale)
            else
                display(current_frequency)
                display(frequencies[end])
                display(line_idx)
                error("Frequency decreased at line $line_idx: $current_frequency < $(frequencies[end])")
    
            end
    
            append!(block, parsed_data)
        
        
        elseif length(line) == 5 && N_ports == 2 # noise data
            current_noise_frequency = parse(Float64, line[1])
            if current_noise_frequency ≤ frequencies[end]
                append!(noise_frequencies,current_noise_frequency)
                append!(NF_min,parse(Float64, line[2]))
                append!(Gamma_opt,parse(Float64, line[3])*cis(parse(Float64, line[4])))
                append!(r_norm_noise,parse(Float64, line[5]))
            end
    
        else # Continuation line
            append!(block, parse.(Float64, line))
        end
    
        # Once we have a full block, process it
        if length(block) == block_length
            tmp_mtrx = if format == "MA"
                reshape([block[i] * cis(block[i+1] * π / 180) for i in 1:2:block_length-1], N_ports, N_ports)
            elseif format == "RI"
                reshape([block[i] + im*block[i+1] for i in 1:2:block_length-1],N_ports,N_ports)
            elseif format == "DB"
                reshape([10^(block[i]/20) * cis(block[i+1] * π / 180) for i in 1:2:block_length-1], N_ports, N_ports)
            else
                error("Format $format is not supported")
            end
    
            push!(params, N_ports == 2 ? transpose(tmp_mtrx) : tmp_mtrx)
            empty!(block)
        end
    end
    params_3d = cat(params...; dims=3)

    return NetworkParameters(N_ports, 
                                R_ref,
                                frequencies,
                                params_3d,
                                parameter_type,
                                noise_frequencies,
                                NF_min,
                                Gamma_opt,
                                r_norm_noise,
                                Dict())



end




function s2z(S::Matrix,R_ref::Float64)
    G = I * np.R_ref
    F = I * (1/(2*sqrt(abs(np.R_ref))))
    return F^-1*(I-S)^-1 * (S*G+conj(G))*F
end

function z2s(Z::Matrix,R_ref::Float64)
    G = I * np.R_ref
    F = I * (1/(2*sqrt(abs(np.R_ref))))
    return F*(Z-conj(G))*(Z+G)^-1*F^-1
end