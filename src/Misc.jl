##
##
##  import from various sources
##
##

function parse_mdf(filename::String)
    reading_block = false
    data = Dict{String, Any}()
    block_meta = Dict{String, Any}()
    current_meta = Dict{String, Any}()
    num_lines = 0
    num_entries = 0
    block_name = String[]
    vars = []
    block_data = []

    file = open(filename)
    for line in readlines(file)
    if startswith(line,"!") || isempty(line) 
        # continue
    else
        tokens = split(strip(line))

        if tokens[1] == "VAR"
            var_name_type = tokens[2]
            var_value = parse(Float64, tokens[4])  # tokens[3] is "="
            if endswith(var_name_type, "(real)")
                block_meta[var_name_type] = var_value
            elseif endswith(var_name_type, "(complex)")  # not seen, but handle just in case
                error("Complex VAR not supported yet")
            elseif endswith(var_name_type, "(integer)")
                block_meta[var_name_type] = Int(var_value)
            else
                error("Unknown VAR type: $var_name_type")
            end
        

        elseif  ~reading_block && tokens[1] == "BEGIN"
            reading_block = true
            block_name = tokens[2]
            num_lines = 0
            num_entries = 0
            vars = []
            block_data = []
            current_meta = copy(block_meta)  # store associated VARs
            empty!(block_meta)  # clear for next block

        elseif reading_block && tokens[1] == "%"
            num_lines += 1
            append!(vars,tokens[2:end])

        elseif reading_block && tokens[1] ∉ ["%", "END", "BEGIN"]
            if num_entries == 0
                for var in vars
                    if occursin("(complex)",var)
                        num_entries += 2
                    elseif occursin("(real)",var)
                        num_entries += 1
                    elseif occursin("(integer)",var)
                        num_entries += 1
                    else
                        error("Unknown variable type: $var")
                    end
                end
            end
            append!(block_data,parse.(Float64,tokens))

        elseif reading_block && tokens[1] == "END"
            reading_block = false
    
            if length(block_data)%num_entries ≠ 0
                error("wrong number of data entries")
            end

            n_rows = length(block_data) ÷ num_entries
            matrix = transpose(reshape(block_data, (num_entries, n_rows)))
            col_data = Dict{String, Any}()
            i = 1
            for var in vars
                if endswith(var, "(real)")
                    col_data[var] = matrix[:, i]
                    i += 1
                elseif endswith(var, "(integer)")
                    col_data[var] = matrix[:, i]
                    i += 1
                elseif endswith(var, "(complex)")
                    real_part = matrix[:, i]
                    imag_part = matrix[:, i+1]
                    col_data[var] = complex.(real_part, imag_part)
                    i += 2
                else
                    error("Unknown var type")
                end
            end
            # data[block_name] = col_data
            if isnothing(current_meta)
                data[block_name] = (
                    vars = vars,
                    values = col_data
                )
            else
                if !haskey(data, block_name)
                    data[block_name] = []
                end

                push!(data[block_name], (
                    vars = copy(vars),
                    values = col_data,
                    meta = current_meta
                ))
            end
        end

    end
    end
    return data
end




function parse_ffs(filename::String)
    lines = readlines(filename)

    meta = Dict{String, Any}()
    data = Dict{String, Any}()
    data["phi"] = Float64[]
    data["theta"] = Float64[]
    data["E_Theta"] = ComplexF64[]
    data["E_Phi"] = ComplexF64[]
    num_entries = Inf
    vars = []
    reading_data = false

    phi = Float64[]
    theta = Float64[]
    E_Theta = ComplexF64[]
    E_Phi = ComplexF64[]

    for (i,line) in enumerate(strip.(lines))
        if occursin("Version", line)
            meta["version"] = strip(lines[i+1])

        elseif occursin("Data Type", line)
            meta["data_type"] = strip(lines[i+1])

        elseif occursin("#Frequencies", line)
            meta["num_frequencies"] = parse(Int, strip(lines[i+1]))

        elseif occursin("Position", line)
            meta["position"] = parse.(Float64, split(strip(lines[i+1])))

        elseif occursin("zAxis", line)
            meta["z_axis"] = parse.(Float64, split(strip(lines[i+1])))

        elseif occursin("xAxis", line)
            meta["x_axis"] = parse.(Float64, split(strip(lines[i+1])))

        elseif occursin("Radiated/Accepted/Stimulated Power", line)
            meta["power_radiated"] = parse(Float64, strip(lines[i+1]))
            meta["power_accepted"] = parse(Float64, strip(lines[i+2]))
            meta["power_stimulated"] = parse(Float64, strip(lines[i+3]))
            meta["frequency"] = parse(Float64, strip(lines[i+4]))

        elseif occursin(">> Total #phi samples", line)
            nums = split(strip(lines[i+1]))
            meta["phi_samples"] = parse(Int, nums[end - 1])
            meta["theta_samples"] = parse(Int, nums[end])
        elseif occursin(">> Phi, Theta,",line)
            vars = split(line)[3:end]
            for var in vars
                data[var] = []
            end
            num_entries = length(vars)
            reading_data = true
        end

        if reading_data == true && length(split(line)) == num_entries
            vals = parse.(Float64,split(line))

            # hardcoded order, sorry
            push!(data["phi"],vals[1])
            push!(data["theta"],vals[2])
            push!(data["E_Theta"],vals[3]+im*vals[4])
            push!(data["E_Phi"],vals[5]+im*vals[6])

        end


    end
    return data
end

