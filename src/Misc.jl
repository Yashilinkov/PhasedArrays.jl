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
