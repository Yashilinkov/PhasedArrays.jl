#######################
##                   ##
##      pattern      ##
##      plots        ##
##                   ##
#######################

#=
struct Pattern
theta
phi
Υθ
Υϕ
Υ
end
=#

@enum PatternType begin
    PatternLin
    PatterndB
    DirectivityLin
    DirectivitydB
    PowerLin
    PowerdB
end

function get_pattern_values(patt::Pattern,pattern_type::PatternType)
    if pattern_type == PatternLin
        return patt.Υ
    elseif pattern_type == PatterndB
        return 20*log10.(abs.(patt.Υ))
    elseif pattern_type == PowerLin
        return get_power_pattern(patt)
    elseif pattern_type == PowerdB
        return 10*log10.(get_power_pattern(patt))
    elseif pattern_type == DirectivityLin
        return get_directivity(patt)
    elseif pattern_type == DirectivitydB
        return 20*log10(get_directivity(patt))
    else
        error("Unsupported PatternType: $(pattern_type)")
    end

    
end

function plot_pattern_cuts(patt::Pattern,patter_type::PatternType; theta=nothing,phi=[0.0,90.0])

    values = get_pattern_values(patt,patter_type)

    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Angle [°]", ylabel="Pattern")
    if isnothing(theta)
        ϕ_idx = [findfirst(isapprox.(patt.phi,ϕ)) for ϕ in phi]
        θ_vals = patt.theta
        curves = [lines!(ax, θ_vals, values[:, idx], label="ϕ = $(phi[i])°") for (i, idx) in enumerate(ϕ_idx) if idx !== nothing]
    elseif !isnothing(theta) 
        θ_idx = [findfirst(isapprox.(patt.phi,θ)) for θ in theta]
        ϕ_vals = patt.phi
        curves = [lines!(ax, ϕ_vals, values[idx, :], label="θ = $(theta[i])°") for (i, idx) in enumerate(θ_idx) if idx !== nothing]
    end
    
    # axislegend(ax)
    legend = Legend(fig[1, 2], ax, framevisible = false)
    return fig

end






function plot_pattern_uv(patt::Pattern, patter_type::PatternType)
    values = get_pattern_values(patt, patter_type)
    theta = patt.theta
    phi = patt.phi
    upper_mask = theta .≤ 90.0
    lower_mask = theta .≥ 90.0
    theta_u = theta[upper_mask]
    theta_l = theta[lower_mask]
    values_u = values[upper_mask, :]
    values_l = values[lower_mask, :]
    u_u = zeros(Float64,length(theta_u),length(phi))
    u_l = zeros(Float64,length(theta_l),length(phi))
    v_u = zeros(Float64,length(theta_u),length(phi))
    v_l = zeros(Float64,length(theta_l),length(phi))
    for (i,θ) in enumerate(theta_u)
        for (j,ϕ) in enumerate(phi)
            u_u[i,j] = sind(θ)*cosd(ϕ)
            v_u[i,j] = sind(θ)*sind(ϕ)
        end
    end
    for (i,θ) in enumerate(theta_l)
        for (j,ϕ) in enumerate(phi)
            u_l[i,j] = sind(θ)*cosd(ϕ)
            v_l[i,j] = sind(θ)*sind(ϕ)
        end
    end
    
    fig = Figure()
    ax1 = Axis(fig[1, 1],
        xlabel = "u", ylabel = "v",
        xgridvisible = false, ygridvisible = false,
        backgroundcolor = :white
    )
    ax2 = Axis(fig[1, 2],
        xlabel = "u", ylabel = "v",
        xgridvisible = false, ygridvisible = false,
        backgroundcolor = :white
    )
    
    surface!(ax1,u_u,v_u,values_u;shading = false)
    surface!(ax2,u_l,v_l,values_l;shading = false)
    return fig
end