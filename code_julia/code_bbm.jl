# Install packages
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


# Load packages
using DelimitedFiles
using LinearAlgebra

using SummationByPartsOperators
using Plots: Plots, plot, plot!, savefig, gr


# Physical setup of a traveling wave solution with speed `c()`
xmin() = -90.0
xmax() = -xmin()
c() = 1.2
function usol(t, x)
  A = 3 * (c() - 1)
  K = 0.5 * sqrt(1 - 1 / c())
  x_t = mod(x - c() * t - xmin(), xmax() - xmin()) + xmin()

  return A / cosh(K * x_t)^2
end


# Numerical setup of the semidiscretization with Fourier pseudospectral methods
function bbm_rhs_cubic(u, param, t)
  D1_invImD2 = param.D1_invImD2

  # this semidiscretization conserves the linear and cubic invariants
  tmp = @. -(0.5 * u^2 + u)

  return D1_invImD2 * tmp
end

function bbm_setup(nnodes)
  D1 = fourier_derivative_operator(xmin(), xmax(), nnodes)
  D1op = D1
  D2 = D1^2
  ImD2 = I - D2
  D1_invImD2 = D1 / (I - D2)

  tspan = (0.0, 40 * (xmax() - xmin()) / c())
  x = grid(D1)
  u0 = usol.(tspan[1], x)

  param = (; D1op, D1, D2, ImD2, D1_invImD2, usol)

  return u0, param, tspan
end


# Functionals to save during the time integration
function bbm_functionals(u, params, t)
  #(; D2) = params
  D2 = params.D2
  D1 = params.D1op
  x = grid(D1)

  linear = integrate(u, D1)

  tmp1 = D2 * u
  @. tmp1 = u^2 - u * tmp1
  quadratic = integrate(tmp1, D1)

  @. tmp1 = (u + 1)^3
  cubic = integrate(tmp1, D1)

  @. tmp1 = u - usol(t, x)
  error_l2 = integrate(abs2, tmp1, D1) |> sqrt

  return (; linear, quadratic, cubic, error_l2)
end



# The main time integration loop
function bbm_solver(; rhs, alg, nnodes = 2^6, dt = 0.25)
  u0, params, tspan = bbm_setup(nnodes)

  time      = Vector{Float64}()
  linear    = Vector{Float64}()
  quadratic = Vector{Float64}()
  cubic     = Vector{Float64}()
  error_l2  = Vector{Float64}()

  t = first(tspan)
  u = copy(u0)
  unew = similar(u)
  step = 0
  while t < last(tspan)
    step += 1

    if alg === :midpoint
      k1 = rhs(u, params, t)
      y2 = @. u + (dt / 2) * k1
      k2 = rhs(y2, params, t + dt / 2)
      @. unew = u + dt * k2
    elseif alg === :ep5
      k1 = rhs(u, params, t)
      y2 = @. u + (1/10) * dt * k1
      k2 = rhs(y2, params, t + (1/10) * dt)
      y3 = @. u + (56795/35721) * dt * k2 - (35816/35721) * dt * k1
      k3 = rhs(y3, params, t + (19/20) * dt)
      y4 = @. u + (215846127/181744000) * dt * k3 + (11994761/5328000)* dt * k1 - (11002961/4420800) * dt * k2
      k4 = rhs(y4, params, t + (37/63) * dt)
      @. unew = u + (dt) * ((-17/222) * k1 + (6250/15657) * k2 + (5250987/10382126) * k3 + (4000/23307) * k4)
    end

    t += dt
    u .= unew

    if (step % 100 == 0) && (t >= 1)
      functionals = bbm_functionals(u, params, t)
      push!(time,      t)
      push!(linear,    functionals.linear)
      push!(quadratic, functionals.quadratic)
      push!(cubic,     functionals.cubic)
      push!(error_l2,  functionals.error_l2)
    end
  end

  x = grid(params.D1op)
  return (; time, linear, quadratic, cubic, error_l2, x, u, u0)
end

function plot_kwargs()
  fontsizes = (
    xtickfontsize = 14, ytickfontsize = 14,
    xguidefontsize = 16, yguidefontsize = 16,
    legendfontsize = 14)
  return (; linewidth = 4, gridlinewidth = 2, fontsizes...)
end

# Code to generate the plots in function of a Boolean 'pep'.
# This returns the arrays for time steps and the error and saves
# them in a text file.
# By default, we set it to false, so that the function
# 'generate_arrays' works with the midpoint method.

function generate_arrays(pep = false)
  if pep == true
    res = bbm_solver(; nnodes = 2^6, dt = 0.2, rhs = bbm_rhs_cubic, alg = :ep5)
  else
    res = bbm_solver(; nnodes = 2^6, dt = 0.1, rhs = bbm_rhs_cubic, alg = :midpoint)
  end
  # Save the error and time arrays only.
  error = res.error_l2
  t_array = res.time

  # Save the arrays in a text file (because we want to make the plots
  # in Python).
  file_path = joinpath(@__DIR__, "arrays_data_$(pep).txt")
  open(file_path, "w") do file
    # Export arrays to the text file
    writedlm(file, [t_array error])
  end

  return t_array, error, file_path
end

# This function calls the previous one to make the plot.
function bbm_plot(pep = false)
  t, error, _ = generate_arrays(pep)

  if pep == true
    plt = plot(t, error, label = "PEP Method")
  else
    plt = plot(t, error, label = "Explicit-Midpoint")
  end
  # Show the plot
  display(plt)
end


# Main function to generate all data
function main()
  fig = plot(; xlabel = "t", ylabel = "Error", plot_kwargs()...,
               title = "Error with respect to time", legend = :topleft,
               xscale = :log10, yscale = :log10)

  t, error, file_path = generate_arrays(true)
  plot!(fig, t, error; label = "PEP(4, 2, 5)", color = "red", plot_kwargs()...)
  @info "Data using the PEP method saved in $(file_path)"

  t, error, file_path = generate_arrays(false)
  plot!(fig, t, error; label = "RK(2, 2)", color = "blue", plot_kwargs()...)
  @info "Data using the standard method saved in $(file_path)"

  return fig
end
