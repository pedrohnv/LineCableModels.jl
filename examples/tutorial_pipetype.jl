#=
# Tutorial 2.1 - Calculating a cable's electrical parameters

This tutorial calculates the series impedance and shunt admittance of the cable as a function of the frequency. Then the modal decomposition is performed.
=#

using LineCableModels
using Plots

steel = Material(13.8e-8, 1.0, 300.0, 20.0, 0.00450)
copper = Material(1.835e-8, 1.0, 0.999994, 20.0, 0.00393)
insulator1 = Material(1e15, 3.31, 1.0, 20.0, 0.005)
insulator2 = Material(1e15, 2.3, 1.0, 20.0, 0.005)
insulator3 = Material(1e15, 10.0, 1.0, 20.0, 0.005)

rc = radii_single_core = [0.0, 9.6, 17.054, 18.054, 19.50] .* 1e-3
ra = radii_armor = [48.0, 59.0, 65.0] .* 1e-3

# Single-core jacket
core = ConductorGroup(Tubular(rc[1], rc[2], copper))
core_insulation = InsulatorGroup(Insulator(rc[2], rc[3], insulator1))
core_cc = CableComponent("core", core, core_insulation)
jacket_design = CableDesign("Jacket", core_cc)

sheath = ConductorGroup(Tubular(rc[3], rc[4], copper))
sheath_insulation = InsulatorGroup(Insulator(rc[4], rc[5], insulator2))
sheath_cc = CableComponent("sheath", sheath, sheath_insulation)
add!(jacket_design, sheath_cc)

# Pipe (armor)
pipe_armor = ConductorGroup(Tubular(ra[1], ra[2], steel))
pipe_sheath = InsulatorGroup(Insulator(ra[2], ra[3], insulator3))
pipe_cc = CableComponent("pipe", pipe_armor, pipe_sheath)
pipe_design = CableDesign("Pipe", pipe_cc)

# Cable system
x0, y0 = 0.0, -0.10
# xa, ya, xb, yb, xc, yc = trifoil_formation(x0, y0, rc[end])  # FIXME is this correct?
di = rc[end] / cos(deg2rad(30))
xa = x0
ya = y0 + di * sin(deg2rad(90))

xb = x0 + di * cos(deg2rad(210))
yb = y0 + di * sin(deg2rad(210))

xc = x0 + di * cos(deg2rad(330))
yc = y0 + di * sin(deg2rad(330))

# Initialize the `LineCableSystem` with the first cable (jacket A):
cablepos = CablePosition(jacket_design, xa, ya, Dict("core" => 1, "sheath" => 2))
cable_system = LineCableSystem("pipe_type_trifoil", 1000.0, cablepos)

# Add remaining jackets (B and C):
add!(cable_system, jacket_design, xb, yb, Dict("core" => 3, "sheath" => 4))
add!(cable_system, jacket_design, xc, yc, Dict("core" => 5, "sheath" => 6))

# Add the pipe surrounding the three cables:
add!(cable_system, pipe_design, x0, y0, Dict("pipe" => 7))

# Visualize the cross-section of the cable system:
cable_preview = preview(cable_system, zoom_factor=0.020)

# Frequency range [Hz]
nf = 100  # number of frequency points
freq = exp10.(range(0, 9, nf))  # logspace
# on 1 GHz the quasi-TEM assumption for transmission lines may no longer be valid, but we will calculate it nonetheless

# The computation functions expect the angular frequency in rad/s. This allows to easily use them with numerical Laplace transforms.
complex_frequencies = 2im * pi * freq  # jω

# extract the parameters from the cable system
rho = [copper.rho, copper.rho, steel.rho]
mu_r = [copper.mu_r, copper.mu_r, steel.mu_r]
epsilon_r = [insulator1.eps_r, insulator2.eps_r, insulator3.eps_r]
center_distances = fill(radii_single_core[end] / cos(deg2rad(30)), 3)
theta = deg2rad.([0, 120, -120])
loss_factor = zeros(3)

# Compute the series impedance Z and shunt admittance Y matrices
Y = comp_pipe_admittance(
    complex_frequencies,
    radii_single_core,
    radii_armor,
    center_distances,
    theta,
    epsilon_r,
    loss_factor
)

Z = comp_pipe_impedance(
    complex_frequencies,
    radii_single_core,
    radii_armor,
    center_distances,
    theta,
    rho,
    mu_r,
)

# Plot Z and Y
zplot = plot(
    freq,
    [abs.(Z[1, 1, :]), abs.(Z[1, 2, :]), abs.(Z[1, 3, :]), abs.(Z[1, 4, :]), abs.(Z[1, 7, :])],
    label = ["core A - core A"  "core A - sheath A"  "core A - core B"  "core A - sheath B"  "core A - armor"],
    xlabel = "Frequency [Hz]",
    ylabel = "Impedance Magnitude [Ω/m]",
    legend = :topleft,
    axis = :log10,
    minorgrid = true,
    xticks = exp10.(0:9),
    yticks = exp10.(-4:3),
)

yplot = plot(
    freq,
    [abs.(Y[1, 1, :]), abs.(Y[1, 2, :]), abs.(Y[1, 3, :]), abs.(Y[1, 4, :]), abs.(Z[1, 7, :])],
    label = ["core A - core A"  "core A - sheath A"  "core A - core B"  "core A - sheath B"  "core A - armor"],
    xlabel = "Frequency [Hz]",
    ylabel = "Admittance Magnitude [S/m]",
    legend = :right,
    axis = :log10,
    majorgrid = true,
    xticks = exp10.(0:9),
    xminorticks = true,
    yticks = [1e-21, 1e-9, 1e-3, 1e0],
)

# Modal Decomposition
# note, velocityy and attenuation are calculated from the propagation constant.
propagation, velocity, attenuation = calc_propagation_modes(Z, Y, complex_frequencies)

vel_plot = plot(
    freq,
    velocity,
    label = "",
    xlabel = "Frequency [Hz]",
    ylabel = "Modal Velocities [m/μs]",
    xaxis = :log10,
    minorgrid = true,
    xticks = exp10.(0:9),
)

att_plot = plot(
    freq,
    attenuation,
    label = "",
    xlabel = "Frequency [Hz]",
    ylabel = "Modal Attenuation [dB/m]",
    axis = :log10,
    minorgrid = true,
    xticks = exp10.(0:9),
    yticks = exp10.(-8:-1),
)
