"""
	LineCableModels.DataModel.PipeType

The [`PipeType`](@ref) submodule provides functions for determining the longitudinal impedance Z and transversal admittance Y per unit length in frequency-domain of Single-Core and Pipe-Type cables.

# Overview
- Calculation of series impedance (Z) and shunt admittance (Y) per unit length for single-core and pipe-type cables, including skin effect and dielectric losses.
- Computation of impedance and Maxwell potential coefficient matrices for both the internal conductors and the metallic sheath (armor) of pipe-type cables.
- Assembly of full multi-phase impedance and admittance matrices for cables with multiple single-core conductors inside a common pipe.
- Calculation of equivalent nodal admittance matrices for transmission line sections.
- Eigenvalue-based analysis of propagation modes, including velocity and attenuation, for multi-conductor systems.
- Utility functions for eigenvector rotation and Newton-Raphson refinement to ensure robust modal decomposition.
- These functions enable accurate frequency-dependent modeling of complex cable systems for electromagnetic transient and steady-state studies.

# Dependencies

(IMPORTS)

# Exports

(EXPORTS)
"""
module PipeType

# Load common dependencies
include("CommonDeps.jl")
using ...Utils

const μ₀ = 4π * 1e-7
const ε₀ = 8.8541878128e-12

# Module-specific dependencies
using SpecialFunctions
using LinearAlgebra

# ==============================================================================
# Pipe-Type related functions
# ==============================================================================

# TODO should we bother to adapt the formulas for N conductive sheaths? See, for example:
# Ghosh, S., & Das, S.K. (2021). A generalized approach to deduce the impedance and admittance of underground cable containing N semiconductor screens. Engineering Science and Technology, an International Journal.

# TODO maybe we get a performance boost if we change all `AbstractVector` to `StaticArrays`, see: https://github.com/JuliaArrays/StaticArrays.jl

# FIXME there is no sanity check if the `center_distances` is consistent with the given radii. If the cables overlap, the results will be nonsensical without throwing an error (silent bug).

"""
(TYPEDSIGNATURES)

Calculates the series impedance (`Z`) per unit length matrix of a coaxial cable compromised of a core conductor (hollow or solid), a dielectric layer surrounding the core and, optionally, a conductive sheath and second dieletric layer. Skin-effects are considered [ametani1980general](@cite).

# Arguments

- `complex_frequency`: complex angular frequency `s = c + jω` \\[rad/s\\].
- `radii`: array of the radius of each layer of the tubular conductor \\[m\\]. This must be either a length 3 or 5 array.
    `radii[1]`: inner radius of the core, use `0` if solid conductor.
    `radii[2]`: outer radius of the core.
    `radii[3]`: outer radius of the first insulation.
    `radii[4]`, Optional: outer radius of the conductive sheath.
    `radii[5]`, Optional: outer radius of the second insulation.
- `rho`: array of the conductivity of the conductor materials \\[dimensionless\\].
    `rho[1]`: core conductor conductivity.
    `rho[2]`, Optional: conductive sheath conductivity.
- `mu_r`: array of the Relative permeability of the conductor materials \\[dimensionless\\].
    `mu_r[1]`: core conductor permeability.
    `mu_r[2]`, Optional: conductive sheath permeability.

# Returns

- Series impedance matrix per unit length of the cable \\[Ω/m\\].

# Examples

```julia
Z = (FUNCTIONNAME)(...)
# Output: Z
```

# See also

- [`calc_tubular_resistance`](@ref)
"""
function calc_single_core_impedance(
    complex_frequency::Number,
    radii::AbstractVector{<:Number},
    rho::AbstractVector{<:Number},
    mu_r::AbstractVector{<:Number},
)
    r = radii
    # Some sanity checks
    if length(r) == 3
        if !(length(rho) == 1 && length(mu_r) == 1)
            throw(ArgumentError("For length(radii) == 3, rho and mu_r must have length 1."))
        end
    elseif length(r) == 5
        if !(length(rho) == 2 && length(mu_r) == 2)
            throw(ArgumentError("For length(radii) == 5, rho and mu_r must have length 2."))
        end
    else
        throw(ArgumentError("It is expected the radii to have length 3 or 5."))
    end

    jw = complex_frequency
    jwu_2pi = jw * μ₀ / (2π)
    eta_1 = sqrt((jw * μ₀ * mu_r[1]) / rho[1])
    w1 = eta_1 * r[1]
    w2 = eta_1 * r[2]
    if r[1] ≈ 0.0
        num = besselix(0, w2)
        den = besselix(1, w2)
        z1 = jw * μ₀ * mu_r[1] / (2π * w2) * num / den
    else
        scale = exp(abs(real(w1)) - abs(real(w2)) + w1 - w2)

        num1 = besselix(0, w2) * besselkx(1, w1)
        num2 = besselix(1, w1) * besselkx(0, w2)
        num = num1 + scale * num2

        den1 = besselix(1, w2) * besselkx(1, w1)
        den2 = besselix(1, w1) * besselkx(1, w2)
        den = den1 - scale * den2
    end
    z1 = jwu_2pi * mu_r[1] / w2 * num / den
    z2 = jwu_2pi * log(r[3] / r[2])  # external insulation

    if length(r) == 3  # no sheath
        Z = zeros(Complex, (1, 1))
        Z[1, 1] = z1 + z2
    elseif length(r) == 5  # with conductive sheath
        eta_2 = sqrt((jw * μ₀ * mu_r[2]) / rho[2])
        w1 = eta_2 * r[3]
        w2 = eta_2 * r[4]
        scale = exp(abs(real(w1)) - w2 - abs(real(w2)) + w1)

        den1 = besselix(1, w2) * besselkx(1, w1)
        den2 = besselix(1, w1) * besselkx(1, w2)
        den = den1 - scale * den2

        n1 = besselix(0, w1) * besselkx(1, w2)
        n2 = besselix(1, w2) * besselkx(0, w1)
        num_1 = scale * n1 + n2
        z3 = jwu_2pi * mu_r[2] / w1 * num_1 / den

        z4 = rho[2] * exp(w1 - abs(real(w2))) / (2π * r[3] * r[4] * den)

        n1 = besselix(0, w2) * besselkx(1, w1)
        n2 = besselix(1, w1) * besselkx(0, w2)
        num_2 = n1 + scale * n2
        z5 = jwu_2pi * mu_r[2] / w2 * num_2 / den

        z6 = jwu_2pi * log(r[5] / r[4])
        Z = zeros(Complex, (2, 2))
        Z[1, 1] = z1 + z2 + z3 + z5 + z6 - 2 * z4
        Z[1, 2] = z5 + z6 - z4
        Z[2, 1] = z5 + z6 - z4
        Z[2, 2] = z5 + z6
    else
        throw(ArgumentError("It is expected the radii to have length 3 or 5."))
    end
    return Z
end


"""
(TYPEDSIGNATURES)

Calculates the shunt Maxwell's potential coefficient (`P`) per unit length matrix of a coaxial cable compromised of a core conductor (hollow or solid), a dielectric layer surrounding the core and, optionally, a conductive sheath and second dieletric layer. Skin-effects are considered [ametani1980general](@cite).

# Arguments

- `complex_frequency`: complex angular frequency `s = c + jω` \\[rad/s\\].
- `radii`: array of the radius of each layer of the tubular conductor \\[m\\]. This must be either a length 3 or 5 array.
    `radii[1]`: inner radius of the core, use `0` if solid conductor.
    `radii[2]`: outer radius of the core.
    `radii[3]`: outer radius of the first insulation.
    `radii[4]`, Optional: outer radius of the conductive sheath.
    `radii[5]`, Optional: outer radius of the second insulation.
- `epsilon_r`: array of the Relative permittivity of the dieletric materials \\[dimensionless\\].
    `epsilon_r[1]`: first insulation permittivity.
    `epsilon_r[2]`, Optional: second insulation permittivity.
- `loss_factor`, Optional: array of the loss factor (tangent) of the dielectric materials \\[dimensionless\\].
    `loss_factor[1]`: first insulation loss factor.
    `loss_factor[2]`, Optional: second insulation loss factor.

# Returns

- Maxwell's potential coefficient matrix per unit length of the cable \\[V/c/m\\].

# Examples

```julia
P = (FUNCTIONNAME)(...)
# Output: P
```

# See also

- [`calc_tubular_resistance`](@ref)
"""
function calc_single_core_potential(
    radii::AbstractVector{<:Number},
    epsilon_r::AbstractVector{<:Number},
    loss_factor::AbstractVector{<:Number} = [0.0, 0.0],
)
    r = radii
    # Some sanity checks
    if length(r) == 3
        if !(length(epsilon_r) == 1 && length(loss_factor) == 1)
            throw(ArgumentError("For length(radii) == 3, epsilon_r and loss_factor must have length 1."))
        end
    elseif length(r) == 5
        if !(length(epsilon_r) == 2 && length(loss_factor) == 2)
            throw(ArgumentError("For length(radii) == 5, epsilon_r and loss_factor must have length 2."))
        end
    else
        throw(ArgumentError("It is expected the radii to have length 3 or 5."))
    end

    td1 = loss_factor[1]
    p1 = log(r[3] / r[2]) / (2π * epsilon_r[1] * ε₀) * (1 + 1im * td1) / (1 + td1^2)

    if length(r) == 3  # no sheath
        P = zeros(Complex, (1, 1))
        P[1, 1] = p1
    elseif length(r) == 5  # with conductive sheath
        td2 = loss_factor[2]
        p2 = log(r[5] / r[4]) / (2π * epsilon_r[2] * ε₀) * (1 + 1im * td2) / (1 + td2^2)
        P = zeros(Complex, (2, 2))
        P[1, 1] = p1 + p2
        P[1, 2] = p2
        P[2, 1] = p2
        P[2, 2] = p2
    else
        throw(ArgumentError("It is expected the radii to have length 3 or 5."))
    end
    return P
end


"""
(TYPEDSIGNATURES)

Calculates the impedance (`Z`) per unit length matrix for a Pipe-Type (PT) cable as a function of the metallic sheath.

The cables inside the pipe are considered to be identical single-core (SC), each composed of a central conductor, internal insulation, and optionally, a shield and external insulation.

# Arguments

- `complex_frequency`: Complex angular frequency `s = c + jω` \\[rad/s\\].
- `radii_single_core`: array of the radius of each layer of the SC cables \\[m\\]. This must be either a length 3 or 5 array.
    `radii_single_core[1]`: inner radius of the core, use `0` if solid conductor.
    `radii_single_core[2]`: outer radius of the core.
    `radii_single_core[3]`: outer radius of the first insulation.
    `radii_single_core[4]`, Optional: outer radius of the conductive sheath.
    `radii_single_core[5]`, Optional: outer radius of the second insulation.
- `radii_armor`: array of the armor radii \\[m\\]. Must have length 3.
    `radii_armor[1]`: inner radius of the armor.
    `radii_armor[2]`: outer radius of the armor.
    `radii_armor[3]`: outer radius of the external insulation surrounding the armor.
- `center_distances`: array (length `Nd`) of the distances from the center of the Pipe to the center of each SC cable \\[m\\].
- `theta`: array of the angles of each SC cable inside the Pipe \\[rad\\].
- `rho_armor`: Electrical resistivity of the armor \\[Ω·m\\].
- `mur_armor`: Relative magnetic permeability of the armor \\[dimensionless\\].

# Returns

- `Z`: Series impedance matrix \\[Ω/m\\], shape `(Nd * ncx + 1, Nd * ncx + 1)`.
"""
function calc_pipe_internal_impedance(
    complex_frequency::Number,
    radii_single_core::AbstractVector{<:Number},
    radii_armor::AbstractVector{<:Number},
    center_distances::AbstractVector{<:Number},
    theta::AbstractVector{<:Number},
    rho_armor::Number,
    mur_armor::Number,
)
    rc = radii_single_core
    ra = radii_armor
    di = center_distances
    rho_a = rho_armor
    mur_a = mur_armor
    jw = complex_frequency
    if length(radii_armor) != 3
        throw(ArgumentError("It is expected the radii_armor to have length 3 or 5."))
    end
    # ncx: number of "active" conductors in each single-core cable
    if length(rc) == 3
        ncx = 1
    elseif length(rc) == 5
        ncx = 2
    else
        throw(ArgumentError("radii_single_core must have length 3 or 5"))
    end
    re = rc[end]
    rp1 = ra[1]
    y1 = rp1 * sqrt((jw * μ₀ * mur_a) / rho_a)
    q = zeros(Complex, length(di), length(di))
    Ntrunc = 26  # number of terms to truncate the infinite series
    kve_n_y1 = [besselkx(n, y1) for n in 0:Ntrunc]
    kve_n_y1_ratio = kve_n_y1[1:end-1] ./ kve_n_y1[2:end]
    rp1_sq = rp1^2
    di_sq = di .^ 2

    # impedance
    for k in eachindex(di)
        for i in eachindex(di)
            theta_jk = theta[i] - theta[k]
            sum1 = 0.0
            sum2 = 0.0
            if k == i
                for n in 1:Ntrunc
                    sum1 += (di[k] / rp1)^(2n) / ((mur_a + 1) * n + y1 * kve_n_y1_ratio[n])
                end
                sum1 = 2 * mur_a * sum1 + log((rp1 * (1 - (di[k] / rp1)^2)) / re)
                q[k, i] = sum1
            else
                for n in 1:Ntrunc
                    sum1 -= (cos(theta_jk * n) * ((di[k] * di[i]) / rp1_sq)^n) / n
                    sum2 += (cos(theta_jk * n) * ((di[k] * di[i]) / rp1_sq)^n) / ((mur_a + 1) * n + y1 * kve_n_y1_ratio[n])
                end
                q[k, i] = sum1 + 2 * mur_a * sum2 + log(rp1 / sqrt(-(2 * cos(theta_jk) * di[k] * di[i]) + di_sq[k] + di_sq[i]))
            end
        end
    end

    zpp = ((jw * μ₀) * (q .+ mur_a / y1 * kve_n_y1_ratio[1])) / (2π)
    auxZ = zeros(Complex, size(zpp) .* ncx)
    for k in axes(zpp, 1)
        i1 = (k - 1) * ncx + 1
        i2 = i1 + ncx - 1
        for i in axes(zpp, 2)
            k1 = (i - 1) * ncx + 1
            k2 = k1 + ncx - 1
            auxZ[i1:i2, k1:k2] .= zpp[k, i]
        end
    end

    Zpipe = zeros(Complex, size(auxZ) .+ 1)
    for k in axes(auxZ, 1)
        for i in axes(auxZ, 2)
            Zpipe[k, i] = auxZ[k, i]
        end
    end

    return Zpipe
end


"""
(TYPEDSIGNATURES)

Calculates the Maxwell potential coefficient (`P`) per unit length matrix for a Pipe-Type (PT) cable as a function of the metallic sheath.

The cables inside the pipe are considered to be identical single-core (SC), each composed of a central conductor, internal insulation, and optionally, a shield and external insulation.

# Arguments

- `radii_single_core`: array of the radius of each layer of the SC cables \\[m\\]. This must be either a length 3 or 5 array.
    `radii_single_core[1]`: inner radius of the core, use `0` if solid conductor.
    `radii_single_core[2]`: outer radius of the core.
    `radii_single_core[3]`: outer radius of the first insulation.
    `radii_single_core[4]`, Optional: outer radius of the conductive sheath.
    `radii_single_core[5]`, Optional: outer radius of the second insulation.
- `radii_armor`: array of the armor radii \\[m\\]. Must have length 3.
    `radii_armor[1]`: inner radius of the armor.
    `radii_armor[2]`: outer radius of the armor.
    `radii_armor[3]`: outer radius of the external insulation surrounding the armor.
- `center_distances`: array (length `Nd`) of the distances from the center of the Pipe to the center of each SC cable \\[m\\].
- `theta`: array of the angles of each SC cable inside the Pipe \\[rad\\].
- `epsr_internal`: Relative permittivity of the internal insulation layer surrounding the SC cables \\[dimensionless\\].
- `loss_factor`, Optional: array of the loss factor (tangent) of the dielectric materials \\[dimensionless\\].

# Returns

- `P`: Maxwell potential coefficient matrix \\[V/c/m\\], shape `(Nd * ncx + 1, Nd * ncx + 1)`.
"""
function calc_pipe_internal_potential(
    radii_single_core::AbstractVector{<:Number},
    radii_armor::AbstractVector{<:Number},
    center_distances::AbstractVector{<:Number},
    theta::AbstractVector{<:Number},
    epsr_internal::Number,
    loss_factor::Number = 0,
)
    rc = radii_single_core
    ra = radii_armor
    di = center_distances
    epsr_b = epsr_internal
    if length(radii_armor) != 3
        throw(ArgumentError("It is expected the radii_armor to have length 3 or 5."))
    end
    # ncx: number of "active" conductors in each single-core cable
    if length(rc) == 3
        ncx = 1
    elseif length(rc) == 5
        ncx = 2
    else
        throw(ArgumentError("radii_single_core must have length 3 or 5"))
    end
    re = rc[end]
    rp1 = ra[1]
    Ntrunc = 26  # number of terms to truncate the infinite series
    rp1_sq = rp1^2
    di_sq = di .^ 2
    # potential coefficients
    qp = zeros(Complex, length(di), length(di))
    for k in eachindex(di)
        for i in eachindex(di)
            theta_jk = theta[i] - theta[k]
            sum1 = 0.0
            if k == i
                qp[k, i] = log((rp1 * (1 - (di[k] / rp1)^2)) / re)
            else
                for n in 1:Ntrunc
                    sum1 -= (cos(theta_jk * n) * ((di[k] * di[i]) / rp1_sq)^n) / n
                end
                qp[k, i] = sum1 + log(rp1 / sqrt(-(2 * cos(theta_jk) * di[k] * di[i]) + di_sq[k] + di_sq[i]))
            end
        end
    end

    ppp = qp ./ (2π * ε₀ * epsr_b) * (1 + 1im * loss_factor) / (1 + loss_factor^2)
    auxP = zeros(Complex, size(ppp) .* ncx)
    for k in axes(ppp, 1)
        i1 = (k - 1) * ncx + 1
        i2 = i1 + ncx - 1
        for i in axes(ppp, 2)
            k1 = (i - 1) * ncx + 1
            k2 = k1 + ncx - 1
            auxP[i1:i2, k1:k2] .= ppp[k, i]
        end
    end

    Ppipe = zeros(Complex, size(auxP) .+ 1)
    for k in axes(auxP, 1)
        for i in axes(auxP, 2)
            Ppipe[k, i] = auxP[k, i]
        end
    end

    return Ppipe
end


"""
(TYPEDSIGNATURES)

Calculates the impedance (`Z`) per unit length matrix for the metallic sheath (armor) of a Pipe-Type (PT) cable considering all internal single-core (SC) cables identical to each other.

# Arguments

- `complex_frequency`: Complex angular frequency `s = c + jω` \\[rad/s\\].
- `radii_single_core`: array of the radius of each layer of the tubular conductors \\[m\\]. This must be either a length 3 or 5 array.
    `radii_single_core[1]`: inner radius of the core, use `0` if solid conductor.
    `radii_single_core[2]`: outer radius of the core.
    `radii_single_core[3]`: outer radius of the first insulation.
    `radii_single_core[4]`, Optional: outer radius of the conductive sheath.
    `radii_single_core[5]`, Optional: outer radius of the second insulation.
- `radii_armor`: array of the armor radii \\[m\\]. Must have length 3.
    `radii_armor[1]`: inner radius of the armor.
    `radii_armor[2]`: outer radius of the armor.
    `radii_armor[3]`: outer radius of the external insulation surrounding the armor.
- `center_distances`: array (length `Nd`) of the distances from the center of the Pipe to the center of each single-core cables \\[m\\].
- `rho_armor`: Conductivity of the armor \\[dimensionless\\].
- `mu_r`: Relative permeability of the armor \\[dimensionless\\].

# Returns

- `Z`: Series impedance matrix \\[Ω/m\\], shape `(Nd * ncx + 1, Nd * ncx + 1)`.
"""
function calc_pipe_armor_impedance(
    complex_frequency::Number,
    radii_single_core::AbstractVector{<:Number},
    radii_armor::AbstractVector{<:Number},
    center_distances::AbstractVector{<:Number},
    rho_armor::Number,
    mur_armor::Number,
)
    rc = radii_single_core  # only used to determine ncx
    # ncx: number of "active" conductors in each single-core cable
    if length(rc) == 3
        ncx = 1
    elseif length(rc) == 5
        ncx = 2
    else
        throw(ArgumentError("radii_single_core must have length 3 or 5"))
    end

    ra = radii_armor
    di = center_distances
    jw = complex_frequency
    mur_a = mur_armor
    rho_a = rho_armor
    rai = ra[1]
    rae = ra[2]
    ma = sqrt((jw * μ₀ * mur_a) / rho_a)
    eta_i = ma * rai
    eta_e = ma * rae
    scale = exp(abs(real(eta_i)) - eta_e - abs(real(eta_e)) + eta_i)
    den = (besselix(1, eta_e) * besselkx(1, eta_i)) - (besselix(1, eta_i) * besselkx(1, eta_e)) * scale
    numi = (besselix(0, eta_i) * besselkx(1, eta_e)) * scale + (besselix(1, eta_e) * besselkx(0, eta_i))
    nume = (besselix(0, eta_e) * besselkx(1, eta_i)) + (besselix(1, eta_i) * besselkx(0, eta_e)) * scale
    zai = (ma * rho_a) / (2π * rai) * numi / den
    zae = (ma * rho_a) / (2π * rae) * nume / den
    se = exp(abs(real(eta_e)) - eta_i)
    zam = rho_a / (2π * rai * rae * den * se)

    rf = ra[3]
    zins = ((jw * μ₀) / (2π)) * log(rf / rae)
    zc1 = zai + zae - 2 * zam + zins
    zc2 = zae - zam + zins
    zc3 = zae + zins
    n = length(di) * ncx + 1
    Zp = zc1
    Zpipe = zeros(Complex, n, n)
    for i in 1:length(ra)
        i1 = 2 * (i - 1) + 1
        i2 = 2 * i
        for j in 1:length(ra)
            j1 = 2 * (j - 1) + 1
            j2 = 2 * j
            Zpipe[i1:i2, j1:j2] .= Zp
        end
    end

    Zpipe[end, :] .= zc2
    Zpipe[:, end] .= zc2
    Zpipe[end, end] = zc3
    return Zpipe
end


"""
(TYPEDSIGNATURES)

Calculates the Maxwell potential coefficient (`P`) per unit length matrix for the metallic sheath (armor) of a Pipe-Type (PT) cable considering all internal single-core (SC) cables identical to each other.

# Arguments

- `radii_single_core`: array of the radius of each layer of the tubular conductors \\[m\\]. This must be either a length 3 or 5 array.
    `radii_single_core[1]`: inner radius of the core, use `0` if solid conductor.
    `radii_single_core[2]`: outer radius of the core.
    `radii_single_core[3]`: outer radius of the first insulation.
    `radii_single_core[4]`, Optional: outer radius of the conductive sheath.
    `radii_single_core[5]`, Optional: outer radius of the second insulation.
- `radii_armor`: array of the armor radii \\[m\\]. Must have length 3.
    `radii_armor[1]`: inner radius of the armor.
    `radii_armor[2]`: outer radius of the armor.
    `radii_armor[3]`: outer radius of the external insulation surrounding the armor.
- `center_distances`: array (length `Nd`) of the distances from the center of the Pipe to the center of each single-core cables \\[m\\].
- `epsr_external`: Relative permittivity of the external dieletric material, surrounding the armor \\[dimensionless\\].
- `loss_factor`, Optional: array of the loss factor (tangent) of the dielectric materials \\[dimensionless\\].

# Returns

- `P`: Maxwell potential coefficient matrix \\[V/c/m\\], shape `(Nd * ncx + 1, Nd * ncx + 1)`.
"""
function calc_pipe_armor_potential(
    radii_single_core::AbstractVector{<:Number},
    radii_armor::AbstractVector{<:Number},
    center_distances::AbstractVector{<:Number},
    epsr_external::Number,
    loss_factor::Number = 0,
)
    rc = radii_single_core  # only used to determine ncx
    # ncx: number of "active" conductors in each single-core cable
    if length(rc) == 3
        ncx = 1
    elseif length(rc) == 5
        ncx = 2
    else
        throw(ArgumentError("radii_single_core must have length 3 or 5"))
    end

    ra = radii_armor
    di = center_distances
    epsr_a = epsr_external
    rae = ra[2]
    rf = ra[3]
    n = length(di) * ncx + 1
    pins = log(rf / rae) / (2π * ε₀ * epsr_a) * (1 + 1im * loss_factor) / (1 + loss_factor^2)
    Ppipe = fill(pins, n, n)
    return Ppipe
end


"""
(TYPEDSIGNATURES)

Computes the series impedance (`Z`) per unit length matrices of a Pipe-Type (PT) cable that has `N` identical single-core (SC) cables inside it.

The number of SC cables is determined by the length of the `center_distances` array.

The SC cables may be compromise of core conductor and, optionally, a metallic sheath. In the first case, the arrays `rho, mu_r, epsilon_r, loss_factor` must all have length 2. In the second case, they must all have length 3.

# Arguments

- `complex_frequencies`: Vector (length `Nf`) of the complex angular frequency `s = c + jω` \\[rad/s\\].
- `radii_single_core`: array of the radius of each layer of the SC cables \\[m\\]. This must be either a length 3 or 5 array.
    `radii_single_core[1]`: inner radius of the core, use `0` if solid conductor.
    `radii_single_core[2]`: outer radius of the core.
    `radii_single_core[3]`: outer radius of the first insulation.
    `radii_single_core[4]`, Optional: outer radius of the conductive sheath.
    `radii_single_core[5]`, Optional: outer radius of the second insulation.
- `radii_armor`: array of the armor radii \\[m\\]. Must have length 3.
    `radii_armor[1]`: inner radius of the armor.
    `radii_armor[2]`, Optional: outer radius of the armor.
    `radii_armor[3]`: outer radius of the external insulation surrounding the armor.
- `center_distances`: array (length `Nd`) of the distances from the center of the Pipe to the center of each single-core cables \\[m\\].
- `theta`: array of the angles of each SC cable inside the Pipe \\[rad\\].
- `rho`: array of the conductivity of the conductor materials \\[dimensionless\\]. This must be either a length 2 or 3 array.
    `rho[1]`: core conductor conductivity.
    `rho[2]`, Optional: conductive sheath conductivity.
- `mu_r`: array of the Relative permeability of the conductor materials \\[dimensionless\\]. This must be either a length 2 or 3 array.
    `mu_r[1]`: core conductor permeability of the SC cables.
    `mu_r[2]`, Optional: conductive sheath permeability of the SC cables.
    `mu_r[3]`: metallic armor permeability of the PT cable.

# Returns

- `Z`: Series impedance matrix \\[Ω/m\\], shape `(Nd * ncx + 1, Nd * ncx + 1, Nf)`.

Where ncx = 1 or 2, depending if the SC cables have a metallic sheath or not.
"""
function comp_pipe_impedance(
    complex_frequencies::AbstractVector{<:Number},
    radii_single_core::AbstractVector{<:Number},
    radii_armor::AbstractVector{<:Number},
    center_distances::AbstractVector{<:Number},
    theta::AbstractVector{<:Number},
    rho::AbstractVector{<:Number},
    mu_r::AbstractVector{<:Number},
)
    rc = radii_single_core
    # ncx: number of "active" conductors in each single-core cable
    if length(rc) == 3
        ncx = 1
    elseif length(rc) == 5
        ncx = 2
    else
        throw(ArgumentError("radii_single_core must have length 3 or 5"))
    end
    ra = radii_armor
    di = center_distances
    nd = length(di)
    nc = ncx * nd + 1
    nf = length(complex_frequencies)
    Z = Array{Complex}(undef, nc, nc, nf)
    for f = 1:nf
        jw = complex_frequencies[f]
        if length(rc) == 3 && all(x -> length(x) == 2, (rho, mu_r))
            Zin = calc_single_core_impedance(jw, rc, rho[1], mu_r[1])
            Zpipe_int = calc_pipe_internal_impedance(jw, rc, ra, di, theta,rho[2], mu_r[2])
            Zpipe_ext = calc_pipe_armor_impedance(jw, rc, ra, di, rho[2], mu_r[2])
        elseif length(rc) == 5 && all(x -> length(x) == 3, (rho, mu_r))
            Zin = calc_single_core_impedance(jw, rc, rho[1:2], mu_r[1:2])
            Zpipe_int = calc_pipe_internal_impedance(jw, rc, ra, di, theta, rho[3], mu_r[3])
            Zpipe_ext = calc_pipe_armor_impedance(jw, rc, ra, di, rho[3], mu_r[3])
        else
            msg = "`rho, mu_r, epsilon_r, loss_factor` must be all of equal length, either 2 when `length(radii_single_core) == 3` or 3 when `length(radii_single_core) == 5`."
            throw(ArgumentError(msg))
        end

        nn = nd * ncx + 1
        Zint = zeros(Complex, (nn, nn))
        for i = 1:nd
            i1 = 2 * i - 1
            i2 = 2 * i
            Zint[i1:i2, i1:i2] = Zin
        end
        Zf = @view Z[:, :, f]
        Zf .= Zint + Zpipe_int + Zpipe_ext

        # force explicit symmetry because of float arithmetic errors
        Zf .= (Zf + transpose(Zf)) / 2
    end
    return Z
end


"""
(TYPEDSIGNATURES)

Computes the Maxwell potential coefficient (`P`) per unit length matrix of a Pipe-Type (PT) cable that has `N` identical single-core (SC) cables inside it.

The number of SC cables is determined by the length of the `center_distances` array.

The SC cables may be compromise of core conductor and, optionally, a metallic sheath. In the first case, the arrays `rho, mu_r, epsilon_r, loss_factor` must all have length 2. In the second case, they must all have length 3.

# Arguments

- `radii_single_core`: array of the radius of each layer of the SC cables \\[m\\]. This must be either a length 3 or 5 array.
    `radii_single_core[1]`: inner radius of the core, use `0` if solid conductor.
    `radii_single_core[2]`: outer radius of the core.
    `radii_single_core[3]`: outer radius of the first insulation.
    `radii_single_core[4]`, Optional: outer radius of the conductive sheath.
    `radii_single_core[5]`, Optional: outer radius of the second insulation.
- `radii_armor`: array of the armor radii \\[m\\]. Must have length 3.
    `radii_armor[1]`: inner radius of the armor.
    `radii_armor[2]`, Optional: outer radius of the armor.
    `radii_armor[3]`: outer radius of the external insulation surrounding the armor.
- `center_distances`: array (length `Nd`) of the distances from the center of the Pipe to the center of each single-core cables \\[m\\].
- `theta`: array of the angles of each SC cable inside the Pipe \\[rad\\].
- `epsilon_r`: array of the Relative permittivity of the dieletric materials \\[dimensionless\\]. This must be either a length 2 or 3 array.
    `epsilon_r[1]`: permittivity of the first insulation of the SC cables.
    `epsilon_r[2]`, Optional: permittivity of the second insulation of the SC cables.
    `epsilon_r[3]`: permittivity of the outer insulation of the PT cable, surrounding the armor.
- `loss_factor`, Optional: array of the loss factor (tangent) of the dielectric materials \\[dimensionless\\]. This must be either a length 2 or 3 array.
    `loss_factor[1]`: loss factor of the first insulation of the SC cables.
    `loss_factor[2]`, Optional: loss factor of the second insulation of the SC cables.
    `loss_factor[3]`: loss_factor of the outer insulation of the PT cable, surrounding the armor (unused for now).

# Returns

- `P`: Maxwell potential coefficient matrix \\[V/c/m\\], shape `(Nd * ncx + 1, Nd * ncx + 1)`.

Where ncx = 1 or 2, depending if the SC cables have a metallic sheath or not.
"""
function comp_pipe_potential(
    radii_single_core::AbstractVector{<:Number},
    radii_armor::AbstractVector{<:Number},
    center_distances::AbstractVector{<:Number},
    theta::AbstractVector{<:Number},
    epsilon_r::AbstractVector{<:Number},
    loss_factor::AbstractVector{<:Number},
)
    rc = radii_single_core
    # ncx: number of "active" conductors in each single-core cable
    if length(rc) == 3
        ncx = 1
    elseif length(rc) == 5
        ncx = 2
    else
        throw(ArgumentError("radii_single_core must have length 3 or 5"))
    end
    ra = radii_armor
    di = center_distances
    nd = length(di)
    if length(rc) == 3 && all(x -> length(x) == 2, (epsilon_r, loss_factor))
        Pin = calc_single_core_potential(rc, epsilon_r[1], loss_factor[1])
        Ppipe_int = calc_pipe_internal_potential(rc, ra, di, theta, epsilon_r[1])
        Ppipe_ext = calc_pipe_armor_potential(rc, ra, di, epsilon_r[2])
    elseif length(rc) == 5 && all(x -> length(x) == 3, (epsilon_r, loss_factor))
        Pin = calc_single_core_potential(rc, epsilon_r[1:2], loss_factor[1:2])
        Ppipe_int = calc_pipe_internal_potential(rc, ra, di, theta, epsilon_r[2])
        Ppipe_ext = calc_pipe_armor_potential(rc, ra, di, epsilon_r[3])
    else
        msg = "`rho, mu_r, epsilon_r, loss_factor` must be all of equal length, either 2 when `length(radii_single_core) == 3` or 3 when `length(radii_single_core) == 5`."
        throw(ArgumentError(msg))
    end
    pint = zeros(Complex, size(Ppipe_ext))
    for i = 1:nd
        i1 = 2 * i - 1
        i2 = 2 * i
        pint[i1:i2, i1:i2] = Pin
    end
    P = pint + Ppipe_int + Ppipe_ext
    # force explicit symmetry because of float arithmetic errors
    P = (P + transpose(P)) / 2
    return P
end


"""
(TYPEDSIGNATURES)

Computes the shunt admittance (`Y`) per unit length matrices of a Pipe-Type (PT) cable that has `N` identical single-core (SC) cables inside it.

The number of SC cables is determined by the length of the `center_distances` array.

The SC cables may be compromise of core conductor and, optionally, a metallic sheath. In the first case, the arrays `rho, mu_r, epsilon_r, loss_factor` must all have length 2. In the second case, they must all have length 3.

# Arguments

- `complex_frequencies`: Vector (length `Nf`) of the complex angular frequency `s = c + jω` \\[rad/s\\].
- `radii_single_core`: array of the radius of each layer of the SC cables \\[m\\]. This must be either a length 3 or 5 array.
    `radii_single_core[1]`: inner radius of the core, use `0` if solid conductor.
    `radii_single_core[2]`: outer radius of the core.
    `radii_single_core[3]`: outer radius of the first insulation.
    `radii_single_core[4]`, Optional: outer radius of the conductive sheath.
    `radii_single_core[5]`, Optional: outer radius of the second insulation.
- `radii_armor`: array of the armor radii \\[m\\]. Must have length 3.
    `radii_armor[1]`: inner radius of the armor.
    `radii_armor[2]`, Optional: outer radius of the armor.
    `radii_armor[3]`: outer radius of the external insulation surrounding the armor.
- `center_distances`: array (length `Nd`) of the distances from the center of the Pipe to the center of each single-core cables \\[m\\].
- `theta`: array of the angles of each SC cable inside the Pipe \\[rad\\].
- `epsilon_r`: array of the Relative permittivity of the dieletric materials \\[dimensionless\\]. This must be either a length 2 or 3 array.
    `epsilon_r[1]`: permittivity of the first insulation of the SC cables.
    `epsilon_r[2]`, Optional: permittivity of the second insulation of the SC cables.
    `epsilon_r[3]`: permittivity of the outer insulation of the PT cable, surrounding the armor.
- `loss_factor`, Optional: array of the loss factor (tangent) of the dielectric materials \\[dimensionless\\]. This must be either a length 2 or 3 array.
    `loss_factor[1]`: loss factor of the first insulation of the SC cables.
    `loss_factor[2]`, Optional: loss factor of the second insulation of the SC cables.
    `loss_factor[3]`: loss_factor of the outer insulation of the PT cable, surrounding the armor (unused for now).

# Returns

- `Y`: Shunt admittance matrix per unit length of the cable \\[S/m\\], shape `(Nd * ncx + 1, Nd * ncx + 1, Nf)`.

Where ncx = 1 or 2, depending if the SC cables have a metallic sheath or not.
"""
function comp_pipe_admittance(
    complex_frequencies::AbstractVector{<:Number},
    radii_single_core::AbstractVector{<:Number},
    radii_armor::AbstractVector{<:Number},
    center_distances::AbstractVector{<:Number},
    theta::AbstractVector{<:Number},
    epsilon_r::AbstractVector{<:Number},
    loss_factor::AbstractVector{<:Number},
)
    rc = radii_single_core
    # ncx: number of "active" conductors in each single-core cable
    if length(rc) == 3
        ncx = 1
    elseif length(rc) == 5
        ncx = 2
    else
        throw(ArgumentError("radii_single_core must have length 3 or 5"))
    end
    ra = radii_armor
    di = center_distances
    nd = length(di)
    nc = ncx * nd + 1
    nf = length(complex_frequencies)
    P = comp_pipe_potential(rc, ra, di, theta, epsilon_r, loss_factor)
    inv_P = inv(P)
    # force explicit symmetry because of float arithmetic errors
    inv_P = (inv_P + transpose(inv_P )) / 2
    Y = Array{Complex}(undef, nc, nc, nf)
    for f = 1:nf
        Y[:, :, f] .= complex_frequencies[f] * inv_P
    end
    return Y
end


# ==============================================================================
# Other utility functions for transmission lines
# ==============================================================================


"""
(TYPEDSIGNATURES)

Computes the equivalent nodal admittance matrix of a transmission line from its per unit length series impedance (`Z`) and shunt admittance (`Y`) matrices.

# Arguments

- `Z`: Series impedance matrix per unit length of the cable \\[Ω/m\\].
- `Y`: Shunt admittance matrix per unit length of the cable \\[S/m\\].
- `line_length`: the total length \\[m\\] of the transmission line.

# Returns

- `YN`: equivalent nodal admittance matrix \\[S\\].
"""
function comp_equivalent_nodal_admittance(
    Z::AbstractMatrix{<:Number},
    Y::AbstractMatrix{<:Number},
    line_length::Number,
)
    n = size(Z, 1)
    sqrtzy = sqrt(Z * Y)
    Yc = inv(Z) * sqrtzy
    H = exp(-line_length * sqrtzy)
    H2 = H * H
    inv_IH2 = inv(I - H2)
    YL1 = Yc * (I + H2) * inv_IH2
    YL2 = -2.0 * Yc * H * inv_IH2
    yn = [[YL1  YL2]; [YL2  YL1]]
    # force explicit symmetry to cancel float rounding errors
    yn = (yn + transpose(yn)) / 2
    return yn
end


"""
(TYPEDSIGNATURES)

Eliminates eigenvector crossings using the Newton-Raphson method.

# Arguments

- `A`: Matrix for which the eigenproblem is to be solved.
- `V0`: Initial guess for the eigenvectors.
- `d0`: Initial guess for the eigenvalues.
- `tol`: Tolerance for convergence (default: 1e-9).
- `maxiter`: Maximum number of iterations (default: 10000).

# Returns

- `outd`: Computed eigenvalues.
- `outv`: Computed eigenvectors.
"""
function calc_eigen_NR(
    A::AbstractMatrix{<:Number},
    V0::AbstractMatrix{<:Number},
    d0::AbstractVector{<:Number},
    tol::Number = 1e-9,
    maxiter::Int = 10000,
    verbose::Bool = false,
)
    nc = size(A, 1)
    scale = norm(A)
    A /= scale
    outd = zeros(Complex, nc)
    outv = zeros(Complex, (nc, nc))
    jac = zeros(Complex, (nc + 1, nc + 1))
    vec1 = zeros(Complex, nc + 1)
    res1 = zeros(Complex, nc + 1)
    for i = 1:nc
        v = V0[:,i]
        d = d0[i]
        resi = maximum(abs.(A * v - v .* d))
        itr = 0
        while resi > tol && itr < maxiter
            itr += 1
            jac[1:nc, 1:nc] = A - d * I
            jac[1:nc, (nc+1)] = -v
            jac[(nc+1), 1:nc] = 2 * v
            # `dot(v, w)` conjugates v in julia...
            res1[:] = [A * v - v .* d ; dot(conj(v), v) - 1.0]
            for i = 1:(nc + 1)
                jac[i, i] += 1e-12  # to avoid 1/0
            end
            vec1[:] = qr(jac) \ res1
            resi = maximum(abs.(vec1))
            v = v - vec1[1:nc]
            d = d - vec1[nc+1]
        end
        if verbose && itr == maxiter
            println("Maximum iterations reached in calc_eigen_NR")
        end
        outd[i] = d * scale
        outv[:, i] = v ./ norm(v)
    end
    return outd, outv
end


"""
(TYPEDSIGNATURES)

Rotates the eigenvector matrix to minimize the imaginary part.  
This function complements the Newton-Raphson method.

# Arguments

- `S`: Eigenvector matrix.

# Returns

- Rotated eigenvector matrix with minimized imaginary part.
"""
function calc_eigen_rot(S::AbstractMatrix{<:Number})
    Nc = size(S)[1]
    scale1 = zeros(Complex, Nc)
    scale2 = copy(scale1)
    scale = copy(scale1)
    ang = copy(scale1)
    err1 = copy(scale1)
    err2 = copy(scale1)
    numerator = zeros(Nc)
    denominator = copy(numerator)
    SA = zeros(Complex, (Nc, Nc))
    SB = copy(SA)
    for col = 1:Nc
        for j = 1:Nc
            numerator[col] += imag(S[j, col]) * real(S[j, col])
            denominator[col] += real(S[j, col])^2 - imag(S[j, col])^2
        end
        numerator[col] *= -2
        ang[col] = 0.5 * atan(denominator[col], numerator[col])
        scale1[col] = cos(ang[col]) + 1.0im * sin(ang[col])
        scale2[col] = cos(ang[col] + pi / 2) + 1.0im * sin(ang[col] + pi / 2)
        aaa = bbb = ccc = ddd = eee = fff = 0.0
        
        for j = 1:Nc
            SA[j, col] *= scale1[col]
            SB[j, col] *= scale2[col]
            aaa += imag(SA[j, col])^2
            bbb += real(SA[j, col]) * imag(SA[j, col])
            ccc += real(SA[j, col])^2
            ddd += imag(SB[j, col])^2
            eee += real(SB[j, col]) * imag(SB[j, col])
            fff += real(SB[j, col])^2
        end
        err1[col] = aaa * cos(ang[col])^2 + bbb * sin(2 * ang[col]) + ccc * sin(ang[col])^2
        err2[col] = ddd * cos(ang[col])^2 + eee * sin(2 * ang[col]) + fff * sin(ang[col])^2
        if abs(err1[col]) < abs(err2[col])
            scale[col] = scale1[col]
        else
            scale[col] = scale2[col]
        end
        S[:, col] *= scale[col]
    end
    return S
end


"""
(TYPEDSIGNATURES)

Calculates the propagation velocity and attenuation of the propagation modes.

The attenuation `a` \\[Np/m\\] is the natural logarithm ratio of the wave at the start `v0` to the wave 1 m ahead `v1`:

```math
a = \\ln(v_0 / v_1)
```
```math
v_1 = v_0 / \\exp(a)
```

# Arguments

- `zc`: Complex array of shape (Nc, Nc, Nf). Longitudinal impedance matrix per unit length \\[Ω/m\\].
- `yc`: Complex array of shape (Nc, Nc, Nf). Shunt admittance matrix per unit length \\[S/m\\].
- `complex_frequencies`: Complex array of shape (Nf,). Vector of complex angular frequencies in the format `c + jω` \\[rad/s\\].

# Returns

- `propagation`: Complex array of shape (Nf, Nc). Propagation modes `a + jb`, where `a` is \\[Np/m\\] and `b` is \\[rad/s\\].
- `velocity`: Real array of shape (Nf, Nc). Propagation mode velocities in \\[m/μs\\].
- `attenuation`: Real array of shape (Nf, Nc). Propagation mode attenuation in \\[dB/m\\].

# Notes

The relation between Neper and decibels is: `1 dB = log(10) / 20 Np`.
"""
function calc_propagation_modes(
    zc::AbstractArray{<:Number, 3},
    yc::AbstractArray{<:Number, 3},
    complex_frequencies::AbstractVector{<:Number},
)
    nf = length(complex_frequencies)
    Nc = size(zc, 1)
    YZ = yc[:, :, 1] * zc[:, :, 1]
    evals, evecs = eigen(YZ)
    evecs = calc_eigen_rot(evecs)
    propagation = zeros(Complex, Nc, nf)
    velocity = zeros(Nc, nf)
    attenuation = similar(velocity)
    Np_to_db = log(10) / 20  # conversion factor from Neper to dB
    for k = 1:nf
        YZ = yc[:, :, k] * zc[:, :, k]
        evals, evecs = calc_eigen_NR(YZ, evecs, evals, 1e-6, 20000)
        propagation[:, k] = sqrt.(evals)
        velocity[:, k] = 1e-6 .* imag(complex_frequencies[k]) ./ imag.(propagation[:, k])
        attenuation[:, k] = real.(propagation[:, k]) * Np_to_db
    end
    return transpose.((propagation, velocity, attenuation))
end

Utils.@_autoexport

end