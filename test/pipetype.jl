using Test
using LineCableModels
using CSV, DataFrames


@testset "3-phase Pipe-Type cable parameters" begin
    freq = exp10.(range(0, 9, 200))
    nf = length(freq)
    complex_frequencies = 2im * pi * freq
    rc = radii_single_core = [0.0, 9.6, 17.054, 18.054, 19.50] .* 1e-3
    ra = radii_armor = [48.0, 59.0, 65.0] .* 1e-3
    rho = [1.7241e-8, 2.2e-7, 2.86e-8]
    mu_r = [1.0, 1.0, 300.0]
    epsilon_r = [3.31, 2.3, 10.0]
    loss_factor = [0.0, 0.0, 0.0]  # FIXME tests for nonzero loss factor
    di = center_distances = fill(radii_single_core[end] / cos(deg2rad(30)), 3)
    theta = deg2rad.([0, 120, -120])

    @testset "Single-Core potential coefficient" begin
        df_p = CSV.read("data/single_core_P.csv", DataFrame)
        Pin = calc_single_core_potential(rc, epsilon_r[1:2], loss_factor[1:2])
        for i = 1:2
            for k = 1:2
                colp = "P_$(i)_$(k)"
                p_target = tryparse.(ComplexF64, replace.(df_p[!, colp], "(" => "", ")" => ""))
                @test Pin[k, i] ≈ p_target[1]
            end
        end
    end

    @testset "Single-Core impedance" begin
        df_z = CSV.read("data/single_core_Z.csv", DataFrame)
        Zin = zeros(ComplexF64, 2, 2, nf)
        for i = 1:nf
            jw = complex_frequencies[i]
            Zin[:, :, i] = calc_single_core_impedance(jw, rc, rho[1:2], mu_r[1:2])
        end
        for i = 1:2
            for k = 1:2
                colz = "Z_$(i)_$(k)"
                z_target = tryparse.(ComplexF64, replace.(df_z[!, colz], "(" => "", ")" => ""))
                @test Zin[k, i, :] ≈ z_target
            end
        end
    end

    @testset "Pipe internal potential coefficient" begin
        df_p = CSV.read("data/pipe_int_P_three_phase.csv", DataFrame)
        Ppipe_int = calc_pipe_internal_potential(rc, ra, di, theta, epsilon_r[2], loss_factor[2])
        for i = 1:2
            for k = 1:2
                colp = "P_$(i)_$(k)"
                p_target = tryparse.(ComplexF64, replace.(df_p[!, colp], "(" => "", ")" => ""))
                @test Ppipe_int[k, i] ≈ p_target[1]
            end
        end
    end

    @testset "Pipe internal impedance" begin
        df_z = CSV.read("data/pipe_int_Z_three_phase.csv", DataFrame)
        Zpipe_int = zeros(ComplexF64, 7, 7, nf)
        for i = 1:nf
            jw = complex_frequencies[i]
            Zpipe_int[:, :, i] = calc_pipe_internal_impedance(jw, rc, ra, di, theta, rho[3], mu_r[3])
        end
        for i = 1:7
            for k = 1:7
                colz = "Z_$(i)_$(k)"
                z_target = tryparse.(ComplexF64, replace.(df_z[!, colz], "(" => "", ")" => ""))
                @test Zpipe_int[k, i, :] ≈ z_target
            end
        end
    end

    @testset "Pipe external potential coefficient" begin
        df_p = CSV.read("data/pipe_ext_P_three_phase.csv", DataFrame)
        Ppipe_ext = calc_pipe_armor_potential(rc, ra, di, epsilon_r[3], loss_factor[3])
        for i = 1:7
            for k = 1:7
                colp = "P_$(i)_$(k)"
                p_target = tryparse.(ComplexF64, replace.(df_p[!, colp], "(" => "", ")" => ""))
                @test Ppipe_ext[k, i] ≈ p_target[1]
            end
        end
    end

    @testset "Pipe external impedance" begin
        df_z = CSV.read("data/pipe_ext_Z_three_phase.csv", DataFrame)
        Zpipe_ext = zeros(ComplexF64, 7, 7, nf)
        for i = 1:nf
            jw = complex_frequencies[i]
            Zpipe_ext[:, :, i] = calc_pipe_armor_impedance(jw, rc, ra, di, rho[3], mu_r[3])
        end
        for i = 1:7
            for k = 1:7
                colz = "Z_$(i)_$(k)"
                z_target = tryparse.(ComplexF64, replace.(df_z[!, colz], "(" => "", ")" => ""))
                @test Zpipe_ext[k, i, :] ≈ z_target
            end
        end
    end

    @testset "Full Y per unit length matrix" begin
        # as P is calculated internally for Y, it is tested implicitly
        Y = comp_pipe_admittance(
            complex_frequencies,
            radii_single_core,
            radii_armor,
            center_distances,
            theta,
            epsilon_r,
            loss_factor,
        )
        df = CSV.read("data/ZY_pipe_type_three_phase.csv", DataFrame)
        @test freq ≈ df[!, "Frequency_Hz"]
        for i = 1:7
            for k = 1:7
                coly = "Y_$(i)_$(k)"
                y_target = tryparse.(ComplexF64, replace.(df[!, coly], "(" => "", ")" => ""))
                @test Y[k, i, :] ≈ y_target atol = 1e-9
            end
        end
    end

    @testset "Full Z per unit length matrix" begin
        Z = comp_pipe_impedance(
            complex_frequencies,
            radii_single_core,
            radii_armor,
            center_distances,
            theta,
            rho,
            mu_r,
        )
        df = CSV.read("data/ZY_pipe_type_three_phase.csv", DataFrame)
        @test freq ≈ df[!, "Frequency_Hz"]
        for i = 1:7
            for k = 1:7
                colz = "Z_$(i)_$(k)"
                z_target = tryparse.(ComplexF64, replace.(df[!, colz], "(" => "", ")" => ""))
                @test Z[k, i, :] ≈ z_target
            end
        end
    end

    @testset "Equivalent Nodal Admittance matrix" begin
        df = CSV.read("data/ZY_pipe_type_three_phase.csv", DataFrame)
        freq = df[!, "Frequency_Hz"]
        nf = length(freq)
        Z = zeros(ComplexF64, 7, 7, nf)
        Y = similar(Z)
        for i = 1:7
            for k = 1:7
                colz = "Z_$(i)_$(k)"
                z_target = tryparse.(ComplexF64, replace.(df[!, colz], "(" => "", ")" => ""))
                Z[k, i, :] = z_target
                
                coly = "Y_$(i)_$(k)"
                y_target = tryparse.(ComplexF64, replace.(df[!, coly], "(" => "", ")" => ""))
                Y[k, i, :] = y_target
            end
        end
        YN = zeros(ComplexF64, 14, 14, nf)
        for f = 1:nf
            YN[:, :, f] = comp_equivalent_nodal_admittance(Z[:, :, f], Y[:, :, f], 2500.0)
        end
        df = CSV.read("data/YN_pipe_type_three_phase.csv", DataFrame)
        for i = 1:14
            for k = 1:14
                col = "YN_$(i)_$(k)"
                yn_target = tryparse.(ComplexF64, replace.(df[!, col], "(" => "", ")" => ""))
                @test YN[k, i, :] ≈ yn_target
            end
        end
    end

    @testset "Modal decomposition" begin
        df = CSV.read("data/ZY_pipe_type_three_phase.csv", DataFrame)
        freq = df[!, "Frequency_Hz"]
        nf = length(freq)
        Z = zeros(ComplexF64, 7, 7, nf)
        Y = similar(Z)
        for i = 1:7
            for k = 1:7
                colz = "Z_$(i)_$(k)"
                z_target = tryparse.(ComplexF64, replace.(df[!, colz], "(" => "", ")" => ""))
                Z[k, i, :] = z_target
                
                coly = "Y_$(i)_$(k)"
                y_target = tryparse.(ComplexF64, replace.(df[!, coly], "(" => "", ")" => ""))
                Y[k, i, :] = y_target
            end
        end
        propagation, velocity, attenuation = calc_propagation_modes(Z, Y, complex_frequencies)
        df = CSV.read("data/gamma_pipe_type_three_phase.csv", DataFrame)
        target = zeros(ComplexF64, nf, 7)
        for i = 1:7
            col = "gamma_$(i)"
            target[:, i] = tryparse.(ComplexF64, replace.(df[!, col], "(" => "", ")" => ""))
        end
        # rearrange the lines of A according to the smallest distance to B, column by column
        A = transpose(propagation)
        B = transpose(target)
        @test size(A) == size(B)
        n, m = size(A)
        A_reordered = similar(A)
        for col in 1:m
            A_col = A[:, col]
            B_col = B[:, col]
            # Find the permutation that aligns A_col with B_col
            permutation = sortperm(real.(A_col))[sortperm(sortperm(real.(B_col)))]
            A_reordered[:, col] = A_col[permutation]
        end
        @test A_reordered ≈ B
    end
end
