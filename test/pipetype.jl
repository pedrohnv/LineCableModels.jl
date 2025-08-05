using Test
using LineCableModels
using CSV, DataFrames

# Define a tolerance for floating-point comparisons

@testset "Pipe-Type formulas" begin
    # 3-phase Pipe-Type cable parameters
    freq = exp10.(range(0, 9, 200))
    nf = length(freq)
    complex_frequencies = 2im * pi * freq
    rc = radii_single_core = [0.0, 9.6, 17.054, 18.054, 19.50] .* 1e-3
    ra = radii_armor = [48.0, 59.0, 65.0] .* 1e-3
    rho = [1.7241e-8, 2.2e-7, 2.86e-8]
    mu_r = [1.0, 1.0, 300.0]
    epsilon_r = [3.31, 2.3, 10.0]
    loss_factor = [0.0, 0.0, 0.0]
    di = center_distances = fill(radii_single_core[end] / cos(deg2rad(30)), 3)
    
    @testset "Single-Core formulas" begin
        Zin = zeros(ComplexF64, 2, 2, nf)
        Yin = similar(Zin)
        for i = 1:nf
            jw = complex_frequencies[i]
            Zin[:, :, i], Yin[:, :, i] = calc_single_core_impedance_admittance(jw, rc, rho[1:2], mu_r[1:2], epsilon_r[1:2], loss_factor[1:2])
        end
        df_y = CSV.read("data/single_core_Y.csv", DataFrame)
        df_z = CSV.read("data/single_core_Z.csv", DataFrame)
        for i = 1:2
            for k = 1:2
                coly = "Y_$(i)_$(k)"
                y_target = tryparse.(ComplexF64, replace.(df_y[!, coly], "(" => "", ")" => ""))
                @test Yin[k, i, :] ≈ y_target
                colz = "Z_$(i)_$(k)"
                z_target = tryparse.(ComplexF64, replace.(df_z[!, colz], "(" => "", ")" => ""))
                @test Zin[k, i, :] ≈ z_target
            end
        end
    end

    @testset "Internal parameters" begin
        Zpipe_int = zeros(ComplexF64, 7, 7, nf)
        Ppipe_int = similar(Zpipe_int)
        for i = 1:nf
            jw = complex_frequencies[i]
            Zpipe_int[:, :, i], Ppipe_int[:, :, i] = calc_pipe_internal_impedance_potential(jw, rc, ra, di, rho[3], mu_r[3], epsilon_r[2])
        end
        df_p = CSV.read("data/pipe_int_P_three_phase.csv", DataFrame)
        df_z = CSV.read("data/pipe_int_Z_three_phase.csv", DataFrame)
        for i = 1:7
            for k = 1:7
                colp = "P_$(i)_$(k)"
                p_target = tryparse.(ComplexF64, replace.(df_p[!, colp], "(" => "", ")" => ""))
                @test Ppipe_int[k, i, :] ≈ p_target
                colz = "Z_$(i)_$(k)"
                z_target = tryparse.(ComplexF64, replace.(df_z[!, colz], "(" => "", ")" => ""))
                @test Zpipe_int[k, i, :] ≈ z_target
            end
        end
    end
    
    @testset "External parameters" begin
        Zpipe_ext = zeros(ComplexF64, 7, 7, nf)
        Ppipe_ext = similar(Zpipe_ext)
        for i = 1:nf
            jw = complex_frequencies[i]
            Zpipe_ext[:, :, i], Ppipe_ext[:, :, i] = calc_pipe_armor_impedance_potential(jw, rc, ra, di, rho[3], mu_r[3], epsilon_r[3])
        end
        df_p = CSV.read("data/pipe_ext_P_three_phase.csv", DataFrame)
        df_z = CSV.read("data/pipe_ext_Z_three_phase.csv", DataFrame)
        for i = 1:7
            for k = 1:7
                colp = "P_$(i)_$(k)"
                p_target = tryparse.(ComplexF64, replace.(df_p[!, colp], "(" => "", ")" => ""))
                @test Ppipe_ext[k, i, :] ≈ p_target
                colz = "Z_$(i)_$(k)"
                z_target = tryparse.(ComplexF64, replace.(df_z[!, colz], "(" => "", ")" => ""))
                @test Zpipe_ext[k, i, :] ≈ z_target
            end
        end
    end

    @testset "Full Z and Y per unit length matrices, equivalent Nodal Admittance and Modal Decomposition" begin
        Z, Y = comp_pipe_impedance_admittance(
            complex_frequencies,
            radii_single_core,
            radii_armor,
            center_distances,
            rho,
            mu_r,
            epsilon_r,
            loss_factor
        )
        df = CSV.read("data/ZY_pipe_type_three_phase.csv", DataFrame)
        @test freq ≈ df[!, "Frequency_Hz"]
        for i = 1:7
            for k = 1:7
                colz = "Z_$(i)_$(k)"
                z_target = tryparse.(ComplexF64, replace.(df[!, colz], "(" => "", ")" => ""))
                @test Z[k, i, :] ≈ z_target
                coly = "Y_$(i)_$(k)"
                y_target = tryparse.(ComplexF64, replace.(df[!, coly], "(" => "", ")" => ""))
                @test Y[k, i, :] ≈ y_target atol = 1e-9
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

end # PipeType Module testset
