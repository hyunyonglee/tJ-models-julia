using Base.Filesystem

# including files
include("tJ_models.jl")
include("run_model.jl")

# Creating directory function
function create_directory(dir_path)
    if !isdir(dir_path)
        mkpath(dir_path)
    end
end

# main function
let

    # Defining params
    Nx = 3
    Ny = 2
    N = Nx * Ny
    chi = 20
    nsweeps = 5
    sym = "d_wave"

    # Defining directory path
    dir_path = pwd() * "/data/$(sym)_$(Nx)x$(Ny)_chi_$(chi)"
    create_directory(dir_path * "/mps")
    create_directory(dir_path * "/logs")
    create_directory(dir_path * "/observables")

    # Defining angles
    θ₍min₎ = 0.0
    θ₍max₎ = 0.4
    θ₍step₎ = 0.1

    ϕ₍min₎ = -0.4
    ϕ₍max₎ = 0.4
    ϕ₍step₎ = 0.1

    # Defining params dictionary
    model_params = Dict{String,Any}(
        "Nx" => Nx,
        "Ny" => Ny,
        "t" => 1.0,
        "Dx" => 0.0,
        "Dy" => 0.0,
        "mu" => 0.0
    )

    # Defining dmrg params dictionary
    dmrg_params = Dict{String,Any}(
        "nsweeps" => nsweeps,
        "maxdim" => [50, 100, 150, 200, chi],
        "noise" => [1E-3, 1E-3, 1E-4, 1E-4, 1E-5, 1E-5, 1E-6, 1E-6, 1E-7, 1E-7, 0],
        "cutoff" => 1E-9
    )

    # Loop over angles
    for θ = θ₍min₎:θ₍step₎:θ₍max₎
        for ϕ = ϕ₍min₎:ϕ₍step₎:ϕ₍max₎

            # Define model params
            t = cos(θ * π) * cos(ϕ * π)
            Δ = sin(θ * π) * cos(ϕ * π)
            μ = cos(θ * π) * sin(ϕ * π)

            if sym == "s_wave"
                c = 1.0
            elseif sym == "d_wave"
                c = -1.0
            end

            # Define Hamiltonian
            H, sites = PBCS_H(Nx, Ny, t, Δ, c * Δ, μ)

            # Run DMRG
            psi, e, e_var, dmrg_observer = run_dmrg(H, dmrg_params, sites)

            # Calculate observables
            Nups = expect(psi, "Nup")
            Ndns = expect(psi, "Ndn")
            Cups = expect(psi, "Cup")
            Cdns = expect(psi, "Cdn")

            # Write observables to the file
            writing_observable(Nups, dir_path * "/observables/Nups_theta_$(θ)_phi_$(ϕ).txt")
            writing_observable(Ndns, dir_path * "/observables/Ndns_theta_$(θ)_phi_$(ϕ).txt")
            writing_observable(Cups, dir_path * "/observables/Cups_theta_$(θ)_phi_$(ϕ).txt")
            writing_observable(Cdns, dir_path * "/observables/Cdns_theta_$(θ)_phi_$(ϕ).txt")

            # Write mps
            writing_mps(psi, dir_path * "/mps/psi_theta_$(θ)_phi_$(ϕ).h5")

            # Write logs
            model_params["t"] = t
            model_params["Dx"] = Δ
            model_params["Dy"] = c * Δ
            model_params["mu"] = μ
            writing_logs(model_params, dmrg_params, dmrg_observer, dir_path * "/logs/log_theta_$(θ)_phi_$(ϕ).txt")

            # Write observables to the file
            f1 = open(dir_path * "/observables.txt", "a+")
            @printf(f1, "%.3f %.3f %.8f %.8f %.8f %.8f %.8f %.8f\n", θ, ϕ, e, e_var, sum(Nups) / N, sum(Ndns) / N, sum(Cups) / N, sum(Cdns) / N)
            close(f1)

            # println("θ = ", θ, ", ϕ = ", ϕ, ", t = ", t, ", Δ = ", Δ, ", μ = ", μ)

        end
    end
end


