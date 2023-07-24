using Base.Filesystem
using ITensors

ITensors.Strided.disable_threads()

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
    Nx = 4
    Ny = 4
    N = Nx * Ny
    chi = 100
    nsweeps = 50
    t = 1.0
    J = 0.5
    sym = "d_wave"

    # Defining directory path
    dir_path = pwd() * "/data/tJ_$(Nx)x$(Ny)_J_$(J)_chi_$(chi)"
    create_directory(dir_path * "/mps")
    create_directory(dir_path * "/logs")
    create_directory(dir_path * "/observables")

    # Defining N scan params
    Nf_min = 6
    Nf_max = 10
    Nf_step = 2

    # Defining params dictionary
    model_params = Dict{String,Any}(
        "Nx" => Nx,
        "Ny" => Ny,
        "t" => t,
        "J" => J,
        # "mu" => μ
    )

    # Defining dmrg params dictionary
    dmrg_params = Dict{String,Any}(
        "nsweeps" => nsweeps,
        "maxdim" => min.([100, 200, 400, 800, 2000, 3000, chi], chi),
        "noise" => [1E-6, 1E-7, 1E-8, 0.0],
        "cutoff" => 1E-8 # 1E-9
    )

    # Loop over J
    for Nf = Nf_min:Nf_step:Nf_max

        # Print params
        print_params(model_params, dmrg_params)

        # Define Hamiltonian
        H, sites = tJ_H(Nx, Ny, t, J)

        # Define initial state
        state = ["0" for n=1:N]
        for n = 1:Nf
            if n % 2 == 0
                state[n] = "Dn"
            else
                state[n] = "Up"
            end
        end
       
        @show(state)
        psi0 = randomMPS(sites, state, 100)
        # psi0 = productMPS(sites,state)

        # Run DMRG
        # psi, e, e_var, dmrg_observer = run_dmrg(H, dmrg_params, psi0)
        psi, e, dmrg_observer = run_dmrg(H, dmrg_params, psi0)
        @show( expect(psi, "N") )

        # Write mps
        writing_mps(psi, dir_path * "/mps/psi_t_$(t)_Nf_$(Nf).h5")

        # Write logs
        model_params["Nf"] = Nf
        writing_logs(model_params, dmrg_params, dmrg_observer, dir_path * "/logs/log_t_$(t)_Nf_$(Nf).txt")

        # Write observables to the file
        f = open(dir_path * "/observables.txt", "a+")
        @printf(f, "%.3f %.3f %.8f\n", J, e)
        close(f)

    end

end
