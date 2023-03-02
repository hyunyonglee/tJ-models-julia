using ITensors
using Printf
using HDF5
using DelimitedFiles

# Running DMRG function
function run_dmrg(H, dmrg_params, sites)

    # Reading dmrg params
    nsweeps = dmrg_params["nsweeps"]
    maxdim = dmrg_params["maxdim"]
    noise = dmrg_params["noise"]
    cutoff = dmrg_params["cutoff"]

    # Run DMRG
    dmrg_observer = DMRGObserver(; energy_tol=1E-8)
    psi0 = randomMPS(sites, linkdims=20)
    e, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff, noise, observer=dmrg_observer, outputlevel=1)

    # Energy variance for sanity check
    e_var = inner(H, psi, H, psi) - e^2

    return psi, e, e_var, dmrg_observer

end


# Printing command line arguments
function print_params(model_params, dmrg_params)
    println("====================================")
    println("Input Model Params:")
    for (arg, val) in model_params
        println(" - ", arg, "  =>  ", val)
    end
    println("Input DMRG Params:")
    for (arg, val) in dmrg_params
        println(" - ", arg, "  =>  ", val)
    end
    println("====================================")
end


# Writing logs
function writing_logs(model_params, dmrg_params, dmrg_observer, log_file_name)

    # Write logs to the file
    f = open(log_file_name, "a+")
    @printf(f, "====================================\n")
    @printf(f, "Input Model Params:\n")
    for (arg, val) in model_params
        @printf(f, " - %s  =>  %s\n", arg, val)
    end
    @printf(f, "Input DMRG Params:\n")
    for (arg, val) in dmrg_params
        @printf(f, " - %s  =>  %s\n", arg, val)
    end
    @printf(f, "====================================\n")
    writedlm(f, transpose(dmrg_observer.energies))
    writedlm(f, transpose(dmrg_observer.truncerrs))
    close(f)

end


# Writing mps
function writing_mps(psi, mps_file_name)
    f = h5open(mps_file_name, "w")
    write(f, mps_file_name, psi)
    close(f)
end


# Writing observable
function writing_observable(observable, observable_file_name)
    writedlm(observable_file_name, transpose(observable))
end