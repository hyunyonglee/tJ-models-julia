using ITensors
using ITensors.HDF5
using Printf

dir_pbcs = "/home/hylee/tJ-models-julia/data/d_wave_5x4_chi_1000"
dir_tJ = "/home/hylee/tJ-models-julia/data/tJ_5x4_chi_1000"

# Loading tJ-model ground state
f = h5open( dir_tJ * "/mps/psi_t_1.0_J_0.5.h5", "r")
psi_tJ = read(f, "psi", MPS)
close(f)
sites_tJ = siteinds(psi_tJ)
        
# Defining angles
θ₍min₎ = 0.00
θ₍max₎ = 0.45
θ₍step₎ = 0.05

ϕ₍min₎ = -0.45
ϕ₍max₎ = 0.45
ϕ₍step₎ = 0.05

# Write observables to the file
f1 = open("overlaps.txt", "a+")

# Loop over angles
for θ = θ₍min₎:θ₍step₎:θ₍max₎
    for ϕ = ϕ₍min₎:ϕ₍step₎:ϕ₍max₎

        # Loading pbcs-model ground state
        f2 = h5open(dir_pbcs * "/mps/psi_theta_$(θ)_phi_$(ϕ).h5", "r")
        psi_pbcs = read(f2, "psi", MPS)
        close(f2)

        # Calculate the overlap between the two MPS
        sites_pbcs = siteinds(psi_pbcs)

        for j=1:length(sites_pbcs)
            replaceind!(psi_pbcs[j],sites_pbcs[j],sites_tJ[j])
        end

        # Calculate the overlap between the two MPS
        overlap = inner(psi_tJ, psi_pbcs)
        @show(overlap)
        @printf(f1, "%.3f %.3f %.8f\n", θ, ϕ, overlap)

    end
end

close(f1)
