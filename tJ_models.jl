using ITensors
using MKL

# Defining PBCS Hamiltonian
function PBCS_H(Nx, Ny, t, Dx, Dy, mu)

    # Define sites
    sites = siteinds("tJ", Nx * Ny; conserve_sz=true, conserve_nfparity=true)
    # sites = siteinds("tJ", Nx * Ny)

    # Defining MPO
    os = OpSum()
    for x = 1:Nx

        for y = 1:Ny

            s = (x - 1) * Ny + y
            su = s + 1
            sr = s + Ny

            if y == Ny
                su -= Ny
            end

            # print("x = ", x, ", y = ", y, ", s = ", s, ", su = ", su, ", sr = ", sr, "\n")

            # On-site terms
            os .+= -mu, "Ntot", s

            # Hopping
            os .+= -t, "Cdagup", s, "Cup", su
            os .+= -t, "Cdagdn", s, "Cdn", su
            os .+= -t, "Cdagup", su, "Cup", s
            os .+= -t, "Cdagdn", su, "Cdn", s

            # BCS
            os .+= Dy, "Cdagup", s, "Cdagdn", su
            os .+= -Dy, "Cdagdn", s, "Cdagup", su

            os .+= Dy, "Cdn", su, "Cup", s
            os .+= -Dy, "Cup", su, "Cdn", s

            if x < Nx
                # Hopping
                os .+= -t, "Cdagup", s, "Cup", sr
                os .+= -t, "Cdagdn", s, "Cdn", sr
                os .+= -t, "Cdagup", sr, "Cup", s
                os .+= -t, "Cdagdn", sr, "Cdn", s

                # BCS
                os .+= Dx, "Cdagup", s, "Cdagdn", sr
                os .+= -Dx, "Cdagdn", s, "Cdagup", sr

                os .+= Dx, "Cdn", sr, "Cup", s
                os .+= -Dx, "Cup", sr, "Cdn", s
            end

        end
    end

    return MPO(os, sites; splitblocks=true), sites
    # return MPO(os, sites), sites
end


# Defining t-J Hamiltonian
function tJ_H(Nx, Ny, t, J)

    # Define sites
    sites = siteinds("tJ", Nx * Ny; conserve_nf=true; conserve_sz=true, conserve_nfparity=true)
    
    # Defining MPO
    os = OpSum()
    for x = 1:Nx

        for y = 1:Ny

            s = (x - 1) * Ny + y
            su = s + 1
            sr = s + Ny

            if y == Ny
                su -= Ny
            end

            # On-site terms
            # os .+= -mu, "Ntot", s

            # Hopping: y-direction
            os .+= -t, "Cdagup", s, "Cup", su
            os .+= -t, "Cdagdn", s, "Cdn", su
            os .+= -t, "Cdagup", su, "Cup", s
            os .+= -t, "Cdagdn", su, "Cdn", s

            # Heisenberg: y-direction
            os .+= J, "Sz", s, "Sz", su
            os .+= J/2., "S+", s, "S-", su
            os .+= J/2., "S-", s, "S+", su
            os .+= -J/4, "Ntot", s, "Ntot", su
            
            if x < Nx
                # Hopping: x-direction
                os .+= -t, "Cdagup", s, "Cup", sr
                os .+= -t, "Cdagdn", s, "Cdn", sr
                os .+= -t, "Cdagup", sr, "Cup", s
                os .+= -t, "Cdagdn", sr, "Cdn", s

                # Heisenberg: x-direction
                os .+= J, "Sz", s, "Sz", sr
                os .+= J/2., "S+", s, "S-", sr
                os .+= J/2., "S-", s, "S+", sr
                os .+= -J/4, "Ntot", s, "Ntot", sr
            end
        end
    end

    return MPO(os, sites; splitblocks=true), sites
end