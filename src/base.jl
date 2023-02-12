module SINDy_Base
    export Modified_SINDy_Problem, SINDy_Alg, Default_SINDy # types
    export SINDy_Problem # constructors
    export solve_SINDy, sparsify, CalDerivative # functions

    using LinearAlgebra

    
    # New Data types
    # ----------------------------------------------------------------------------------
    abstract type SINDy_Alg end

    struct Default_SINDy <: SINDy_Alg 
    end

    mutable struct Modified_SINDy_Problem
        u::Matrix{Float64}
        du::Matrix{Float64}
        dt::Float64
        Lib::Function
        alg::SINDy_Alg
        λs::Vector{Float64}
        iter::Int16
        ρ::Float64
        abstol::Float64
        reltol::Float64
        Θ::Matrix{Float64}
    end


    # Constructors for convenience
    # ----------------------------------------------------------------------------------
    SINDy_Problem(u, du, dt, Lib, λs, iter, alg) = Modified_SINDy_Problem(u, du, dt, Lib, alg, λs, iter, 1, 0.000001, 0.000001, zeros(1,1))
    SINDy_Problem(u, du, Lib, λs, iter) = Modified_SINDy_Problem(u, du, 1.0, Lib, Default_SINDy(), λs, iter, 1, 0.000001, 0.000001, zeros(1,1))


    # Functions 
    # ----------------------------------------------------------------------------------

    function solve_SINDy(prob::Modified_SINDy_Problem)
        # Takes in any Modified_SINDy_Problem and calls the correct function to solve the problem

        if prob.alg isa Default_SINDy
            prob.Θ = prob.Lib(prob.u) # for default SINDy, you don't need to input a Θ
            return sparsify(prob)
        else
            return solve_SINDy(prob, prob.alg) # call the corresponding solve function for the algorithm type
        end
        return 0
    end

    function sparsify(Θ,du,λs,iter; ρ=1, abstol=0.000001, reltol = 0.000001)
        # standard sparsify function 
        # Input: 
        #   Θ: The library matrix with size n x p, where n is the data length, p is the number of nonlinear basis. 
        #   du: Estimated or measured derivative of dynamics.
        #   λ: Thresholding parameter.
        #   iter: Number of regressions you would like to perform.
        # return ues: Estimate past state.
        #---------------------------
        # First get the number of states
        n_state=size(du,2)
        
        # Initialize comparison values
        min_loss = 9999999999
        best_iter = 0
        prev_smallinds = [1]

        # Next, perform regression, get an estimate of the selection matrix Ξes
        Ξes = Θ \ du
        X_prev = Θ * Ξes # our first SINDy prediction

        for i=1:iter
            # At each iteration, try all λ values to find the best one
            for λ in λs

                # Get the index of values whose absolute value is smaller than λ
                smallinds = (abs.(Ξes).<λ)

                # If the effect of λ is the same as the previous one, no need to do calculations again
                if smallinds == prev_smallinds
                    continue
                end

                # Make a temporary Ξes matrix to test out the effect of λ
                temp_Ξes = copy(Ξes)
                
                # Set the parameter value of library term whose absolute value is smaller than λ as zero
                temp_Ξes[smallinds].=0
                
                # Regress the dynamics to the remaining terms
                for ind=1:n_state
                    biginds = .!smallinds[:,ind]
                    temp_Ξes[biginds,ind] = Θ[:,biginds]\du[:,ind]
                end

                # Save the current small indices
                prev_smallinds = smallinds

                # calculate the loss and compare it to our best loss
                loss = SINDy_loss(du, Θ, temp_Ξes, ρ)
                if loss < min_loss
                    Ξes = copy(temp_Ξes)
                    min_loss = loss
                end
            end

            X = Θ * Ξes # make new SINDy prediction

            # If nothing, or very little, changed in one iteration, then we have converged
            if _is_converged(X, X_prev, abstol, reltol)
                break
            end
            
            X_prev = X
        end

        return Ξes, min_loss
    end

    SINDy_loss(du, Θ, Ξes, ρ) = sum(abs.(Θ * Ξes - du)) + ρ*count(x-> abs.(x)>0, Ξes)

    function _is_converged(X, X_prev, abstol, reltol)::Bool
        Δ = norm(X .- X_prev)
        Δ < abstol && return true
        δ = Δ / norm(X)
        δ < reltol && return true
        return false
    end

    function sparsify(prob::Modified_SINDy_Problem)
        # standard sparsify function that accepts a sindy problem 
        return sparsify(prob.Θ,prob.du,prob.λs,prob.iter; ρ=prob.ρ, abstol=prob.abstol, reltol=prob.reltol)
    end

    function CalDerivative(u,dt)
        # Simplest numerical differentiation algorithm (should probably move it ADO_OLE.jl)
        # Input:
        #  u: Measurement data we wish to approxiamte the derivative. 
        #  It should be of size n x m, where n is the number of measurement, m is the number of states.
        # dt: Time step
        # return du: The approximated derivative. 
        #--------------------------- 

        # First, get the information of the data length.
        len=size(u,1);
        n=size(u,2)
        
        # Define the coeficient for different orders of derivative
        p1=1/12;p2=-2/3;p3=0;p4=2/3;p5=-1/12;
        
        du=(p1*u[1:end-4,:]+p2*u[2:end-3,:]+p3*u[3:end-2,:]+p4*u[4:end-1,:]+p5*u[5:end,:])/dt;
            
        return du
    end
end