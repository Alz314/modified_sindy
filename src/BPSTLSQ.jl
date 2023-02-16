module BPSTLSQ_Module
    export BPSTLSQ, BandedBasisTerm
    export SINDy_Problem, term 
    export solve_SINDy 

    # Must include these two lines (with the dots)
    import ..SINDy_Base: SINDy_Alg, Modified_SINDy_Problem, solve_SINDy, SINDy_Problem, SINDy_loss, _is_converged


    # New Data types
    # ----------------------------------------------------------------------------------
    mutable struct BPSTLSQ <: SINDy_Alg
        basis::AbstractArray{Function}
        upper_bounds::Vector{Float32}
        lower_bounds::Vector{Float32}
    end

    mutable struct BandedBasisTerm
        theta::Function
        upper_bound::Float32
        lower_bound::Float32
    end

    
    # Constructors for convenience
    # ----------------------------------------------------------------------------------
    function BPSTLSQ(library::AbstractArray{BandedBasisTerm})
        n = length(library)
        basis = Vector{Function}(undef, n)
        upper_bounds = zeros(n)
        lower_bounds = zeros(n)
        for i in 1:n
            basis[i] = library[i].theta
            upper_bounds[i] = library[i].upper_bound
            lower_bounds[i] = library[i].lower_bound
        end
        
        return BPSTLSQ(basis, upper_bounds, lower_bounds)
    end

    term(theta, ub, lb) = BandedBasisTerm(theta, ub, lb) # shorter name 

    function SINDy_Problem(u, du, dt, iter, alg::BPSTLSQ)
        function Lib(u)
            return reduce(hcat,[u |> entry for entry in alg.basis])
        end

        return Modified_SINDy_Problem(u, du, dt, Lib, alg, [0.0], iter, 1, 0.000001, 0.000001, Lib(u))
    end


    # Functions 
    # ----------------------------------------------------------------------------------
    """
        solve_SINDy(prob::Modified_SINDy_Problem, params::BPSTLSQ)

    Performs STLSQ algorithm with band passes. Each term in the basis has both an upper and lower bound it can be within.
    """
    function solve_SINDy(prob::Modified_SINDy_Problem, params::BPSTLSQ)
        # standard sparsify function 
        # Input: 
        #   Θ: The library matrix with size n x p, where n is the data length, p is the number of nonlinear basis. 
        #   du: Estimated or measured derivative of dynamics.
        #   λ: Thresholding parameter.
        #   iter: Number of regressions you would like to perform.
        # return ues: Estimate past state.
        #---------------------------
        # First get the number of states

        n_state=size(prob.du,2)
        
        # Initialize comparison values
        min_loss = 9999999999

        # Next, perform regression, get an estimate of the selection matrix Ξes
        Ξes = prob.Θ \ prob.du
        X_prev = prob.Θ * Ξes # our first SINDy prediction

        for i=1:prob.iter
            # Get the index of values that are between bounds
            inactive_inds = (abs.(Ξes) .> params.lower_bounds) .|| (abs.(Ξes) .< params.upper_bounds)

            # Make a temporary Ξes matrix to test out the effect of λ
            temp_Ξes = copy(Ξes)
            
            # Set the parameter value of library term whose absolute value is smaller than λ as zero
            temp_Ξes[inactive_inds].=0
            
            # Regress the dynamics to the remaining terms
            for ind=1:n_state
                active_inds = .!inactive_inds[:,ind]
                temp_Ξes[active_inds,ind] = prob.Θ[:,active_inds]\prob.du[:,ind]
            end

            # calculate the loss and compare it to our best loss
            loss = SINDy_loss(prob.du, prob.Θ, temp_Ξes, prob.ρ)
            if loss < min_loss
                Ξes = copy(temp_Ξes)
                min_loss = loss
            end

            X = prob.Θ * Ξes # make new SINDy prediction

            # If nothing, or very little, changed in one iteration, then we have converged
            if _is_converged(X, X_prev, prob.abstol, prob.reltol)
                break
            end
            X_prev = X
        end

        

        return Ξes, min_loss
    end
end