"""
This file is deprecated and will be deleted later
Use include("src/Modified_SINDy.jl") instead
"""


using LinearAlgebra
using DifferentialEquations, ForwardDiff
using Statistics
using Random
using Flux
using Distributions
using Optimization , OptimizationOptimJL

abstract type Differentiation_Alg end

abstract type SINDy_Alg end

struct Default_SINDy <: SINDy_Alg 
end

mutable struct ADO <: SINDy_Alg
    x::Vector{Float64}
    q::Int16
    ω::Float64
    Nloop::Int16
    Ntrain::Int16
    opt::Flux.Optimise.AbstractOptimiser
    verbose::Bool
end

mutable struct OLE <: SINDy_Alg
    x::Vector{Float64}
    ADType::SciMLBase.AbstractADType
    opt::Optim.AbstractOptimizer
    lb::Vector{Float64}
    ub::Vector{Float64}
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

SINDy_Problem(u, du, dt, Lib, λs, iter, alg) = Modified_SINDy_Problem(u, du, dt, Lib, alg, λs, iter, 1, 0.000001, 0.000001, zeros(1,1))
SINDy_Problem(u, du, Lib, λs, iter) = Modified_SINDy_Problem(u, du, 1.0, Lib, Default_SINDy(), λs, iter, 1, 0.000001, 0.000001, zeros(1,1))
OLE(x, ADType, opt) = OLE(x, ADType, opt, [0], [0])


function solve_SINDy(prob::Modified_SINDy_Problem)
    if prob.alg isa Default_SINDy
        prob.Θ = prob.Lib(prob.u)
        return sparsify(prob)
    else
        return solve_SINDy(prob, prob.alg)
    end
    return 0
end

function CalDerivative(u,dt)
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

function sparsify(Θ,du,λs,iter; ρ=1, abstol=0.000001, reltol = 0.000001)
    # Θ: The library matrix with size n x p, where n is the data length, p is the number of nonlinear basis. 
    # du: Estimated or measured derivative of dynamics.
    # λ: Thresholding parameter.
    # iter: Number of regressions you would like to perform.
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

function sparsify(prob::Modified_SINDy_Problem)
    # standard sparsify function that accepts a sindy problem 
    return sparsify(prob.Θ,prob.du,prob.λs,prob.iter; ρ=prob.ρ, abstol=prob.abstol, reltol=prob.reltol)
end

SINDy_loss(du, Θ, Ξes, ρ) = sum(abs.(Θ * Ξes - du)) + ρ*count(x-> abs.(x)>0, Ξes)

function _is_converged(X, X_prev, abstol, reltol)::Bool
    Δ = norm(X .- X_prev)
    Δ < abstol && return true
    δ = Δ / norm(X)
    δ < reltol && return true
    return false
end

function ODE_RHS(u,x,Ξ)
    # u: The measuremt state matrix with size n x m, where n is the data length, m is the number of states. 
    # Ξ: Selection matrix with size p x m, where p is the number of candidate terms in the library, m is the number of state variables. 
    #--------------------------- 
    du=Lib(u,x)*Ξ
     
    return du
     
 end

 function RK_45_Forward(u,x,Ξ,dt)
    # u: The measuremt state matrix with size n x m, where n is the data length, m is the number of states. 
    # Ξ: Selection matrix with size p x m, where p is the number of candidate terms in the library, m is the number of state variables.
    # dt: Time step.
    # return ues: Estimate future state.
    #---------------------------
    Kf1=ODE_RHS(u,x,Ξ)*dt

    dum1=u+0.5*Kf1
    Kf2=ODE_RHS(dum1,x,Ξ)*dt

    dum2=u+0.5*Kf2
    Kf3=ODE_RHS(dum2,x,Ξ)*dt

    dum3=u+Kf3
    Kf4=ODE_RHS(dum3,x,Ξ)*dt
    #
    ues=u.+(1/6)*(Kf1.+2*Kf2.+2*Kf3.+Kf4)
    
    return ues
end

function RK_45_Backward(u,x,Ξ,dt)
    # u: The measuremt state matrix with size n x m, where n is the data length, m is the number of states. 
    # Ξ: Selection matrix with size p x m, where p is the number of candidate terms in the library, m is the number of state variables.
    # q: backward step.
    # dt: Time step.
    # return ues: Estimate past state.
    #---------------------------
    Kb1=ODE_RHS(u,x,Ξ)*dt

    dum1=u-0.5*Kb1
    Kb2=ODE_RHS(dum1,x,Ξ)*dt

    dum2=u-0.5*Kb2
    Kb3=ODE_RHS(dum2,x,Ξ)*dt

    dum3=u-Kb3
    Kb4=ODE_RHS(dum3,x,Ξ)*dt
    #
    ues=u.-(1/6)*(Kb1.+2*Kb2.+2*Kb3.+Kb4)

    return ues
end

function ADO_loss(Measurement,NoiseEs,x,Ξes,Ξactive,q,dt,ω)
    # Measurement: Measurment matrix with size of n x m, where n is the data length and m is the number state variables.
    # NoiseEs: Estimate of noise added to signal.
    # Ξes: Estimated selection vector.
    # x: Extended variable estimates
    # Ξactive: Activation matrix. If the elements in Ξactive is zero, then it means the correpsonding basis in the library is set to inactive.
    # q: Forward/backward simulation step.
    # dt: Time step.
    # ω: Decaying factor of loss value.
    # return LossVal: Value of loss function.
    #--------------------------- 
    # The selection vector that defines our ODE's right hand side should exclude the inactive nonlinear terms
    Ξdummy=Ξes.*Ξactive
    
    # Define matrix to store value
    FutureEstimate=[]
    PastEstimate=[]
    lossFuture=0
    lossPast=0
    lossEstimate=0
    
    # Simulate system based on intial conditions
    for i=1:q
        # If this is the first prediction
        if i==1
            # Estimate the future state, the initial condition's index ranges from t=1 to t=n-q, thus the estimated states'
            # index will range from t=2 to t=n-q+1
            FutureEstimate=RK_45_Forward(Measurement[1:end-q,:]-NoiseEs[1:end-q,:],x,Ξdummy,dt)
            
            # Estimate the past state, the initial condition's index ranges from t=q+1 to t=n, thus the estimated states'
            # index will range from t=q to t=n-1
            PastEstimate=RK_45_Backward(Measurement[q+1:end,:]-NoiseEs[q+1:end,:],x,Ξdummy,dt)
            
            # Calculate the future estimate loss
            lossFuture=norm(Measurement[1+i:end-q+i,:]-NoiseEs[1+i:end-q+i,:]-FutureEstimate).^2
            
            # Calculate the past estimate loss
            lossPast=norm(Measurement[q+1-i:end-i,:]-NoiseEs[q+1-i:end-i,:]-PastEstimate).^2
            
            # Add future and past loss together
            lossEstimate=(ω^(i-1))*(lossFuture+lossPast)
        else
            # Again estimate the future state one step forward based on previous simulated state value
            # The estimated states' index will range from t=1+i to t=n-q+i this time.
            FutureEstimate=RK_45_Forward(FutureEstimate,x,Ξdummy,dt) 

            # Similarly, estimate the past state one step backward based on previous simulated state
            # The estimated states' index will range from t=q+1-i to t=n-i this time.
            PastEstimate=RK_45_Backward(PastEstimate,x,Ξdummy,dt) 

            # Calculate the future estimate loss
            lossFuture=norm(Measurement[1+i:end-q+i,:]-NoiseEs[1+i:end-q+i,:]-FutureEstimate).^2
            
            # Calculate the past estimate loss
            lossPast=norm(Measurement[q+1-i:end-i,:]-NoiseEs[q+1-i:end-i,:]-PastEstimate).^2
            
            # Add future and past loss together
            lossEstimate=lossEstimate+(ω^(i-1))*(lossFuture+lossPast)
        end
    end
            
    # After calculated the simulation error, calculate the derivative error
    #Calculate the numerical derivative on the estimated true state to get the approximation of true derivative
    EstimatedDerivative_LHS=CalDerivative(Measurement-NoiseEs,dt)
    
    #Calculate the derivative estimate by ODE right hand side
    EstimatedDerivative_RHS=Lib(Measurement[3:end-2,:]-NoiseEs[3:end-2,:],x)*Ξdummy
    
    #Calculate the loss of two derivative estimates
    lossDerivative=norm(EstimatedDerivative_LHS-EstimatedDerivative_RHS).^2
    
    # Finally, add all the loss together
    lossVal=lossEstimate+lossDerivative
    
    return lossVal
end

function trainLoss(Measurement,NoiseEs,x,Ξes,Ξactive,q,dt,ω,Ntrain,θ,opt, verbose)
    # Measurement: Measurment matrix with size of n x m, where n is the data length and m is the number state variables.
    # NoiseEs: Estimate of noise added to signal.
    # Ξes: Estimated selection vector.
    # Ξactive: Activation matrix. If the elements in Ξactive is zero, then it means the correpsonding basis in the library is set to inactive.
    # q: Forward/backward simulation step.
    # dt: Time step.
    # ω: Decaying factor of loss value.
    # Ntrain: Iteration of the optimizer.
    # θ: Optimization variable
    # opt: Optimizer.
    # return: none
    #--------------------------- 
    for i=1:Ntrain
        # Get the gradient of loss function
        grads=gradient(()->ADO_loss(Measurement,NoiseEs,x,Ξes,Ξactive,q,dt,ω),θ)
        
        # Update the optimization variable
        for p in θ
            Flux.Optimise.update!(opt, p, grads[p])
        end
        
        # Print training process
        if verbose && mod(i,Int(floor(Ntrain/5)))==0
            println("\t $((i/(Ntrain/5))*20) percent finished...Current cost value is $(ADO_loss(Measurement,NoiseEs,x,Ξes,Ξactive,q,dt,ω))")
        end
    end
        
end

function ADO_SINDy(Measurement,x,q,dt,ω,λms,iter,Nloop,Ntrain,opt; ρ=1)
    # Measurement: Measurment matrix with size of n x m, where n is the data length and m is the number state variables.
    # q: Forward/backward simulation step.
    # dt: Time step.
    # ω: Decaying factor of loss value.
    # λms: Thresholding parameter for modified SINDy.
    # iter: Regression iteration of SINDy algorithm.
    # Nloop: Number of optimization loop.
    # Ntrain: Iteration of the optimizer.
    # opt: Optimizer.
    # c_or_g: Choose wether to use GPU, available input: cpu, gpu
    # return: NoiseEs: Estimated noise. Ξes: Estimated parameters of nonlinear basis
    #---------------------------              
    # First estimate the noise, here we do not use any pre-smoothing and set the initial estimate of noise as zeros
    NoiseEs=zeros(size(Measurement))
    
    # After estimating the noise, we can estimate the true state by subtracting estimated noise from measurement
    StateEs0=Measurement-NoiseEs
    
    # Get the number of states
    n_state=size(Measurement,2)
    
    # Based on the estimated state, we can estimate the derivative of the state
    StateDerivativeEs=CalDerivative(StateEs0,dt)

    # Once we get the estimated state and its corresponding derivative, we can estimate the value of Ξes
    Θms=Lib(StateEs0[3:end-2,:], x)
    Ξes, loss_SINDy = sparsify(Θms,StateDerivativeEs,λms,iter; ρ=ρ, abstol=0.000001, reltol = 0.000001)
    X_prev = Lib(StateEs0[4:end-3,:], x) * Ξes
    
    # When first started, set the activation matrix as ones
    Ξactive=ones(size(Ξes))


    # Now start training loop
    for j=1:Nloop
        # Define optimization variables we would like to train
        θ=Flux.params([NoiseEs,x,Ξes])
        
        println("Loop $(j)...")
        
        # Training statrt!
        trainLoss(Measurement,NoiseEs,x,Ξes,Ξactive,q,dt,ω,Ntrain,θ,opt,true)
        
        println("Loop $(j) finished!")

        # Get the new estimate of state, and transfer this to cpu
        # Note: Here we discarded first and last q points since those point are only 
        # constrained using forward or backward dynamics, and not the both.
        StateEs=Measurement[q+1:end-q-1,:]-NoiseEs[q+1:end-q-1,:]
        
        # Bsed on the updated state estimate, calculate the estimated derivative
        StateDerivativeEs=CalDerivative(StateEs,dt)
        
        # Get the new library matrix based on updated state estimate 
        Θms=Lib(StateEs[3:end-2,:], x)

        Ξes, _ = sparsify(Θms,StateDerivativeEs,λms,iter; ρ=ρ, abstol=0.000001, reltol = 0.000001)
        

        # Display the new parameter matrix
        println("Current value of the Ξ is shown below...\n")
        display(Ξes)
        println("Current value of x is shown below...\n")
        display(x)
    end
    
    return NoiseEs, Ξes, x
end

function solve_SINDy(prob::Modified_SINDy_Problem, ADO_params::ADO)
    # ADO_SINDy function fixed to work with Modified_SINDy_Problem type
    #---------------------------
    x = ADO_params.x

    # First estimate the noise, here we do not use any pre-smoothing and set the initial estimate of noise as zeros
    NoiseEs=zeros(size(prob.du))
    
    # After estimating the noise, we can estimate the true state by subtracting estimated noise from measurement
    StateEs0=prob.u-NoiseEs
    
    # Get the number of states
    n_state=size(prob.du,2)
    
    # Based on the estimated state, we can estimate the derivative of the state
    StateDerivativeEs=CalDerivative(StateEs0,prob.dt)

    # Once we get the estimated state and its corresponding derivative, we can estimate the value of Ξes
    prob.Θ=prob.Lib(StateEs0[3:end-2,:], x)
    Ξes, loss_SINDy = sparsify(prob)
    X_prev = prob.Lib(StateEs0[4:end-3,:], x) * Ξes
    
    # When first started, set the activation matrix as ones
    Ξactive=ones(size(Ξes))


    # Now start training loop
    for j=1:ADO_params.Nloop
        # Define optimization variables we would like to train
        θ=Flux.params([NoiseEs,x,Ξes])
        
        ADO_params.verbose && println("Loop $(j)...")
        
        # Training statrt!
        trainLoss(prob.du, NoiseEs, x, Ξes, Ξactive, ADO_params.q, prob.dt, ADO_params.ω, ADO_params.Ntrain, θ, ADO_params.opt, ADO_params.verbose)
        
        ADO_params.verbose && println("Loop $(j) finished!")

        # Get the new estimate of state, and transfer this to cpu
        # Note: Here we discarded first and last q points since those point are only 
        # constrained using forward or backward dynamics, and not the both.
        StateEs=prob.du[ADO_params.q+1:end-ADO_params.q-1,:]-NoiseEs[ADO_params.q+1:end-ADO_params.q-1,:]
        
        # Bsed on the updated state estimate, calculate the estimated derivative
        StateDerivativeEs=CalDerivative(StateEs,prob.dt)
        
        # Get the new library matrix based on updated state estimate 
        prob.Θ=prob.Lib(StateEs[3:end-2,:], x)

        Ξes, _ = sparsify(prob)
        

        # Display the new parameter matrix
        if ADO_params.verbose
            println("Current value of the Ξ is shown below...\n")
            display(Ξes)
            println("Current value of x is shown below...\n")
            display(x)
        end
    end
    ADO_params.x = x
    return NoiseEs, Ξes, x
end

function solve_SINDy(prob::Modified_SINDy_Problem, OLE_params::OLE)
    prob.Θ = prob.Lib(prob.u, OLE_params.x)

    function lor_loss(x, p)
        prob.Θ = prob.Lib(prob.u, x)
        Ξes, loss = sparsify(prob)
        return loss
    end
    lor_optf = OptimizationFunction(lor_loss, OLE_params.ADType)
    if OLE_params.lb == [0] && OLE_params.ub == [0]
        lor_optprob = OptimizationProblem(lor_optf, OLE_params.x)
    else
        lor_optprob = OptimizationProblem(lor_optf, OLE_params.x, lb = OLE_params.lb, ub = OLE_params.ub)
    end
    OLE_params.x = solve(lor_optprob, OLE_params.opt).u

    prob.Θ = prob.Lib(prob.u, OLE_params.x)
    return sparsify(prob)..., OLE_params.x
end