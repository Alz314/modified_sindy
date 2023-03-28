using Revise
include("..\\src\\Modified_Sindy.jl"); 
using .ModifiedSINDy
using DifferentialEquations, ForwardDiff
using Plots
using LinearAlgebra

function Lorenz(du, u, p, t)
    #  p: Parameter vector.
    # du: Derivative of the system (left-hand side).
    #  u: State of the system.
    #  t: Time index.
    # return du: Derivative of system.
    #--------------------------- 
    du[1] = p[1]*(u[2]-u[1])
    du[2] = u[1]*(p[2]-u[3])-u[2]
    du[3] = u[1]*u[2]-p[3]*u[3] 
     
    return du
     
end

function add_noise(NoiseLevel, u, dt)
    # Define the random seed
    Random.seed!(0)

    # Calculate noise magnitude
    NoiseMag=NoiseLevel*std(u,dims=1)/100

    # Generate noise
    Noise=NoiseMag.*randn(size(u))

    # Now, add the noise to the clean data to generate the noisy measurment data
    un=u+Noise;

    dun=CalDerivative(un,dt)

    # Discard the first and last two points
    un=un[3:end-2,:];
    Noise=Noise[3:end-2,:];
    return un, dun, Noise
end
 

# Define time step
dt=0.01

# Define final time
T1=25

# Define time horizon
tspan1=(0.0,T1)

# Define initial condition
u0=[5.0;5.0;25.0]

# Define parameters for simulating the 
p0=[10.0;28.0;8/3]

# Define the ODE problem
prob = ODEProblem(Lorenz,u0,tspan1,p0)

# Solve the ODE problem and get the clean data 
u=transpose(Array(solve(prob,saveat=dt)));

# Calculate the derivative
du=zeros(size(u))

# Calculate ground truth for derivative
for i=1:size(u,1)
    du[i,:]=Lorenz(zeros(1,3), u[i,:], p0, 0)
end
du=du[3:end-2,:]

plot(u)