"""
Follow this template for any new algorithms. 
Make sure that you edit Modified_SINDy.jl to have:
    include("alg.jl")
    using .NewModuleName
Also, if you export any additional functions here, you need to export them in Modified_SINDy.jl too
"""

module NewModuleName
    export alg # need to export new data type
    export solve_SINDy # also export any other functions you would like to be able to call elsewhere

    # Must include these two lines (with the dots)
    include("base.jl")
    using ..SINDy_Base


    # New Data types
    # ----------------------------------------------------------------------------------
    mutable struct alg <: SINDy_Alg
        # Add any parameters for new algorithm
    end


    
    # Constructors for convenience
    # ----------------------------------------------------------------------------------



    # Functions 
    # ----------------------------------------------------------------------------------

    # For new functions, probably best practice to write documentation in this form above the function
    # It's a lot of extra work and not completely necessary, but it would be nice
    """
        solve_SINDy(prob::Modified_SINDy_Problem, params::alg)

    Explanation of function
    
    # Arguments
        * `prob`: Modified_SINDy_Problem type
        * etc.

    # Notes
        * Optional notes regarding the function 

    # Examples
    ```julia
    julia> # optionally write how this function should be called here
    output
    ```
    """
    function solve_SINDy(prob::Modified_SINDy_Problem, params::alg)
        # Implement new algorithm here
    end
end