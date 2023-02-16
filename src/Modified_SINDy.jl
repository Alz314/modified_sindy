module ModifiedSINDy

    export Modified_SINDy_Problem, SINDy_Alg, Default_SINDy, BPSTLSQ, BandedBasisTerm # types
    export SINDy_Problem, term # constructors
    export solve_SINDy, sparsify, CalDerivative # functions
    export OLE, ADO

    include("base.jl")
    include("OLE.jl")
    include("ADO_OLE.jl")
    include("BPSTLSQ.jl")

    using .SINDy_Base
    using .OLE_Module
    using .ADO_OLE_Module
    using .BPSTLSQ_Module
end