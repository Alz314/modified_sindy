module Modified_SINDy

    export solve_SINDy, sparsify, CalDerivative

    include("base.jl")
    include("OLE.jl")
    include("ADO_OLE.jl")

    using .SINDy_Base
    using .OLE_Module
    using .ADO_OLE_Module
end