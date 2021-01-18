
include("StructDomain.jl")
include("GenerateGS.jl")
include("FEA.jl")
include("TrussTop.jl")

using BenchmarkTools

function TrussScript4()

    # ----- INPUT DATA -----

    Nx, Ny = 12, 12  	 	 		            # Number of cells
    Lx, Ly = 12., 12.  	 		                # Domain size
    Supp = [[0.,Ly/2., 1,1],[Lx,Ly/2., 1,1]]    # Support list: position and type
    Load = [[Lx/2.,Ly/2., 0.,-1.]]              # Load list: position and value
    GSLvl = 12 		 	 		                # Ground structure level
    GSColTol = 0.999999  		                # Ground structure colinear tolerance
    VolMax = Lx*Ly/9000. 		                # Volume constraint
    Tol = 1e-8			 		                # Optimization convergence tolerance
    MaxIter = 4000		 		                # Max number of optimization iterations
    E = 7e+7 			 		                # Young's modulus
    Filter = 0.01;                              # Max area fraction threshold (for vieweing)

    # ----- PRE PROCESSING -----

    # Create structured domain
    println("Creating structured domain ..")
    @time Node, Elem, Supp, Load = StructDomain(Nx, Ny, Lx, Ly, Supp, Load)

    # Generate Ground-Structure
    println("Generating Ground structure ..")
    @time Bars = GenerateGS(Node, Elem, GSLvl, GSColTol)

    # Init finite element analysis
    println("Initializing finite element analysis ..")
    @time FEA = InitFEA(Node, Bars, Supp, Load, E)

    # Calc initial areas
    Area = VolMax/sum(FEA.ElemLen)
    XInitial = Area*ones(FEA.NElem)

    # Init optimization
    println("Initializing optimization procedure ..")
    @time Opt = InitOPT(XInitial, VolMax, Tol, MaxIter)

    # ----- OPTIMIZATION -----

    # Run TrussTOP
    println("Running optimization procedure, TrussTOP ..")
    @time XHist, ObjFcnHist, FEA = TrussTop(FEA, Opt)

    # ----- POST PROCESSING -----

    # Show general results
    TrussTopView(XHist, ObjFcnHist, FEA, Filter)
end

TrussScript4()
