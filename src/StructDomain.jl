
struct Support; node::Vector{Int}; type::Array{Bool,2}; end
struct NodalLoad; node::Vector{Int}; val::Array{Float64,2}; end

function StructDomain(Nx::Int, Ny::Int, Lx::Float64, Ly::Float64, Supp, Load)

	# generate NODE using meshgrid
	X, Y = LinRange(0., Lx, Nx+1), LinRange(0., Ly, Ny+1)
	X, Y = MeshGrid2D(X, Y); NODE = [X[:] Y[:]]; Nn = size(NODE, 1)

	# generate 4-node structured ELEM
	k = 0; ELEM = zeros(Int, Nx*Ny, 4)
	for j=1:Ny, i=1:Nx
		n1, n2 = (i-1)*(Ny+1)+j, i*(Ny+1)+j
		k+=1; ELEM[k, :] = [n1 n2 n2+1 n1+1]
	end

	# tolerance for coords search
	tol = 0.01*min(Lx/Nx, Ly/Ny)

	# generate SUPP
	SNode=zeros(Int,Nn); SType=zeros(Bool,Nn,2); NSupp=0
	for s in Supp
		node = FindNodeSet(NODE, s[1:2], tol)
		for j in node; NSupp+=1; SNode[NSupp]=j; SType[NSupp,:]=s[3:4]; end
	end
	SUPP = Support(SNode[1:NSupp], SType[1:NSupp,:])

	# generate LOAD
	LNode=zeros(Int,Nn); LVal=zeros(Nn,2); NLoad=0
	for L in Load
		node = FindNodeSet(NODE, L[1:2], tol)
		for j in node; NLoad+=1; LNode[NLoad]=j; LVal[NLoad,:]=L[3:4]; end
	end
	LOAD = NodalLoad(LNode[1:NLoad], LVal[1:NLoad,:])

	return NODE, ELEM, SUPP, LOAD
end

#------------------------------------------------------------------------------
# Implementation of 2D mesh grid
function MeshGrid2D(x, y)
	nx, ny = length(x), length(y)
	return x' .* ones(ny), ones(nx)' .* y
end

#------------------------------------------------------------------------------
# Find node-set function
function FindNodeSet(Node, pos, tol::Float64)

	x, y = pos
	X, Y = Node[:,1], Node[:,2]

	xMin, xMax, yMin, yMax = -Inf, +Inf, -Inf, +Inf
	if x>=0; xMin=x-tol; xMax=x+tol; end
	if y>=0; yMin=y-tol; yMax=y+tol; end

	return findall((X.>xMin).&(X.<xMax).&(Y.>yMin).&(Y.<yMax))
end
