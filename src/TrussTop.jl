using Statistics, Plots

#------------------------------------------------------------------------------
struct OPT
	xMin::Float64
	xMax::Float64
	xIni::Vector{Float64}
	VolMax::Float64
	Tol::Float64
	MaxIter::Int
	OCMove::Float64
	OCEta::Float64
end

#------------------------------------------------------------------------------
function InitOPT(xIni::Vector{Float64}, VolMax::Float64, Tol::Float64, MaxIter::Int,
	xMin::Float64=1e-4*mean(xIni), xMax::Float64=1e+4*mean(xIni),
	OCMove::Float64=xMax, OCEta::Float64=0.5)::OPT

	OPT(xMin, xMax, xIni, VolMax, Tol, MaxIter, OCMove, OCEta)
end

#------------------------------------------------------------------------------
function TrussTop(FEA::FEM, Opt::OPT)::Tuple{Vector{Vector{Float64}}, Vector{Float64}, FEM}

	Iter, Tol, Change = 0, Opt.Tol, Opt.Tol + 1
	x = Opt.xIni; xTemp = similar(x)

	∂f∂x = zeros(FEA.NElem)
	∂g∂x = zeros(FEA.NElem)
	U = zeros(2FEA.NNode)

	xHist = Vector{Float64}[]; sizehint!(xHist, Opt.MaxIter + 1)
	fHist = Float64[]; sizehint!(fHist, Opt.MaxIter + 1)

	while (Iter < Opt.MaxIter) && (Change > Tol)
		Iter += 1

		# Compute cost functionals and analysis sensitivities
		f = ObjectiveFcn!(∂f∂x, U, x, FEA)
		g = ConstraintFcn!(∂g∂x, x, FEA, Opt)

		# Store design vars. and obj. function history
		push!(xHist, x)
		push!(fHist, f)

		# Update design variable and analysis parameters
		Change = UpdateScheme!(x, xTemp, ∂f∂x, g, ∂g∂x, Opt)
	end

	# Store last design var. and obj. function value
	f = ObjectiveFcn!(∂f∂x, U, x, FEA)
	push!(xHist, x); push!(fHist, f)

	# Return arguments
	xHist, fHist, FEA
end

#------------------------------------------------------------------------------
function ObjectiveFcn!(∂f∂x::Vector{Float64}, U::Vector{Float64}, x::Vector{Float64}, FEA::FEM)::Float64

	FEAnalysis!(U, x, FEA)

	@fastmath @inbounds @. FEA.Sparse.k = -U[FEA.Sparse.i] * FEA.Sparse.k0 * U[FEA.Sparse.j]
	@fastmath cumsum!(FEA.Sparse.k, FEA.Sparse.k)
	∂f∂x[1] = FEA.Sparse.k[FEA.CumSumElemNDof2[1]]
	@fastmath ∂f∂x[2:end] = diff(FEA.Sparse.k[FEA.CumSumElemNDof2])

	@fastmath f = FEA.Force ⋅ U
end

# function ObjectiveFcn!(∂f∂x::Vector{Float64}, U::Vector{Float64}, x::Vector{Float64}, FEA::FEM)::Float64

# 	FEAnalysis!(U, x, FEA)

# 	curr = 0.; prev = 0.; idx = 0
# 	for i in eachindex(FEA.Sparse.i)
# 		curr -= U[FEA.Sparse.i[i]] * FEA.Sparse.k0[i] * U[FEA.Sparse.j[i]]
# 		if (i % FEA.ElemNDof[idx+1]^2) == 0
# 			∂f∂x[idx+=1] = curr - prev
# 			prev = curr
# 		end
# 	end

# 	f = FEA.Force ⋅ U
# end

#------------------------------------------------------------------------------
function ConstraintFcn!(∂g∂x::Vector{Float64}, x::Vector{Float64}, FEA::FEM, Opt::OPT)

	@fastmath ∂g∂x .= FEA.ElemLen
	@fastmath g = FEA.ElemLen ⋅ x - Opt.VolMax
end

#------------------------------------------------------------------------------
function UpdateScheme!(x::Vector{Float64}, xTemp::Vector{Float64}, ∂f∂x::Vector{Float64},
	g::Float64, ∂g∂x::Vector{Float64}, Opt::OPT)::Float64
  
	xMin, xMax, move, η = Opt.xMin, Opt.xMax, Opt.OCMove, Opt.OCEta
  
	@fastmath @. ∂f∂x = min(∂f∂x, 0.)
	@fastmath @. xTemp = -∂f∂x / ∂g∂x
	@fastmath λ₁, λ₂ = 0., 1.2 * maximum(xTemp)
	
	@fastmath @. xTemp = x
	while (λ₂ - λ₁) > (1e-10 * (1 + λ₂))
	  λₘ = (λ₁ + λ₂) / 2.
	
	  @fastmath @. x = xMin + (xTemp - xMin) * ((-∂f∂x / ∂g∂x / λₘ) ^ η)
	  @fastmath @. x = min(x, xTemp + move, xMax)
	  @fastmath @. x = max(x, xTemp - move, xMin)
	
	  @fastmath if (g + ∂g∂x ⋅ (x - xTemp))>0 λ₁ = λₘ else λ₂ = λₘ end
	end
	
	@fastmath @. xTemp = abs(x - xTemp) / (1. + xTemp)
	@fastmath Δ = maximum(xTemp)
end

# function UpdateScheme!(x::Vector{Float64}, xTemp::Vector{Float64}, ∂f∂x::Vector{Float64},
# 	g::Float64, ∂g∂x::Vector{Float64}, Opt::OPT)::Float64

# 	xMin, xMax, move, η = Opt.xMin, Opt.xMax, Opt.OCMove, Opt.OCEta

# 	λ₁, λ₂, = 0., -Inf
# 	for i in eachindex(x)
# 		∂f∂x[i] = min(∂f∂x[i], 0.)
# 		λ₂ = max(λ₂, -∂f∂x[i] / ∂g∂x[i])
# 		xTemp[i] = x[i]
# 	end
# 	λ₂ *= 1.2

# 	while (λ₂ - λ₁) > (1e-10 * (1 + λ₂))
# 		c, λₘ = g, (λ₁ + λ₂) / 2.

# 		for i in eachindex(x)
# 			x[i] = xMin + (xTemp[i] - xMin) * ((-∂f∂x[i] / ∂g∂x[i] / λₘ) ^ η)
# 			x[i] = min(x[i], xTemp[i] + move, xMax)
# 			x[i] = max(x[i], xTemp[i] - move, xMin)

# 			c += ∂g∂x[i] * (x[i] - xTemp[i])
# 		end
# 		if c>0 λ₁ = λₘ else λ₂ = λₘ end
# 	end

# 	Δ = -Inf
# 	for i in eachindex(x)
# 		Δ = max(Δ, abs(x[i] - xTemp[i]) / (1 + xTemp[i]))
# 	end
# 	return Δ
# end

#------------------------------------------------------------------------------
function TrussTopView(xH::Vector{Vector{Float64}}, fH::Vector{Float64}, fem::FEM, Filter::Float64)

	# Apply visualization filter
	x = xH[end]; MaxArea = maximum(x);
	inactive = findall(x .< Filter * MaxArea)
	active = setdiff(1:length(x), inactive); 
	xEff = copy(x); xEff[inactive] .= 0

	# Summary results
	@printf("Number of bars: %d of %d\n", length(active), length(x));
	@printf("Objective function value: %g\n", fH[end]);
	@printf("Number of iterations: %d\n", length(fH) - 1);

	# Calc stress
	U = zeros(2fem.NNode); FEAnalysis!(U, x, fem)
	σ = zeros(fem.NElem); Stress!(σ, U, fem)

	# Final topology
	FinalTopology = plot()
	# println(active)
	for e in active
		X = vec(fem.Node[fem.Elem[e,:],1])
		Y = vec(fem.Node[fem.Elem[e,:],2])
		clr = "blue"; if σ[e]<0; clr="red"; end
		global FinalTopology = plot!(X,Y,color=clr,linewidth=xEff[e]/MaxArea*5,legend=false, title="Final topology",
		aspect_ratio=:equal, grid=false, axis=false, ticks=false)
	end
	ActiveNodes = unique(fem.Elem[active,:])
	X = vec(fem.Node[ActiveNodes,1])
	Y = vec(fem.Node[ActiveNodes,2])
	FinalTopology = scatter!(X,Y,color="white",markersize=5, legend=false)
	display(FinalTopology)
end
