using LinearAlgebra, SparseArrays, Printf
function GenerateGS(Node, Elem, Lvl::Int, ColTol::Float64)

	# connectivity matrix
	Nn, Ne, Nne = size(Node,1), size(Elem,1), size(Elem,2)

	A1 = spzeros(Bool, Nn, Nn)
	for i=1:Ne; A1[Elem[i,:],Elem[i,:]] .= true; end
	A1 = A1 .* .~I(Nn); An = A1

	# Level 1 connectivity
	IdxJ, IdxI = findnz(An)
	Bars = [IdxI IdxJ]
	
	# Bar length and normalized direction
	D = Node[IdxI,:] - Node[IdxJ,:];
	D = D./sqrt.(sum(D.^2, dims=2))

	# Levels 2 and above
	for i = 2:Lvl
		Aold = An; An = map(x->x>0, An*A1); Gn = An - Aold
		IdxJ, IdxI = findnz(Gn .* .~I(Nn))
		if isempty(IdxJ); @printf("No new bars at Level %d\n", Lvl); Lvl = i-1; break; end

		# RmFlag = [] # considering no domain restrictions for now
		# deleteat!(IdxI, RmFlag)
		# deleteat!(IdxJ, RmFlag)

		newD = Node[IdxI,:] - Node[IdxJ,:]
		newD = newD./sqrt.(sum(newD.^2, dims=2))

		# Collinearity check
		p = 1; RmFlag = zeros(Int, length(IdxI)); Nb = size(Bars,1)
		for j = 1:Nn
			# Find IdxI[p:q] - NEW bars starting @ node 'j'
			for P=p:length(IdxI); p=P; if IdxI[p]>=j; break; end; end
			q = p; for Q=p:length(IdxI); q=Q; if IdxI[q]>j; break; end; end
			if IdxI[q]>j; q = q - 1; end

			if IdxI[p]==j
				# Find Bars[m:n] - OLD bars starting @ node 'j'
				m = 1; for M=1:Nb; m=M; if Bars[m,1]>=j; break; end; end
				n = m; for N=m:Nb; n=N; if Bars[n,1]>j; break; end; end
				if Bars[n,1]>j; n = n - 1; end

				if Bars[n,1]==j
					# Dot products of old vs. new bars, if collinear: mark
					C = vec(maximum(D[m:n,:] * newD[p:q,:]', dims=1))
					RmFlag[(p-1).+findall(C.>ColTol)] .= 1
				end
			end
		end

		# Remove collinear bars and make symmetric again.
		# Bars that have one angle marked as collinear but the other not, will be spared
		Ind = findall(RmFlag.==0)
		H = sparse(IdxI[Ind], IdxJ[Ind], true, Nn, Nn)
		IdxJ, IdxI = findnz(H+H')

		BarsIJ = [Bars; IdxI IdxJ]
		Bars = BarsIJ[sortperm(BarsIJ[:,1]),:]

		D = Node[Bars[:,1],:] - Node[Bars[:,2],:];
		D = D./sqrt.(sum(D.^2, dims=2))
	end

	# Only return bars {i,j} with i<j (no duplicated bars)
	A = sparse(Bars[:,1], Bars[:,2], true, Nn, Nn)
	IdxJ, IdxI = findnz(tril(A)); Bars = [IdxI IdxJ]
	return Bars
end
