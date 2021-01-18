using SparseArrays
#------------------------------------------------------------------------------
struct SparseData
    i::Vector{Int}
    j::Vector{Int}
    e::Vector{Int}
    k0::Vector{Float64}
    k::Vector{Float64}
end

#------------------------------------------------------------------------------
struct FEM
    NNode::Int
    Node::Array{Float64,2}
    NElem::Int
    Elem::Array{Int,2}
    ElemLen::Vector{Float64}
    ElemNDof::Vector{Int}
    CumSumElemNDof2::Vector{Int}
    ElemDof::Vector{Vector{Int}}
    Supp::Support
    Load::NodalLoad
    LocalIntMap::Array{Float64,2}
    YoungModulus::Float64
    Force::Vector{Float64}
    FreeDofs::Vector{Int}
    Sparse::SparseData
end

#------------------------------------------------------------------------------
function InitFEA(Node::Array{Float64,2}, Bars::Array{Int,2}, Supp::Support, Load::NodalLoad, E::Float64)::FEM

    NNode, NElem = size(Node, 1), size(Bars, 1)

    D = Node[Bars[:,2],:] - Node[Bars[:,1],:]
    L = vec(sqrt.(sum(D.^2, dims=2))); N = D./L
    LocalInternalMap = [-N N]

    ElemNDof = 4*ones(Int, NElem); CumSumElemNDof2 = cumsum(ElemNDof.^2)
    Ke = zeros(4, 4); Nval::Int = sum(ElemNDof.^2)
    i = zeros(Int, Nval); j = zeros(Int, Nval); e = zeros(Int, Nval); k0 = zeros(Nval)

    idx = 0; ElemDof = Vector{Int}[]; sizehint!(ElemDof, NElem)
    for el = 1:NElem
        LocalK!(Ke, LocalInternalMap[el,:], L[el], E)
        
        NDof = ElemNDof[el]; Idx = idx+1:idx+NDof^2
        EDof = permutedims([2Bars[el,:].-1 2Bars[el,:]])[:]
        push!(ElemDof, EDof)
        
        IdxI = repeat(EDof, outer=(1, NDof))
        IdxJ = permutedims(IdxI)
        i[Idx] = IdxI[:]; j[Idx] = IdxJ[:]; k0[Idx] = Ke[:]; e[Idx] .= el
        idx += NDof^2
    end
    k = zeros(Nval)

    F = zeros(2NNode)
    for iL = 1:length(Load.node)
        F[2Load.node[iL] - 1] += Load.val[iL, 1]
        F[2Load.node[iL]] += Load.val[iL, 2]
    end

    NSupp = length(Supp.node)
    FixedDofs, idx = zeros(Int, 2NSupp), 0
    for iS = 1:NSupp
        if Supp.type[iS,1] == 1; idx += 1; FixedDofs[idx] = 2Supp.node[iS]-1; end
        if Supp.type[iS,2] == 1; idx += 1; FixedDofs[idx] = 2Supp.node[iS];   end
    end
    FreeDofs = setdiff(1:2NNode, FixedDofs[1:idx])

    return FEM(NNode, Node, NElem, Bars, L, ElemNDof, CumSumElemNDof2, ElemDof,
    Supp, Load, LocalInternalMap, E, F, FreeDofs, SparseData(i, j, e, k0, k))
end

#------------------------------------------------------------------------------
function LocalK!(Ke::Array{Float64,2}, LocalIntMap::Vector{Float64}, L::Float64, E::Float64)

    NDof = length(LocalIntMap)
    for j = 1:NDof, i = 1:NDof
        Ke[i, j] = E / L * (LocalIntMap[i] * LocalIntMap[j])
    end
end

#------------------------------------------------------------------------------
function FEAnalysis!(U::Vector{Float64}, x::Vector{Float64}, FEA::FEM)

    NDof, Free = 2FEA.NNode, FEA.FreeDofs
    
    for idx in eachindex(FEA.Sparse.i)
        FEA.Sparse.k[idx] = FEA.Sparse.k0[idx] * x[FEA.Sparse.e[idx]]
    end
    K = sparse(FEA.Sparse.i, FEA.Sparse.j, FEA.Sparse.k, NDof, NDof)

    U[Free] = K[Free, Free] \ FEA.Force[Free]
end

#------------------------------------------------------------------------------
function Stress!(σ::Vector{Float64}, U::Vector{Float64}, FEA::FEM)

    for e = 1:FEA.NElem
        Δ = 0.
        for j = 1:FEA.ElemNDof[e]
            Δ += FEA.LocalIntMap[e, j] * U[FEA.ElemDof[e]][j]
        end
        σ[e] = FEA.YoungModulus / FEA.ElemLen[e] * Δ
    end
end
