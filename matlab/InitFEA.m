function FEA = InitFEA(Node, Bars, Supp, Load, E)

NNode = size(Node, 1); NElem = size(Bars, 1);

D = Node(Bars(:,2),:) - Node(Bars(:,1),:);
L = sqrt(sum(D.^2, 2)); N = D./L;
LocalInternalMap = [-N N];

ElemNDof = 4*ones(NElem, 1); Nval = sum(ElemNDof.^2);
i=zeros(Nval,1); j=zeros(Nval,1); e=zeros(Nval,1); k=zeros(Nval,1);

idx = 0; ElemDof = zeros(NElem,4);
for el = 1:NElem; Ke = LocalK(LocalInternalMap(el,:), E, L(el));
    NDof = ElemNDof(el); Idx = idx+1:idx+NDof^2;
    eDof = reshape([2*Bars(el,:)-1; 2*Bars(el,:)], NDof, 1);
    ElemDof(el,:) = eDof'; IdxI = repmat(eDof, 1, NDof); IdxJ = IdxI'; 
    i(Idx) = IdxI(:); j(Idx) = IdxJ(:); k(Idx) = Ke(:);
    e(Idx) = el*ones(length(Idx),1); idx = idx + NDof^2;
end

NLoad = length(Load.node); F = zeros(2*NNode, 1);
for iL = 1:NLoad
    F(2*Load.node(iL)-1) = F(2*Load.node(iL)-1) + Load.val(iL,1);
    F(2*Load.node(iL)) = F(2*Load.node(iL)) + Load.val(iL, 2);
end

NSupp = length(Supp.node); FixedDofs = zeros(1, 2*NSupp); idx = 0;
for iS = 1:NSupp
    if Supp.type(iS,1)==1; idx=idx+1; FixedDofs(idx) = 2*Supp.node(iS)-1; end
    if Supp.type(iS,2)==1; idx=idx+1; FixedDofs(idx) = 2*Supp.node(iS); end
end
FreeDofs = setdiff(1:2*NNode, FixedDofs(1:idx));

% return FEA dictionary
FEA = struct('NNode',NNode,'Node',Node,'NElem',NElem,'Elem',Bars,...
    'ElemNDof',ElemNDof,'ElemDof',ElemDof,'Supp',Supp,'Load',Load,...
    'LocalIntMap',LocalInternalMap,'L',L,'E',E,'F',F,'FreeDofs',FreeDofs,...
    'i',i,'j',j,'k',k,'e',e);

%--------------------------------------------------------------------------
function K = LocalK(LocalIntMap, E, L)
K = (E/L)*(LocalIntMap'*LocalIntMap);
