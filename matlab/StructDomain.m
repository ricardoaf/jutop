function [NODE, ELEM, SUPP, LOAD]= StructDomain(Nx, Ny, Lx, Ly, Supp, Load)

% generate NODE using meshgrid
X = linspace(0., Lx, Nx+1); Y = linspace(0., Ly, Ny+1);
[X, Y] = meshgrid(X, Y); NODE = [X(:) Y(:)]; Nn = size(NODE, 1);

% generate 4-node structured ELEM
k = 0; ELEM = zeros(Nx*Ny, 4);
for j = 1:Ny
    for i = 1:Nx
        n1 = (i-1)*(Ny+1)+j; n2 = i*(Ny+1)+j;
        k = k+1; ELEM(k,:) = [n1 n2 n2+1 n1+1];
    end
end

% tolerance for coords search
tol = 0.01*min(Lx/Nx, Ly/Ny);

% generate SUPP
SNode=zeros(Nn,1); SType=zeros(Nn,2); NSupp=0;
for s = Supp'
    node = FindNodeSet(NODE, s(1:2), tol);
    for j=node', NSupp=NSupp+1; SNode(NSupp)=j; SType(NSupp,:)=s(3:4); end
end
SUPP = struct('node',SNode(1:NSupp),'type',SType(1:NSupp,:));

% generate LOAD
LNode=zeros(Nn,1); LVal=zeros(Nn,2); NLoad=0;
for L = Load'
    node = FindNodeSet(NODE, L(1:2), tol);
    for j=node'; NLoad=NLoad+1; LNode(NLoad)=j; LVal(NLoad,:)=L(3:4); end
end
LOAD = struct('node',LNode(1:NLoad),'val',LVal(1:NLoad,:));

%--------------------------------------------------------------------------
% Find node-set function
function node = FindNodeSet(Node, pos, tol)

x = pos(1); y = pos(2);
X = Node(:,1); Y = Node(:,2);

xMin = -inf; xMax = inf; yMin = -inf; yMax = inf;
if x>=0, xMin=x-tol; xMax=x+tol; end
if y>=0, yMin=y-tol; yMax=y+tol; end

node = find((X>xMin)&(X<xMax)&(Y>yMin)&(Y<yMax));
