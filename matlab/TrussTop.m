function [xHist, fHist, FEA] = TrussTop (FEA, Opt)

Iter = 0; Tol = Opt.Tol; Change = Tol+1; x = Opt.XInitial;
xHist = zeros(FEA.NElem, Opt.MaxIter+1);
fHist = zeros(1, Opt.MaxIter+1);
while (Iter < Opt.MaxIter) && (Change>Tol)
    Iter = Iter + 1;
    % Compute cost functionals and analysis sensitivities
    [f, dfdx, FEA] = ObjectiveFcn(FEA, x);
    [g, dgdx, FEA] = ConstraintFcn(FEA, x, Opt.VolMax);
    % Store design vars. and obj. function history
    xHist(:,Iter) = x; fHist(Iter) = f;
    % Update design variable and analysis parameters
    [x, Change] = UpdateScheme(dfdx, g, dgdx, x, Opt);
end
% Calculate stress and internal force
[U, FEA] = FEAnalysis(FEA, x); FEA = FEStress(FEA, U);
% Store last design var. and obj. function value, then resize arrays
xHist(:,Iter+1) = x; [fHist(Iter+1), ~, FEA] = ObjectiveFcn(FEA,x);
xHist = xHist(:,1:Iter+1); fHist = fHist(1:Iter+1);

%--------------------------------------------------------------------------
function [f, dfdx, FEA] = ObjectiveFcn(FEA, x)

[U, FEA] = FEAnalysis(FEA, x);
f = dot(FEA.F, U);

temp = cumsum(-U(FEA.i) .* FEA.k .* U(FEA.j));
temp = temp(cumsum(FEA.ElemNDof.^2));
dfdx = [temp(1); temp(2:end)-temp(1:end-1)];

%--------------------------------------------------------------------------
function [g, dgdx, FEA] = ConstraintFcn(FEA, x, VolMax)
g = dot(FEA.L, x) - VolMax;
dgdx = FEA.L;

%--------------------------------------------------------------------------
function [xNew, Change] = UpdateScheme(dfdx, g, dgdx, x, Opt)
xMin = Opt.XBounds(1); xMax = Opt.XBounds(2);
move = Opt.OCMove; eta = Opt.OCEta;
Bm = -dfdx./dgdx; l1 = 0; l2 = 1.2*max(Bm);
while (l2-l1) > (1.0e-10*(1+l2))
    lmid = 0.5*(l1+l2); B = Bm./lmid;
    xCnd = xMin +(x-xMin).*(B.^eta);
    xNew = max(max(min(min(xCnd,x+move),xMax),x-move),xMin);
    if (g+dot(dgdx,xNew-x)>0), l1=lmid; else; l2=lmid; end
end
Change = max(abs(xNew-x)./(1+x));

%--------------------------------------------------------------------------
function [U, FEA] = FEAnalysis(FEA, x)
NDof = 2*FEA.NNode; FreeDofs = FEA.FreeDofs;
K = sparse(FEA.i, FEA.j, FEA.k.*x(FEA.e), NDof, NDof);
U = zeros(NDof, 1); U(FreeDofs) = K(FreeDofs, FreeDofs) \ FEA.F(FreeDofs);

%--------------------------------------------------------------------------
function FEA = FEStress(FEA, U)
FEA.Stress = zeros(FEA.NElem,1);
for el = 1:FEA.NElem
    Ue = U(FEA.ElemDof(el,:)); Delta = FEA.LocalIntMap(el,:)*Ue;
    FEA.Stress(el) = FEA.E/FEA.L(el) * Delta;
end