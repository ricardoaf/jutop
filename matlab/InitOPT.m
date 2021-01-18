function Opt = InitOPT(XInitial, VolMax, Tol, MaxIter, ...
    XBounds, OCMove, OCEta)

if nargin<5 || isempty(XBounds), XBounds=mean(XInitial)*[1e-4 1e+4]; end
if nargin<6 || isempty(OCMove), OCMove=mean(XInitial)*1e+4; end
if nargin<7 || isempty(OCEta), OCEta=0.5; end

Opt = struct('XBounds',XBounds,'XInitial',XInitial,'VolMax',VolMax,...
    'Tol',Tol,'MaxIter',MaxIter,'OCMove',OCMove,'OCEta',OCEta);
