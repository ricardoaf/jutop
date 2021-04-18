%% ----- INPUT DATA -----

Nx = 6; Ny = 6;	 	 		 % Number of cells
Lx = 6.; Ly = 6.; 	 		 % Domain size
Supp = [0. -1., 1 1];        % Support list: position and type
Load = [Lx Ly/2., 0. -100.]; % Load list: position and value
GSLvl = 6; 		 	 		 % Ground structure level
GSColTol = 0.999999;  		 % Ground structure colinear tolerance
VolMax = Lx*Ly/9000.; 		 % Volume constraint
Tol = 5e-9;                  % Optimization convergence tolerance
MaxIter = 4000;		 		 % Max number of optimization iterations
E = 7e+7;			 		 % Young's modulus
Filter = 0.01;               % Max area fraction threshold (for vieweing)

%% ----- PRE PROCESSING -----

% Create structured domain
disp('Creating structured domain ..'); tic;
[Node, Elem, Supp, Load] = StructDomain(Nx, Ny, Lx, Ly, Supp, Load); toc

% Generate Ground-Structure
disp('Generating Ground structure ..'); tic;
Bars = GenerateGS(Node, Elem, GSLvl, GSColTol); toc

% Init finite element analysis
disp('Initializing finite element analysis ..'); tic;
FEA = InitFEA(Node, Bars, Supp, Load, E); toc

% Calc initial areas
Area = VolMax/sum(FEA.L);
XInitial = Area*ones(FEA.NElem, 1);

% Init optimization
disp('Initializing optimization procedure ..'); tic;
Opt = InitOPT(XInitial, VolMax, Tol, MaxIter); toc

%% ----- OPTIMIZATION -----

% Run TrussTOP
disp('Running optimization procedure, TrussTOP ..'); tic;
[XHist, ObjFcnHist, FEA] = TrussTop(FEA, Opt); toc

%% ----- POST PROCESSING -----

% Show general results
TrussTOPView (XHist, ObjFcnHist, FEA, Filter);
