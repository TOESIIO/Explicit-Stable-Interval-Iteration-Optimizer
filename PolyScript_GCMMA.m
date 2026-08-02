%------------------------------- PolyStress ------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyStress: A Matlab implementation%
% for topology optimization with local stress constraints using the       %
% augmented Lagrangian method", Structural and Multidisciplinary          %
% Optimization, DOI 10.1007/s00158-020-02664-7, 2020                      %
%-------------------------------------------------------------------------%
% PolyScript_GCMMA.m: Code with GCMMA (MODIFIED FROM PolyScript.m)
%    Written in May 2026 by
%    Zeng Meng <mengz@hfut.edu.cn>
%    School of Civil Engineering, Hefei University of Technology.
%    Hefei 230009, PR China.
%-------------------------------------------------------------------------%
clear; clc; close all
nelx=400; nely=250; rmin=0.04*nelx;
restoredefaultpath; addpath(genpath('./')); %Use all folders and subfolders
set(0,'defaulttextinterpreter','latex')
%% ------------------------------------------------------------ CREATE Mesh
[Node,Element,Supp,Load] = Mesh_Cantilever(nelx*nely); 
NElem = size(Element,1); % Number of elements
%% ---------------------------------------------------- CREATE 'fem' STRUCT
E0 = 1; % E0 in MPa
G = E0/2.5; Et = E0; Ec = E0;  % 0<=(Et,Ec)<=3*G; %Material props. (linear)
fem = struct(...
  'NNode',size(Node,1),...      % Number of nodes
  'NElem',size(Element,1),...   % Number of elements
  'Node',Node,...               % [NNode x 2] array of nodes
  'Element',{Element},...       % [NElement x Var] cell array of elements
  'Supp',Supp,...               % Array of supports
  'Load',Load,...               % Array of loads
  'Passive',[],...              % Passive elements  
  'Thickness',1,...             % Element thickness
  'MatModel','Bilinear',...     % Material model ('Bilinear','Polynomial')
  'MatParam',[Et,Ec,G],...      % Material parameters for MatModel
  'SLim',25,...                 % Stress limit
  'TolR', 1e-8, ...             % Tolerance for norm of force residual
  'MaxIter', 15, ...            % Max NR iterations per load step
  'MEX', 'No');                 % Tag to use MEX functions in NLFEM routine         
%% ---------------------------------------------------- CREATE 'opt' STRUCT
R = rmin * 1.6/nelx;  q = 3;
p = 3.5; eta0 = 0.5;
m = @(y,B)MatIntFnc(y,'SIMP-H1',[p,B,eta0]);
P = PolyFilter(fem,R,q);
zIni = 0.5*ones(size(P,2),1);
opt = struct(... 
  'zMin',0.0,...              % Lower bound for design variables
  'zMax',1.0,...              % Upper bound for design variables
  'zIni',zIni,...             % Initial design variables
  'MatIntFnc',m,...           % Handle to material interpolation fnc.
  'contB',[5,1,1,10],...      % Threshold projection continuation params.  
  'P',P,...                   % Matrix that maps design to element vars.
  'Tol',0.002,...             % Convergence tolerance on design vars.
  'TolS',0.003,...            % Convergence tolerance on stress constraints
  'MaxIter',40,...           % Maximum number of AL steps
  'MMA_Iter',5,...            % Number of MMA iterations per AL step
  'lambda0',zeros(NElem,1),...% Initial Lagrange multiplier estimators
  'mu0',10,...                % Initial penalty factor for AL function
  'mu_max',10000,...          % Maximum penalty factor for AL function
  'alpha',1.1,...             % Penalty factor update parameter  
  'Move',0.2,...             % Allowable move step in MMA update scheme
  'Osc',0.2,...               % Osc parameter in MMA update scheme
  'AsymInit',0.2,...          % Initial asymptote in MMA update shecme
  'AsymInc',1.2,...           % Asymptote increment in MMA update scheme  
  'AsymDecr',0.7...           % Asymptote decrement in MMA update scheme
   );
%% ------------------------------------------------------- RUN 'PolyStress'
fem = preComputations(fem); % Run preComputations before running PolyStress
[z,V,fem] = PolyStress_GCMMA(fem,opt);
%-------------------------------------------------------------------------%