%----------------------------- PolyStress --------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyStress: A Matlab implementation%
% for topology optimization with local stress constraints using the       %
% augmented Lagrangian method", Structural and Multidisciplinary          %
% Optimization, DOI 10.1007/s00158-020-02664-7, 2020                      %
%-------------------------------------------------------------------------%
function [Node,Element,Supp,Load] = Mesh_Cantilever(Ne_ap)
L = 1;
% nn = 2*floor(round(sqrt(Ne_ap/4))/2); he = L/nn; 
nn = 2*floor(round(sqrt(5*Ne_ap/8))/2); he = L/nn; 
NElem = round(8/5*nn^2); 
[X1,Y1] = meshgrid(he/2:he:1.6*L-he/2,he/2:he:L-he/2); 
P = [X1(:) Y1(:)]; % Mesh seed
[Node,Element,Supp,Load,~] = PolyMesher(@CantileverDomain,NElem,0,P);
%-------------------------------------------------------------------------%