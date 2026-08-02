%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = CantileverDomain(Demand,Arg)
  BdBox = [0 1.6 0 1];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  Dist = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  LeftEdgeNodes = find(abs(Node(:,1)-BdBox(1))<eps);
  LeftUpperNode = find(abs(Node(:,1)-BdBox(1))<eps & ...
                       abs(Node(:,2)-BdBox(4))<eps);
  RightCenterNode = find(abs(Node(:,1)-BdBox(2))<0.05*BdBox(2) & ...
                       abs(Node(:,2))<1.05*BdBox(4)/2 & ...
                       abs(Node(:,2))>0.95*BdBox(4)/2);
  RigthBottomNode = find(abs(Node(:,1)-BdBox(2))<eps & ...
                         abs(Node(:,2)-BdBox(3))<eps);
  FixedNodes = [LeftEdgeNodes];
  Supp = zeros(length(FixedNodes),3);
  Supp(:,1)=FixedNodes; Supp(:,2)=1; Supp(:,3)=1;
  n = length(RightCenterNode);
  Load = [RightCenterNode,zeros(n,1),-1*ones(n,1)/n];
  x = {Supp,Load};
%----------------------------------------------------  SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%