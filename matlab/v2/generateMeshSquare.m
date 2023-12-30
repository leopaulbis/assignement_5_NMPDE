function [mesh]=generateMeshSquare(doPlot)

g=generateSquareModel();

h = 0.75;

[p,e,t] = initmesh(g,"Hmax",h);
% if(nargin==1 && doPlot==1)
%     pdemesh(p,e,t)
% end

m.p = p;
m.e = e;
m.t = t;
mesh.m = m;
T = t';
T = T(:,1:3);
mesh.T = T;
mesh.X = p';

boundaryNodes=getBoundaryNodes(mesh);
mesh.boundaryNodes = boundaryNodes;

listElems = 1:size(mesh.T,1);
[mesh]=refineListElements(mesh,listElems,doPlot);