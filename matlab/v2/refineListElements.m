function [mesh]=refineListElements(mesh,list_ref,doPlot)


method1 = "regular";
method2 = "longest";
refinementMethod = method1;


g=generateSquareModel();

p=mesh.m.p;
e=mesh.m.e;
t=mesh.m.t;

mark_noRef = 1;
mark_ref   = mark_noRef+1;

t(4,:) = mark_noRef;
t(4,list_ref) = mark_ref;
[p,e,t] = refinemesh(g,p,e,t,mark_ref,refinementMethod);

if(nargin==3 && doPlot==1)
    pdemesh(p,e,t)
end

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


