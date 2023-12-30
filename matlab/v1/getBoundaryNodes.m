function [boundaryNodes]=getBoundaryNodes(mesh)

i = pdesde(mesh.m.e);
bEdgeNodes = mesh.m.e(1:2,i);
boundaryNodes = unique(bEdgeNodes(:));

