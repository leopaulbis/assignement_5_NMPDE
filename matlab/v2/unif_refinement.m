function [u,mesh,error_semi,h]=unif_refinement(tau,theReferenceElement)
doPlot=1;
error_semi=[];
h=[];

[mesh]=generateMeshSquare(doPlot);
u=compute_sol(mesh.X,mesh.T,theReferenceElement);
error=compute_H1_semi_error(u,mesh.X,mesh.T,theReferenceElement);

error_semi=[error_semi,error];
num=size(mesh.X,1);
h=[h,1/num^2];

while error^2>tau
    listElems=1:size(mesh.T,1);%numéro des éléments dans le maillage 
    mesh=refineListElements(mesh,listElems,doPlot);
    u=compute_sol(mesh.X,mesh.T,theReferenceElement);
    error=compute_H1_semi_error(u,mesh.X,mesh.T,theReferenceElement);
    
    error_semi=[error_semi,error];
    num=size(mesh.X,1);
    h=[h,1/num^2];
end

end 