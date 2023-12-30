function [u,mesh,error_semi,h]=adaptative_refinement(tau,theReferenceElement,analytic)
doPlot=1;
error_semi=[];
h=[];

[mesh]=generateMeshSquare(doPlot);
u=compute_sol(mesh.X,mesh.T,theReferenceElement);
error=compute_H1_semi_error(u,mesh.X,mesh.T,theReferenceElement);

error_semi=[error_semi,error];
h=[h,1/size(mesh.X,1)^2];

if analytic
    while error^2>tau
        listElems=[];
        for i=1:size(mesh.T,1)
            Te=mesh.T(i,:);
            error_elem=compute_H1_semi_error_elemental(u,mesh.X(Te,:),Te,theReferenceElement);
            if error_elem^2>10^-5
                listElems=[listElems,i];
            end
        end 
        mesh=refineListElements(mesh,listElems,doPlot);
        u=compute_sol(mesh.X,mesh.T,theReferenceElement);
        error=compute_H1_semi_error(u,mesh.X,mesh.T,theReferenceElement);

        error_semi=[error_semi,error];
        h=[h,1/size(mesh.X,1)^2];
    end
else 
     while error^2>tau
        listElems=[];
        error_elem=computeZZelementalErrors(u,mesh.X,mesh.T,theReferenceElement);
        %maximum=max(error);
        for i=1:size(error_elem)
            if error_elem(i)^2>10^-5
                listElems=[listElems,i];
            end
        end
        mesh=refineListElements(mesh,listElems,doPlot);
        u=compute_sol(mesh.X,mesh.T,theReferenceElement);
        error=compute_H1_semi_error(u,mesh.X,mesh.T,theReferenceElement);

        error_semi=[error_semi,error];
        h=[h,1/size(mesh.X,1)^2];
        %disp(error);
     end
end
end 