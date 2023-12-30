function [u]=compute_sol(X,T,theReferenceElement)

%Boundary Dirichlet conditions
x = X(:,1); y = X(:,2); tol=1.e-10;
nodesCCD = find(abs(x)<tol|abs(x-1)<tol|abs(y)<tol|abs(y-1)<tol); %Nodes on the boundary
hold on, plot(x(nodesCCD),y(nodesCCD),'bo','MarkerSize',16); hold off
uCCD=DirichletValue(X(nodesCCD,:));

%System of equations (without BC)
[K,f]=computeSystemLaplace(X,T,theReferenceElement,@sourceTerm);

%Boundary conditions & SYSTEM SOLUTION
methodDBC = 1; % 1=System reduction or 2=Lagrange multipliers 
if methodDBC == 1 %System reduction
    unknowns= setdiff([1:size(X,1)],nodesCCD);
    f = f(unknowns)-K(unknowns,nodesCCD)*uCCD;
    K=K(unknowns,unknowns);
    %System solution
    sol=K\f;
    %Nodal values: system solution and Dirichlet values
    u = zeros(size(X,1),1);
    u(unknowns) = sol; u(nodesCCD) = uCCD;
else %LagrangeMultipliers
    nOfDirichletDOF = length(uCCD); nOfDOF = size(K,1);
    A = spalloc(nOfDirichletDOF,nOfDOF,nOfDirichletDOF);
    A(:,nodesCCD) = eye(nOfDirichletDOF);
    b = uCCD;
    Ktot = [K A'; A spalloc(nOfDirichletDOF,nOfDirichletDOF,0)];
    ftot = [f;b];
    sol = Ktot\ftot;
    u = sol(1:nOfDOF); lambda = sol(nOfDOF+1:end);
end
end