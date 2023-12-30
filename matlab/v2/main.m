clear all;
%PREPROCES
%Reference element
degree = 1; typeOfElement=1; %1=TRI, 0=QUA
theReferenceElement = createReferenceElement(degree,typeOfElement);
nOfElementNodes = size(theReferenceElement.N,2);
%figure(1), drawReferenceElement(theReferenceElement);
%Mesh: regular mesh in a rectangular domain [0,1]x[0,1]
%[X,T] = CreateMesh(1,nOfElementNodes,[0,1,0,1],300,300);

tau=10^-2;
[u_unif,mesh_unif,error_unif,h_unif]=unif_refinement(tau,theReferenceElement);
[u_adap,mesh_adap,error_adap,h_adap]=adaptative_refinement(tau,theReferenceElement,true);
[u_adap_zz,mesh_adap_zz,error_adap_zz,h_adap_zz]=adaptative_refinement(tau,theReferenceElement,false);
disp("number of nodes unif refinement");
disp(size(u_unif,1));
disp("number of elements unif refinement");
disp(size(mesh_unif.T,1));
% 
disp("number of nodes adaptative refinement");
disp(size(u_adap,1));
disp("number of elements adaptative refinement");
disp(size(mesh_adap.T,1));
%
disp("number of nodes adaptative ZZ-refinement");
disp(size(u_adap_zz,1));
disp("number of elements adaptative ZZ-refinement");
disp(size(mesh_adap_zz.T,1));
% % %POSTPROCESS
figure(2);
PlotNodalField(u_adap,mesh_adap.X,mesh_adap.T), title('FEM solution')
% 
figure(3);
 loglog(h_unif, error_unif, '-o', 'DisplayName', 'Uniform refinement');
% 
hold on % Pour superposer le prochain tracé

% Tracer la deuxième série de données (erreur L2 pour les éléments quadratiques)
loglog(h_adap, error_adap, '-s', 'DisplayName', 'Adaptative refinement');
loglog(h_adap_zz,error_adap_zz,'-s','DisplayName','Adaptative refinement ZZ estimator ')

legend('Location', 'Best') % Afficher la légende
xlabel('log(element size)')
ylabel('log(L2 Error)')
legend();
title('Convergence Plot')
hold off 

