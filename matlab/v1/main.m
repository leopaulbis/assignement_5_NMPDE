clear all;
%PREPROCES
%Reference element
degree = 1; typeOfElement=1; %1=TRI, 0=QUA
theReferenceElement = createReferenceElement(degree,typeOfElement);
nOfElementNodes = size(theReferenceElement.N,2);
%figure(1), drawReferenceElement(theReferenceElement);
%Mesh: regular mesh in a rectangular domain [0,1]x[0,1]
%[X,T] = CreateMesh(1,nOfElementNodes,[0,1,0,1],300,300);

% tau=10^-2;
% error=3;
% doPlot=1;
% [mesh]=generateMeshSquare(doPlot);
% u=compute_sol(mesh.X,mesh.T,theReferenceElement);
% error=compute_H1_semi_error(u,mesh.X,mesh.T,theReferenceElement)^2;
% 
% while error>tau
%     listElems=1:size(mesh.T,1);%numéro des éléments dans le maillage 
%     disp(size(mesh.T,1));
%     mesh=refineListElements(mesh,listElems,doPlot);
%     u=compute_sol(mesh.X,mesh.T,theReferenceElement);
%     error=compute_H1_semi_error(u,mesh.X,mesh.T,theReferenceElement)^2;%pour calculer l'erreure sur chaque élément, juste prendre que les x de l'élément, et la ligne correspondante de t .
%     disp(error);
% end


tau=0.5;
doPlot=1;
[mesh]=generateMeshSquare(doPlot);
u=compute_sol(mesh.X,mesh.T,theReferenceElement);
error=compute_H1_semi_error(u,mesh.X,mesh.T,theReferenceElement)^2;
disp(error);
while error>tau
    listElems=[];
    for i=1:size(mesh.T,1)
        Te=mesh.T(i,:);
        error_elem=compute_H1_semi_error_elemental(u,mesh.X(Te,:),Te,theReferenceElement)^2;
        %disp(error_elem^2);
        if error_elem^2>10^-6
            listElems=[listElems,i];
        end
    end 
    mesh=refineListElements(mesh,listElems,doPlot);
    u=compute_sol(mesh.X,mesh.T,theReferenceElement);
    error=compute_H1_semi_error(u,mesh.X,mesh.T,theReferenceElement)^2;
    disp(error);
end

% 
% 
% %POSTPROCESS
figure(3)
PlotNodalField(u,mesh.X,mesh.T), title('FEM solution')
% %[Xpg,Fgp] = gradientElementalField(u,X,T,theReferenceElement);
% % figure(4)
% % PlotMesh(T,X,1,'k-'); hold on, quiver(Xpg(:,1),Xpg(:,2),Fgp(:,1),Fgp(:,1),'LineWidth',2), hold off
% 
% %Comparison with analytical solution
% figure(5)
% [Xfine,Tfine] = CreateMesh(1,nOfElementNodes,[0,1,0,1],41,41);
% PlotNodalField(analytical(Xfine),Xfine,Tfine), title('Analytical solution')
% 
% % L2error=computeL2error(u,X,T,theReferenceElement);
% H1error=computeH1error(u,X,T,theReferenceElement);
% sem_error=compute_H1_semi_error(u,X,T,theReferenceElement);
% disp(sem_error^2);
% disp(H1error);


%PLOT OF THE L2 CONVERGENCE

% %Reference element
% degree = 1; typeOfElement=1; %1=TRI, 0=QUA
% theReferenceElement = createReferenceElement(degree,typeOfElement);
% nOfElementNodes = size(theReferenceElement.N,2);
% disp('H1');
% disp(H1error);
% 
% [node,error_lin_L2]=convergence_l2(5,11,nOfElementNodes,theReferenceElement);
% [node,error_lin_H1]=convergence_H1(5,11,nOfElementNodes,theReferenceElement);


% %Reference element
% degree = 2; typeOfElement=1; %1=TRI, 0=QUA
% theReferenceElement = createReferenceElement(degree,typeOfElement);
% nOfElementNodes = size(theReferenceElement.N,2);
% 
% [node,error_quad_L2]=convergence_l2(5,11,nOfElementNodes,theReferenceElement);
% [node,error_quad_H1]=convergence_H1(5,11,nOfElementNodes,theReferenceElement);
% 
% 
% 
% figure(5)
% 
% % Tracer la première série de données (erreur L2 pour les éléments linéaires)
% loglog(node, error_lin_L2, '-o', 'DisplayName', 'Linear Element L2 Error');
% 
% hold on % Pour superposer le prochain tracé

% % Tracer la deuxième série de données (erreur L2 pour les éléments quadratiques)
% loglog(node, error_quad_L2, '-s', 'DisplayName', 'Quadratic Element L2 Error');
% 
% legend('Location', 'Best') % Afficher la légende
% xlabel('log(element size)')
% ylabel('log(L2 Error)')
% 
% x = node; % Utiliser les mêmes valeurs x que pour les données
% y =x.^2; % L2 Error = c * h^2, où c est une constante
% loglog(x, y, '--', 'DisplayName', 'Reference Line (slope = 2)');
% 
% x = node; % Utiliser les mêmes valeurs x que pour les données
% y =x.^3; % L2 Error = c * h^2, où c est une constante
% loglog(x, y, '--', 'DisplayName', 'Reference Line (slope = 3)');
% legend();
% title('Convergence Plot')
% 
% 
% hold off 
% 
% figure(6)
% 
% loglog(node, error_quad_H1, '-o', 'DisplayName', 'quad Element H1 Error');
% hold on
% loglog(node, error_lin_H1, '-o', 'DisplayName', 'Linear Element H1 Error');
% x = node; % Utiliser les mêmes valeurs x que pour les données
% y =x.^2; % L2 Error = c * h^2, où c est une constante
% loglog(x, y, '--', 'DisplayName', 'Reference Line (slope = 2)');
% x = node; % Utiliser les mêmes valeurs x que pour les données
% y =x; % L2 Error = c * h^2, où c est une constante
% loglog(x, y, '--', 'DisplayName', 'Reference Line (slope = 1)');
% legend('Location', 'Best') % Afficher la légende
% xlabel('log(element size)')
% ylabel('log(L2 Error)')
% title("Convergence Plot")
% legend()



