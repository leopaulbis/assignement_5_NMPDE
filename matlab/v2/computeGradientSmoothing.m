function [g] = computeGradientSmoothing(u,X,T,theReferenceElement) 
% [dudx,dudy] = computeGradientSmoothing(u,X,T,theReferenceElement))
% Computes the nodal field corresponding to the the least-squares fitting
% of the gradient of u^h (given by its nodal values u)
% u^h is assumed to be scalar. Use this routine for each component in the
% case of vector fields.
% 
% u:            vector of nodal values of the FE solution
% X,T:          FE mesh: nodal coordinates and connectivities
% theReferenceElement: reference element

if(nargin<=3)
    degree = 1; typeOfElement=1; %1=TRI, 0=QUA
    theReferenceElement = createReferenceElement(degree,typeOfElement);
end

IPweights = theReferenceElement.IPweights;
IPcoord = theReferenceElement.IPcoord;
N=theReferenceElement.N;
Nxi=theReferenceElement.Nxi;
Neta=theReferenceElement.Neta;
 
[nOfElements,nOfElementNodes] = size(T); 
%numero de nodos total  
nOfNodes= size(X,1); 
 
%Memory allocation 
%M = zeros(nOfNodes,nOfNodes); 
M=spalloc(nOfNodes,nOfNodes,9*nOfNodes);
fx = zeros(nOfNodes,1); 
fy = zeros(nOfNodes,1); 


%Least-Squares fitting: computation of mass matrix and right-hand side vector
%Loop in elements
for ielem = 1:nOfElements 
 Te = T(ielem,:); 
 Xe = X(Te,:); 
 ue = u(Te);
 [Me,fxe,fye] = smoothingElementalComputations(ue,Xe,IPcoord,IPweights,N,Nxi,Neta); 
 M(Te,Te)=M(Te,Te)+Me; 
 fx(Te)=fx(Te)+fxe; 
 fy(Te)=fy(Te)+fye; 
end
%Solution of the Least-Squares fitting system
% dudx = M\fx;
% dudy = M\fy;
g = M\[fx fy];
% dudx = g(:,1);
% dudy = g(:,2);

%_________________________________________________________________ 
%
% Elemental computations
function [Me,fxe,fye] = smoothingElementalComputations(ue,Xe,IPcoordinates,IPweights,N,Nxi,Neta)

nOfElementNodes = length(ue);
Me = zeros(nOfElementNodes,nOfElementNodes); 
fxe = zeros(nOfElementNodes,1); 
fye = zeros(nOfElementNodes,1); 

%x and y coordinates of the nodes in the element
xe = Xe(:,1); ye = Xe(:,2);

%LOOP IN INTEGRATION POINTS
nIP= size(IPcoordinates,1); 
for g = 1:nIP
  %Shape functions and derivatives at the Gauss points
  N_g = N(g,:);  
  Nxi_g = Nxi(g,:);   
  Neta_g = Neta(g,:);
  %Jacobian of the isoparametric transformation
  J = [Nxi_g*xe	  Nxi_g*ye   
       Neta_g*xe  Neta_g*ye];
  %Integration point weight
  dvolu=IPweights(g)*det(J);
  %Derivatives (x,y) of the basis functions
  invJ = inv(J);  
  Nx_g = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
  Ny_g = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;
  %contribution to the mass matrix
  Me = Me + (N_g'*N_g)*dvolu; 
  %contribution to the r.h.s.
  fxe = fxe + N_g'*((Nx_g*ue)*dvolu); 
  fye = fye + N_g'*((Ny_g*ue)*dvolu); 
end 