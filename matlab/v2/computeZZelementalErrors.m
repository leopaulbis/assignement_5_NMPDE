function elemErrors=computeZZelementalErrors(u,X,T,theReferenceElement)

if(nargin<=3)
    degree = 1; typeOfElement=1; %1=TRI, 0=QUA
    theReferenceElement = createReferenceElement(degree,typeOfElement);
end

grad = computeGradientSmoothing(u,X,T,theReferenceElement) ;

IPweights = theReferenceElement.IPweights;
IPcoord = theReferenceElement.IPcoord;
N=theReferenceElement.N;
Nxi=theReferenceElement.Nxi;
Neta=theReferenceElement.Neta;

nOfElements = size(T,1);
elemErrors = zeros(nOfElements,1);
%Loop in elements
for i=1:nOfElements
    Te=T(i,:); Xe=X(Te,:); 
    xe = Xe(:,1); ye=Xe(:,2);
    Ue=u(Te); 
    grade=grad(Te,:);
    for g=1:length(IPweights)
        N_g = N(g,:);
        Nxi_g = Nxi(g,:);
        Neta_g = Neta(g,:);
        J = [Nxi_g*xe	  Nxi_g*ye
            Neta_g*xe  Neta_g*ye];
        Gk = J\[Nxi_g;Neta_g];
        dvolu=IPweights(g)*det(J);
        Xg = N_g*Xe;

        dUe   = Gk*Ue;
        gradg = transpose(N_g*grade);
        elemErrors(i) = elemErrors(i) + (norm(gradg-dUe,2))^2*dvolu;
    end      
    elemErrors(i)=sqrt(elemErrors(i));
end
