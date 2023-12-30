function area=computeArea_elems(X,T,theReferenceElement)

if(nargin<=3)
    degree = 1; typeOfElement=1; %1=TRI, 0=QUA
    theReferenceElement = createReferenceElement(degree,typeOfElement);
end

IPweights = theReferenceElement.IPweights;
IPcoord = theReferenceElement.IPcoord;
N=theReferenceElement.N;
Nxi=theReferenceElement.Nxi;
Neta=theReferenceElement.Neta;

nOfElements = size(T,1);
area = zeros(nOfElements,1);
%Loop in elements
for i=1:nOfElements
    Te=T(i,:); Xe=X(Te,:); 
    xe = Xe(:,1); ye=Xe(:,2);
    for g=1:length(IPweights)
        N_g = N(g,:);
        Nxi_g = Nxi(g,:);
        Neta_g = Neta(g,:);
        J = [Nxi_g*xe	  Nxi_g*ye
            Neta_g*xe  Neta_g*ye];
        dvolu=IPweights(g)*det(J);
        Xg = N_g*Xe;
        area(i) =area(i) + dvolu;
    end      
end
%area = (area);

% nOfElements
% area
% 1/nOfElements
% sqrt(1/nOfElements)
% error('eng')