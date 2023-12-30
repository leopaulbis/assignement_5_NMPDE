function u=analytical(X)

u =h(X(:,1)).*h(X(:,2));


function y=h(x)
y = (4*x.*(1-x)).^10;


