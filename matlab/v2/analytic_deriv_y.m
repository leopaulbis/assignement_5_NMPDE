
function res=analytic_deriv_y(X)
x=X(:,1); y=X(:,2);
res=g(x)*dg(y);
end

function res=g(x)
res=(4*x.*(1-x)).^10;
end 

function res=dg(x)
res=10*(4-8*x).*(4*x.*(1-x)).^9;
end
