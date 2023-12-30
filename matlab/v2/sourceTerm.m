function res=sourceTerm(X)
x=X(:,1); y=X(:,2);
res=-dd_g(x).*g(y)-dd_g(y).*g(x);
end

function y=g(x) %g(x)
y = (4*x.*(1-x)).^10;
end

function res=dd_g(x)
%res=10*(4*x.*(1-x)).^8.*(144-608*x+640*x.*x);
res=10*(4*x.*(1-x)).^8.*(9*(4-8*x).^2-8*4*x.*(1-x));
end 

