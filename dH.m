function f = dH(d,alpha,beta)
% derivative of the soft repulsive potential
    f = alpha*beta*(sech(beta*d).^2)./2;
end