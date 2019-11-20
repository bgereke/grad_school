function [x,y,beta]=generate_data(n, p, nz)
    x=randn(n, p);
    x = x * spdiags(1./sqrt(sum(x.^2))', 0, p, p);
    
    beta = zeros(p, 1);
    nz_idx = randsample(p, nz);
    beta(nz_idx) = 1;
    y = x*beta;
end