lambda = 1:100;
m = linspace(0, 0.5, 51);
[M, L] = meshgrid(m, lambda);
MU = zeros(size(M)); S = MU;
for i = 1 : length(m)
    for j = 1 : length(lambda)
	[mu, s] = intGaussian(m(i), lambda(j));
	MU(j, i) = mu - m(i); S(j, i) = s*lambda(j);
    end
end
figure(1); mesh(M, L, MU);
figure(2); mesh(M, L, S);


function [mu, s] = intGaussian(m, lambda)
    m_int = round(m);
    m = m - m_int;
    n = -30:30;
    weight = exp(-lambda/2*(n-m).^2);
    weight = weight/sum(weight);
    
    mu = sum(n.*weight);
    s = sum(n.^2.*weight) - mu^2;
    mu = mu + m_int;
end
