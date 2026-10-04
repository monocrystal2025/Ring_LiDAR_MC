function [low,high] = wilson(successes,n,alpha)
%WILSON Two-sided binomial interval; alpha can include Bonferroni correction.
arguments
    successes double
    n (1,1) double {mustBePositive}
    alpha (1,1) double {mustBeGreaterThan(alpha,0),mustBeLessThan(alpha,1)} = 0.05
end
z = sqrt(2)*erfcinv(alpha);
p = successes/n;
denominator = 1+z^2/n;
center = (p+z^2/(2*n))/denominator;
half = z*sqrt(p.*(1-p)/n+z^2/(4*n^2))/denominator;
low = max(0,center-half);
high = min(1,center+half);
end
