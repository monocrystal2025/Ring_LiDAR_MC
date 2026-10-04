function [low,high] = binomialBounds(successes,n,alpha)
%BINOMIALBOUNDS Exact Clopper-Pearson interval without a statistics toolbox.
arguments
    successes double
    n (1,1) double {mustBeInteger,mustBePositive}
    alpha (1,1) double {mustBeGreaterThan(alpha,0),mustBeLessThan(alpha,1)} = 0.05
end
successes = round(successes);
low = zeros(size(successes));
high = ones(size(successes));
insideLow = successes>0;
insideHigh = successes<n;
low(insideLow) = betaincinv(alpha/2,successes(insideLow),n-successes(insideLow)+1);
high(insideHigh) = betaincinv(1-alpha/2,successes(insideHigh)+1,n-successes(insideHigh));
end
