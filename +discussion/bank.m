function bank = bank(n,seed)
%BANK Paired random scenarios; does not change MATLAB's global RNG.
arguments
    n (1,1) double {mustBeInteger,mustBePositive}
    seed (1,1) double {mustBeInteger,mustBeNonnegative}
end
stream = RandStream('mt19937ar','Seed',seed);
bank = struct('N',n,'seed',seed,'u',rand(stream,n,8));
end
