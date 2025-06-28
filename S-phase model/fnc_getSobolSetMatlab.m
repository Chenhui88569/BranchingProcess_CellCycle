function X = fnc_getSobolSetMatlab(dim, N)
rng('shuffle')
p = sobolset(dim);
p = scramble(p,'MatousekAffineOwen');
X = net(p,N);