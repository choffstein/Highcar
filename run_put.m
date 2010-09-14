function [prices, z, f_gc_surface] = run_put(s1, s2, s3, r, rho, z, f_gc_surface)
    T = 2; n=10000; N=252*2; m=3;
    dt = T / N;

    if nargin < 6
        [Z, z] = multivariate_gauss(rho, m, n, N);
    else
        [Z, z] = multivariate_gauss(rho, m, n, N, z);
    end

    underlying = s1;
    prices = [4.99 6.09 7.44; ...
              1.58 2.79 4.05; ...
              0.22 0.85 1.72]';
    strikes = [30; 35; 40];
    maturities = [1/12; 3/12; 6/12];
    
    if nargin < 7
        gc_surface = local_vol_surface('Green Co.', prices, maturities, strikes, underlying, r);
        f_gc_surface = @(t,s) local_vol(gc_surface, t, s);
    end
    
    s1 = euler_simulation(underlying, T, f_gc_surface, squeeze(Z(1,:,:)), r);

    prices = exp(-r*N*dt)*max(34 - s1(:, end), 0);
end

function sigma = local_vol(surface, t, s)
    t_index = floor(t / surface.dT);
    s_index = floor(s / surface.dK);
    
    % clamp our volatility on the bottom and top level
    t_index(t_index < 1) = 1;
    t_index(t_index > numel(surface.maturities)) = numel(surface.maturities);
    s_index(s_index < 1) = 1;
    s_index(s_index > numel(surface.strikes)) = numel(surface.strikes);

    sigma = surface.surface(t_index, s_index)';
end