function [prices, z] = run_put(s1, s2, s3, r, rho, z)
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
    gc_surface = local_vol_surface('Green Co.', prices, maturities, strikes, underlying, r);
    f_gc_surface = @(t,s) local_vol(gc_surface, t, s);
    s1 = euler_simulation(underlying, T, f_gc_surface, squeeze(Z(1,:,:)), r);

    prices = exp(-r*N*dt)*max(27.20 - s1(:, end), 0);
end

function sigma = local_vol(surface, t, s)
    t_index = t / surface.dT;
    s_index = s / surface.dK;
    
    % linearly interpolate between our surface points
    alpha = (s_index - floor(s_index)) / surface.dK;
    
    % clamp our volatility on the bottom and top level
    s_index(s_index < 1) = 1;
    s_index(s_index > 200) = 200;
    
    sigma_low = surface.surface(t_index, floor(s_index))';
    sigma_high = surface.surface(t_index, ceil(s_index))';
    
    % linearly interpolate
    sigma = alpha .* sigma_high + (1.0 - alpha) .* sigma_low;
end