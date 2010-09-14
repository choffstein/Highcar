function [prices, z] = run_simulation(s1, s2, s3, r, rho, z, fixed_sigma)
    T = 2; n=10000; N=252*2; m=3;
    dt = T / N;

    if nargin < 6
        [Z, z] = multivariate_gauss(rho, m, n, N);
    else
        [Z, z] = multivariate_gauss(rho, m, n, N, z);
    end
    
    if nargin < 7
        fixed_sigma = NaN;
    end
    
    q1 = 3000000;
    %underlying = 34;
    underlying = s1;
    prices = [4.99 6.09 7.44; ...
              1.58 2.79 4.05; ...
              0.22 0.85 1.72]';
    strikes = [30; 35; 40];
    maturities = [1/12; 3/12; 6/12];
    gc_surface = local_vol_surface('Green Co.', prices, maturities, strikes, underlying, r);
    if ~isnan(fixed_sigma)
        f_gc_surface = @(t,s) fixed_sigma;
    else
        f_gc_surface = @(t,s) local_vol(gc_surface, t, s);
    end
    s1 = euler_simulation(underlying, T, f_gc_surface, squeeze(Z(1,:,:)), r);
    s1_returns = s1 ./ underlying;
    
    q2 = 7000000;
    %underlying = 15;
    underlying = s2;
    prices = [5.10 5.31 5.66; ...
              1.15 1.87 2.47; ...
              0.06 0.39 0.91]';
    strikes = [10; 15; 20];
    maturities = [1/12; 3/12; 6/12];
    so_surface = local_vol_surface('SynerOptics', prices, maturities, strikes, underlying, r);
    
    if ~isnan(fixed_sigma)
        f_so_surface = @(t,s) fixed_sigma;
    else
        f_so_surface = @(t,s) local_vol(so_surface, t, s);
    end
    s2 = euler_simulation(underlying, T, f_so_surface, squeeze(Z(2,:,:)), r);
    s2_returns = s2 ./ underlying;
    
    q3 = 4000000;
    %underlying = 25;
    underlying = s3;
    prices = [5.35 6.03 6.89; ...
              1.60 2.53 3.44; ...
              0.16 0.63 1.31]';
    strikes = [20; 25; 30];
    maturities = [1/12; 3/12; 6/12];
    sw_surface = local_vol_surface('SW Industries', prices, maturities, strikes, underlying, r);
        
    if ~isnan(fixed_sigma)
        f_sw_surface = @(t,s) fixed_sigma;
    else
        f_sw_surface = @(t,s) local_vol(sw_surface, t, s);
    end
    s3 = euler_simulation(underlying, T, f_sw_surface, squeeze(Z(3,:,:)), r);
    s3_returns = s3 ./ underlying;
    
    s1_ret = s1_returns(:, N/2);
    s2_ret = s2_returns(:, N/2);
    s3_ret = s3_returns(:, N/2);
        
    % find the 'best' return for each possible series
    maxes = max([s1_ret s2_ret s3_ret], [], 2);
    
    % compute whether it is the best or not.  If it is not, give a 1 index
    i1 = (s1_ret == maxes);
    i2 = (s2_ret == maxes);
    i3 = (s3_ret == maxes);
    
    s1_ret = s1_returns(:, end);
    s2_ret = s2_returns(:, end);
    s3_ret = s3_returns(:, end);
            
    v = max((i1.*s1_ret + i2.*s2_ret + i3.*s3_ret)-1, 0);
    prices = exp(-r*N*dt)*v;
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