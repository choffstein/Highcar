function result = local_vol_surface(S, r)

    name = S.name;
    option_prices = S.c;
    maturities = S.T;
    strikes = S.K;
    underlying = S.x_0;


    [n, m] = size(option_prices);
   
    % Construct the implied volatility surface
    phi = 1; % given calls
    ivol_surface = zeros(n, m);
    vegas = zeros(n, m);
    for i = 1:n
        for j = 1:m
            ivol_surface(i,j) = max(0, bsm.ivol(phi, option_prices(i,j), ...
                        underlying, strikes(j), maturities(i), 0, r, 0));
             
            vegas(i, j) = bsm.vega(underlying, strikes(j), ...
                            ivol_surface(i, j), maturities(i), 0, r, 0);
        end
    end
    
    fHandle = @(parameters) surface_handle(ivol_surface, maturities, ...
                                    strikes, vegas, parameters);
    guess = randn(numel(maturities), 6);
    options = optimset('MaxFunEvals', Inf, 'TolX', 1e-16, ...
                       'TolFun', 1e-16, 'Display', 'off', ...
                       'Diagnostics', 'off', 'MaxIter', Inf, ...
                        'LargeScale', 'off');
    [parameters, mse] = fminunc(fHandle, guess, options);
                 
    dT = 1/252;
    To = maturities(1):dT:maturities(end);
    dK = 1;
    Ko = strikes(1):dK:strikes(end);
    
    a = parameters(1);
    b = parameters(2);
    c = parameters(3);
    d = parameters(4);
    e = parameters(5);
    f = parameters(6);
    
    smoothed_ivol = zeros(numel(To), numel(Ko));
    for i = 1:numel(To)
        t = To(i);
        
        for j = 1:numel(Ko)
            k = Ko(j);
            linear_shift = a;
            slope_change = b*(t - mean(maturities)) + c*(k - mean(strikes));
            curvature = d*(t - mean(maturities))^2 + ...
                        e*(k - mean(strikes))^2 + ...
                        f*(t - mean(maturities))*(k - mean(strikes));
                    
            smoothed_ivol(i,j) = linear_shift + slope_change + curvature;
        end
    end
    
    %{
    f = figure();
    subplot(2,1,1);
    surf(strikes, maturities, ivol_surface);
    colorbar;
    alpha(.4);
    xlabel('Strikes');
    ylabel('TTM');
    zlabel('Implied Volatility (%)');
    local_header = sprintf('%s Implied Volatility', name);
    title(local_header);

    subplot(3,1,2);
    surf(extended_strikes, extended_maturities, extended_ivol_surface);
    colorbar;
    alpha(.4);
    xlabel('Strikes');
    ylabel('TTM');
    zlabel('Implied Volatility (%)');
    local_header = sprintf('%s Extended Implied Volatility', name);
    title(local_header);

    
    subplot(2,1,2);
    surf(Ko, To, smoothed_ivol);
    colorbar;
    alpha(.4);
    xlabel('Strikes');
    ylabel('TTM');
    zlabel('Implied Volatility (%)');
    local_header = sprintf('%s Interpolated & Extrapolated Implied Volatility', name);
    title(local_header);
    %}
    
    % To construct our local vol models we have to numerically find our 
    % derivatives (we use forward first differences)
    dVol_dT = zeros(size(smoothed_ivol));
    dVol_dT(1:end-1, :) = (smoothed_ivol(2:end, :) - smoothed_ivol(1:end-1, :)) ./ dT;
    
    % cap the edges via replication
    dVol_dT(:, 1) = dVol_dT(:,2);
    dVol_dT(1, :) = dVol_dT(2,:);
    dVol_dT(:, end) = dVol_dT(:, end-1);
    dVol_dT(end, :) = dVol_dT(end-1, :); 
    
    dVol_dK = zeros(size(smoothed_ivol));
    dVol_dK(:, 1:end-1) = (smoothed_ivol(:, 2:end) - smoothed_ivol(:, 1:end-1)) ./ dK;
    dVol_dK(:, 1) = dVol_dK(:,2);
    dVol_dK(1, :) = dVol_dK(2,:);
    dVol_dK(:, end) = dVol_dK(:, end-1);
    dVol_dK(end, :) = dVol_dK(end-1, :); 
    
    % second finite difference
    d2Vol_dK2 = zeros(size(smoothed_ivol));
    d2Vol_dK2(:, 2:end-1) = (dVol_dK(:, 3:end) - 2*dVol_dK(:, 2:end-1) + ...
                                            dVol_dK(:, 1:end-2)) ./ dK^2;
    d2Vol_dK2(:, 1) = d2Vol_dK2(:,2);
    d2Vol_dK2(1, :) = d2Vol_dK2(2,:);
    d2Vol_dK2(:, end) = d2Vol_dK2(:, end-1);
    d2Vol_dK2(end, :) = d2Vol_dK2(end-1, :); 
    
    % Compute our local volatility using Dupire
    % See http://www.fincad.com/derivatives-resources/articles/local-volatility.aspx
    local_vol = zeros(size(smoothed_ivol));
    for i = 1:numel(To)
        for j = 1:numel(Ko)
            sigma = smoothed_ivol(i,j);
            
            tau = To(i);
            k = Ko(j);
            
            d1 = (log(underlying / k) + (r +0.5*sigma^2)*tau) / (sigma*sqrt(tau));
            dSdT = dVol_dT(i,j);
            dSdK = dVol_dK(i,j);
            d2SdK2 = d2Vol_dK2(i,j);
            
            numerator = sigma^2 + 2*tau*sigma*dSdT + ...
                                2*r*k*tau*sigma*dSdK;
            
            denominator = (1 + k*d1*sqrt(tau)*dSdK)^2 + ...
                                (k^2)*tau*sigma*(d2SdK2 - ...
                                        d1*(dSdK^2)*sqrt(tau));
            
            %local vol cutoff at 2 percent (don't go non-negative)
            local_vol(i, j) = max(sqrt(numerator / denominator), 0.02);
        end 
    end
    
    result.surface = local_vol;
    result.maturities = To;
    result.strikes = Ko;
    result.dT = dT;
    result.dK = dK;
    
    %{
    f = figure();
    surf(Ko, To, local_vol);
    colorbar;
    alpha(.4);
    xlabel('Strikes');
    ylabel('TTM');
    zlabel('Local Volatility (%)');
    local_header = sprintf('%s Local Volatility', name);
    title(local_header);
    %}
end

function mse = surface_handle(surface, maturities, strikes, weights, parameters)
    [n m] = size(surface);
    
    fit_surface = zeros(size(surface));
    
    a = parameters(1);
    b = parameters(2);
    c = parameters(3);
    d = parameters(4);
    e = parameters(5);
    f = parameters(6);
    
    for i = 1:n
        t = maturities(i);
            
        for j = 1:m
    
            k = strikes(j);
            linear_shift = a;
            slope_change = b*(t - mean(maturities)) + c*(k - mean(strikes));
            curvature = d*(t - mean(maturities))^2 + ...
                        e*(k - mean(strikes))^2 + ...
                        f*(t - mean(maturities))*(k - mean(strikes));
                    
            fit_surface(i,j) = linear_shift + slope_change + curvature;
        end
    end
    
    mse = norm(weights .* (fit_surface - surface), 'fro');
end

