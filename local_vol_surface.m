function result = local_vol_surface(name, option_prices, maturities, strikes, underlying, r)
    [n, m] = size(option_prices);
   
    % Construct the implied volatility surface
    ivol_surface = zeros(n, m);
    for i = 1:n
        for j = 1:m
            ivol_surface(i,j) = max(0, bsm.ivol(option_prices(i,j), ...
                        underlying, strikes(j), maturities(i), 0, r, 0));
        end
    end
    
    % Extend the ivol surface for when we extrapolate to try to
    % asymptotically approach a constant value
    extended_ivol_surface = zeros(n+2, m+2);
    extended_ivol_surface(2:end-1, 2:end-1) = ivol_surface;
    
    % down the left side
    extended_ivol_surface(1,1) = ivol_surface(1,1);
    extended_ivol_surface(2:end-1, 1) = ivol_surface(:, 1);
    extended_ivol_surface(end, 1) = ivol_surface(end, 1);
    
    % on top
    extended_ivol_surface(1, 2:end-1) = ivol_surface(1, :);
    extended_ivol_surface(1, end) = ivol_surface(1, end);
    
    % on the right
    extended_ivol_surface(2:end-1, end) = ivol_surface(:, end);
    extended_ivol_surface(end, end) = ivol_surface(end, end);
    
    % on the bottom
    extended_ivol_surface(end, 2:end-1) = ivol_surface(end, :);
    
    %%%% CHANGE BACK TO 1/252
    %%%% and dK = 1
    
    
    % define what our extended maturities and strikes are
    extended_maturities = [1/252; maturities; 2];
    extended_strikes = [1; strikes; 200];
    
    % Smoothing the surface using thin-plate splines
    dK = 1;
    Ko = 1:dK:extended_strikes(end);
    dT = 1/252;
    To = 1/252:dT:extended_maturities(end);
    
    NKo = length(Ko); NTo = length(To);
    
    [Kgo, Tgo] = meshgrid(Ko, To);
    
    Kvo = reshape(Kgo, [NKo*NTo, 1]); % vectorize
    Tvo = reshape(Tgo, [NKo*NTo, 1]); 
    
    Tmo = log(Tvo);
    Kmo = log(underlying ./ Kvo) ./ sqrt(Tvo);
    Km = log(underlying ./ extended_strikes) ./ sqrt(extended_maturities);
    Tm = log(extended_maturities);

    sbm = zeros(numel(Km) * numel(Tm), 2);
    for k = 1:numel(Km)
        for t = 1:numel(Tm)
            sbm((k-1) * numel(Km) + t, :) = [Km(k) Tm(t)];
        end
    end
    
    % perform the actual smoothing (Radial Basis Function with Thinplate
    % kernel)
    coef = rbfcreate(sbm', extended_ivol_surface(:)', ...
                    'RBFFunction', 'multiquadric', 'RBFSmooth', 75);
    IVvo = rbfinterp([Kmo'; Tmo'], coef)';
    smoothed_ivol = reshape(IVvo, [NTo, NKo]);
    
    % ensure our volatility does not go negative.
    smoothed_ivol(smoothed_ivol < 0) = 0.02; 
    
    %{
    f = figure();
    surf(strikes, maturities, ivol_surface);
    colorbar;
    alpha(.4);
    xlabel('Strikes');
    ylabel('TTM');
    zlabel('Implied Volatility (%)');
    local_header = sprintf('%s Implied Volatility', name);
    title(local_header);
    
    f = figure();
    surf(extended_strikes, extended_maturities, extended_ivol_surface);
    colorbar;
    alpha(.4);
    xlabel('Strikes');
    ylabel('TTM');
    zlabel('Implied Volatility (%)');
    local_header = sprintf('%s Extended Implied Volatility', name);
    title(local_header);
    
    f = figure();
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
    for i = 1:NTo
        for j = 1:NKo
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
                                        d1*dSdK^2*sqrt(tau));
            %local vol cutoff at 2 percent (don't go non-negative)
            local_vol(i, j) = sqrt(max(numerator / denominator, 0.02));
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

