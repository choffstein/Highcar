function [gammas, standard_errors] = calculate_gamma()
    h = 0.0001;
    gammas = zeros(3, 1);
    standard_errors = zeros(3, 1);

    params = getParam();
    out_original = runRainbowSuspenders(params);
    
    params = getParam();
    params.S(1).x_0 = params.S(1).x_0 + h;
    out_up = runRainbowSuspenders(params, out_original);
    
    params = getParam();
    params.S(1).x_0 = params.S(1).x_0 - h;
    out_down = runRainbowSuspenders(params, out_original);
    
    n = numel(out_original.product);
    curr_gammas = (out_up.product - 2*out_original.product + out_down.product) ./ (h^2);
    gammas(1) = mean(curr_gammas);
    standard_errors(1) = sqrt(var(curr_gammas) / n);
    
    
    
    
    params = getParam();
    params.S(2).x_0 = params.S(2).x_0 + h;
    out_up = runRainbowSuspenders(params, out_original);
    
    params = getParam();
    params.S(2).x_0 = params.S(2).x_0 - h;
    out_down = runRainbowSuspenders(params, out_original);
    
    n = numel(out_original.product);
    curr_gammas = (out_up.product - 2*out_original.product + out_down.product) ./ (h^2);
    gammas(2) = mean(curr_gammas);
    standard_errors(2) = sqrt(var(curr_gammas) / n);
    
    
    
    
    
    params = getParam();
    params.S(3).x_0 = params.S(3).x_0 + h;
    out_up = runRainbowSuspenders(params, out_original);
    
    params = getParam();
    params.S(3).x_0 = params.S(3).x_0 - h;
    out_down = runRainbowSuspenders(params, out_original);
    
    n = numel(out_original.product);
    curr_gammas = (out_up.product - 2*out_original.product + out_down.product) ./ (h^2);
    gammas(3) = mean(curr_gammas);
    standard_errors(3) = sqrt(var(curr_gammas) / n);
end