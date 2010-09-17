function [corr, standard_error] = calculate_correlation()
    h = 0.0001;
    
    params = getParam();
    params.rho = params.rho + h;
    out_up = runRainbowSuspenders(params);
    
    params = getParam();
    params.rho = params.rho - h;
    out_down = runRainbowSuspenders(params, out_up);
    
    n = numel(out_up.rainbow);
    curr_rhos = (out_up.product - out_down.product) / (2*h);
    corr = mean(curr_rhos);
    standard_error = sqrt(var(curr_rhos) / n);
end