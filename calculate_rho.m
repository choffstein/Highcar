function [rho, standard_error] = calculate_rho()
    h = 0.0001;

    params = getParam();
    params.r = params.r + h;
    out_up = runRainbowSuspenders(params);
    
    params = getParam();
    params.r = params.r - h;
    out_down = runRainbowSuspenders(params, out_up);
    
    n = numel(out_up.product);
    curr_rhos = (out_up.product - out_down.product) / (2*h);
    rho = mean(curr_rhos);
    standard_error = sqrt(var(curr_rhos) / n);
end

