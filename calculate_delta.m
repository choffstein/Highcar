function [deltas, standard_errors] = calculate_delta()
    h = 0.0001;
    deltas = zeros(3, 1);
    standard_errors = zeros(3, 1);
    
    params = getParam();
    params.S(1).x_0 = params.S(1).x_0 + h;
    out_up = runRainbowSuspenders(params);
    
    params = getParam();
    params.S(1).x_0 = params.S(1).x_0 - h;
    out_down = runRainbowSuspenders(params, out_up);
    
    n = numel(out_up.product);
    curr_deltas = (out_up.product - out_down.product) / (2*h);
    deltas(1) = mean(curr_deltas);
    standard_errors(1) = sqrt(var(curr_deltas) / n);
    
    
    
    
    params = getParam();
    params.S(2).x_0 = params.S(2).x_0 + h;
    out_up = runRainbowSuspenders(params);
    
    params = getParam();
    params.S(2).x_0 = params.S(2).x_0 - h;
    out_down = runRainbowSuspenders(params, out_up);
    
    n = numel(out_up.product);
    curr_deltas = (out_up.product - out_down.product) / (2*h);
    deltas(2) = mean(curr_deltas);
    standard_errors(2) = sqrt(var(curr_deltas) / n);
    
    
    
    
    params = getParam();
    params.S(3).x_0 = params.S(3).x_0 + h;
    out_up = runRainbowSuspenders(params);
    
    params = getParam();
    params.S(3).x_0 = params.S(3).x_0 - h;
    out_down = runRainbowSuspenders(params, out_up);
    
    n = numel(out_up.product);
    curr_deltas = (out_up.product - out_down.product) / (2*h);
    deltas(3) = mean(curr_deltas);
    standard_errors(3) = sqrt(var(curr_deltas) / n);
end