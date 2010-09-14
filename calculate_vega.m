function [vega, standard_error] = calculate_vega()
    h = 0.0001;
    
    [v1, Z] = run_simulation(34, 15, 25, 0.045, 0.3, [], 0.35+h);
    v2 = run_simulation(34, 15, 25, 0.045, 0.3, Z, 0.35-h);
    
    n = numel(v1(:, end));
    rhos = (v1(:, end) - v2(:, end)) / (2*h);
    vega = mean(rhos);
    standard_error = sqrt(var(rhos) / n);
end