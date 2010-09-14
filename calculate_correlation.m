function [corr, standard_error] = calculate_correlation()
    h = 0.0001;
    [v1, z] = run_simulation(34, 15, 25, 0.045, 0.3-h);
    v2 = run_simulation(34, 15, 25, 0.045, 0.3+h, z);
    
    n = numel(v1(:, end));
    rhos = (v1(:, end) - v2(:, end)) / (2*h);
    corr = mean(rhos);
    standard_error = sqrt(var(rhos) / n);
end