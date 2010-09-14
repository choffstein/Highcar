function [rho, standard_error] = calculate_rho()
    h = 0.0001;
    
    [v1, z] = run_simulation(34, 15, 25, 0.045+h, 0.3);
    v2 = run_simulation(34, 15, 25, 0.045-h, 0.3, z);
    
    n = numel(v1(:, end));
    rhos = (v1(:, end) - v2(:, end)) / (2*h);
    rho = mean(rhos);
    standard_error = sqrt(var(rhos) / n);
end

