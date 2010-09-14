function [model_risks, standard_errors] = calculate_model_risk()
    sigmas = 0.1:0.1:1;
    model_risks = zeros(numel(sigmas), 1);
    standard_errors = zeros(numel(sigmas), 1);
    
    [v1, z] = run_simulation(34, 15, 25, 0.045, 0.3);
    for i = 1:numel(sigmas)
        sigma = sigmas(i);
       
        v2 = run_simulation(34, 15, 25, 0.045, 0.3, z, sigma);
    
        n = numel(v1(:, end));
        mr = v1(:, end) - v2(:, end);
        model_risks(i) = mean(mr);
        standard_errors(i) = sqrt(var(mr) / n);
    end
end