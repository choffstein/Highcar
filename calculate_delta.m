function [deltas, standard_errors] = calculate_delta()
    h = 0.0001;
    deltas = zeros(3, 1);
    standard_errors = zeros(3, 1);

    [v1, z] = run_simulation(34+h, 15, 25, 0.045, 0.3);
    v2 = run_simulation(34-h, 15, 25, 0.045, 0.3, z);
    
    n = numel(v1(:, end));
    curr_deltas = (v1(:, end) - v2(:, end)) / (2*h);
    deltas(1) = mean(curr_deltas);
    standard_errors(1) = sqrt(var(curr_deltas) / n);

    
    v1 = run_simulation(34, 15+h, 25, 0.045, 0.3, z);
    v2 = run_simulation(34, 15-h, 25, 0.045, 0.3, z);
    
    n = numel(v1(:, end));
    curr_deltas = (v1(:, end) - v2(:, end)) / (2*h);
    deltas(2) = mean(curr_deltas);
    standard_errors(2) = sqrt(var(curr_deltas) / n);

    
    v1 = run_simulation(34, 15, 25+h, 0.045, 0.3, z);
    v2 = run_simulation(34, 15, 25-h, 0.045, 0.3, z);
    
    n = numel(v1(:, end));
    curr_deltas = (v1(:, end) - v2(:, end)) / (2*h);
    deltas(3) = mean(curr_deltas);
    standard_errors(3) = sqrt(var(curr_deltas) / n);
end