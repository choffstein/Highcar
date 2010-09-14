function [VaR, CVaR] = calculate_VaR(p)
    values = run_simulation(34, 15, 25, 0.045, 0.3);
    final_values = values(:, end);
    index = floor((1-p) * numel(final_values));
    
    sorted_values = sort(final_values);
    
    mean_value = mean(sorted_values);
    median_value = median(sorted_values);
    VaR = sorted_values(index);
    CVaR = mean(sorted_values(1:index));
    
    f=figure();
    subplot(2,1,1);
    [f,xi] = ksdensity(sorted_values);
    indices = find(xi <= VaR);
    index = indices(end);
    plot(xi, f, 'k');
    hold on;
    area(xi(1:index), f(1:index));
    indices = find(xi <= CVaR);
    index = indices(end);
    plot([xi(index) xi(index)], [0 1], 'r--');
    
    indices = find(xi <= mean_value);
    index = indices(end);
    plot([xi(index) xi(index)], [0 1], 'g--');
    
    indices = find(xi <= median_value);
    index = indices(end);
    plot([xi(index) xi(index)], [0 1], 'b--');
    
    axis([0 xi(end) 0 max(f)]);
    title('Probability Distribution Function');
    ylabel('f(x)');
    xlabel('Value of Product at Maturity');
    
    subplot(2,1,2);
    [f,xi] = ksdensity(sorted_values, 'function', 'cdf');
    indices = find(xi < VaR);
    index = indices(end);
    plot(xi, f, 'k');
    hold on;
    area(xi(1:index), f(1:index));
    indices = find(xi <= mean_value);
    index = indices(end);
    plot([xi(index) xi(index)], [0 1], 'g--');
    
    indices = find(xi <= median_value);
    index = indices(end);
    plot([xi(index) xi(index)], [0 1], 'b--');
    indices = find(xi <= CVaR);
    index = indices(end);
    plot([xi(index) xi(index)], [0 1], 'r--');
    axis([0 xi(end) 0 1]);
    title('Cumulative Distribution Function');
    ylabel('F(x)'); 
    xlabel('Value of Product at Maturity');
end