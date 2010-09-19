%{
  euro put at K_put
  euro call at K_call on each of the two worst performing stocks in 3-basket


  per stock output:
    - complete price path
    - complete delta path
    - complete delta hedge path
    - vector of call price
    - vector of put price

%}

% out = runRainbowSuspenders(param, prev)
function out = runRainbowSuspenders(varargin)

  optargs = {getParam() getEmptyPrev()};
  optargs(1:numel(varargin)) = varargin;
  [param prev] = optargs{:};
  
  getStat = @(x) [mean(x); std(x)/sqrt(numel(x))];

  N = param.T * 252;
  dt = param.T / N;
  m = numel(param.S);
  [Z out.z] = multivariate_gauss(param.rho, m, param.nSim, N, prev.z);

  for i=1:m
    % full stock path
    [x f_surface] = run_stock(param, prev, i, squeeze(Z(i,:,:)));
    
    % vanilla put and call
    callStrike = param.K_call*param.S(i).x_0;
    putStrike = param.K_put*param.S(i).x_0;
    out.S(i).call = exp(-param.r*N*dt) * max(x(:,end) - callStrike, 0);
    out.S(i).callstat = getStat(out.S(i).call);
    out.S(i).put = exp(-param.r*N*dt) * max(putStrike - x(:,end), 0);
    out.S(i).putstat = getStat(out.S(i).put);
    
    out.S(i).f_surface = f_surface;
    out.S(i).x = x;
    
    %{
    sorted_values = sort(out.S(i).x(:, end) ./ param.S(i).x_0 - 1);
    
    f = figure();
    subplot(2,1,1);
    [f,xi] = ksdensity(sorted_values, 'function', 'pdf');
    plot(xi, f, 'k');
    t = sprintf('Density of Stock Returns for %s', param.S(i).name);
    title(t);
    xlabel('Price');
    ylabel('Density');
    subplot(2,1,2);
    hist(sorted_values, 100);
    
    
    f = figure();
    m = mean(out.S(i).x);
    s = std(out.S(1).x);
    plot(1:numel(m), m);
    hold on;
    plot(1:numel(s), m+s, 'r:');
    plot(1:numel(s), m-s, 'r:');
    t = sprintf('Mean and Standard Deviation Over Time for %s', param.S(i).name);
    title(t);
    xlabel('Days');
    ylabel('Price');
    axis([1 numel(s) min(m-s) max(m+s)]);
    %}
    
    disp(i);
    
  end
  
  %{
  f = figure();
  scatterhist(out.S(1).x(:,end) / param.S(1).x_0 - 1, out.S(2).x(:, end) ./ param.S(2).x_0 - 1);
  t = sprintf('Stock Returns for %s v %s', param.S(1).name, param.S(2).name);
  title(t);
  xl = sprintf('%s', param.S(1).name);
  xlabel(xl);
  yl = sprintf('%s', param.S(2).name);
  ylabel(yl);
  
  f = figure();
  scatterhist(out.S(2).x(:,end) / param.S(2).x_0 - 1, out.S(3).x(:, end) ./ param.S(3).x_0 - 1);
  t = sprintf('Stock Returns for %s v %s', param.S(2).name, param.S(3).name);
  title(t);
  xl = sprintf('%s', param.S(2).name);
  xlabel(xl);
  yl = sprintf('%s', param.S(3).name);
  ylabel(yl);
  
  f = figure();
  scatterhist(out.S(3).x(:,end) / param.S(3).x_0 - 1, out.S(2).x(:, end) ./ param.S(2).x_0 - 1);
  t = sprintf('Stock Returns for %s v %s', param.S(3).name, param.S(1).name);
  title(t);
  xl = sprintf('%s', param.S(3).name);
  xlabel(xl);
  yl = sprintf('%s', param.S(1).name);
  ylabel(yl);
  %}
  
  % call on worst two stocks
  s1_returns = out.S(1).x ./ param.S(1).x_0;
  s2_returns = out.S(2).x ./ param.S(2).x_0;
  s3_returns = out.S(3).x ./ param.S(3).x_0;
  s1_ret = s1_returns(:, N/2);
  s2_ret = s2_returns(:, N/2);
  s3_ret = s3_returns(:, N/2);
  
  % find the 'best' return for each possible series
  maxes = max([s1_ret s2_ret s3_ret], [], 2);
  
  % compute whether it is the best or not.  If it is not, give a 1 index
  i1 = (s1_ret ~= maxes);
  i2 = (s2_ret ~= maxes);
  i3 = (s3_ret ~= maxes);
  
  % Stock values at T(2)
  e1 = out.S(1).x(:, end);
  e2 = out.S(2).x(:, end);
  e3 = out.S(3).x(:, end);
  
  % Stock values at T(1)
  h1 = out.S(1).x(:, N/2);
  h2 = out.S(2).x(:, N/2);
  h3 = out.S(3).x(:, N/2);
  
  % Determine the proportion that each stock gives to the basket
  n1 = i1.*(1.0 - h1 ./ (i1.*h1 + i2.*h2 + i3.*h3));
  n2 = i2.*(1.0 - h2 ./ (i1.*h1 + i2.*h2 + i3.*h3));
  n3 = i3.*(1.0 - h3 ./ (i1.*h1 + i2.*h2 + i3.*h3));
  
  % Compute the non-discounted value
  rainbow_v = max((n1.*e1 + n2.*e2 + n3.*e3) - (n1.*h1 + n2.*h2 + n3.*h3), 0);
  
  
  
  K1 = param.K_put * param.S(1).x_0;
  K2 = param.K_put * param.S(2).x_0;
  K3 = param.K_put * param.S(3).x_0;
  
  ks = [K1 K2 K3];
  
  p1 = max(K1 - out.S(1).x(:, end), 0);
  p2 = max(K2 - out.S(2).x(:, end), 0);
  p3 = max(K3 - out.S(3).x(:, end), 0);
  
  %{
  fprintf('Probability of being best performer after year 1\n');
  fprintf('\t%f,%f,%f\n', 1-mean(i1),1-mean(i2),1-mean(i3));
  %}
  out.rainbow = exp(-param.r*N*dt)*rainbow_v;
  
  holdFraction = param.notional / sum([param.S(:).sharesHeld] .* ks);
  putQty = holdFraction .* [param.S(:).sharesHeld];
  qt = min(putQty);
  putPrice = exp(-param.r*N*dt)*[mean(p1) mean(p2) mean(p3)];
  putLeg = sum(putPrice .* putQty);
  
  out.product = exp(-param.r*N*dt)*(-qt*p1 + -qt*p2 + -qt*p3 + qt*rainbow_v);
end