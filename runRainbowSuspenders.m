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
function out = runRainbowSuspenders()

  param = getParam();
  prev.z = [];
  prev.S(1).f_gc_surface = [];
  prev.S(2).f_gc_surface = [];
  prev.S(3).f_gc_surface = [];
  
  
  getStat = @(x) [mean(x); std(x)/sqrt(numel(x))];

  N = param.T * 252;
  dt = param.T / N;
  m = numel(param.S);
  [Z out.z] = multivariate_gauss(param.rho, m, param.nSim, N, prev.z);

  for i=[1 3]
    % full stock path
    [x f_gc_surface] = run_stock(param, prev, i, squeeze(Z(i,:,:)));
    
    %
    
    
    % vanilla put and call
    callStrike = param.K_call*param.S(i).x_0;
    putStrike = param.K_put*param.S(i).x_0;
    out.S(i).call = exp(-param.r*N*dt) * max(x(:,end) - callStrike, 0);
    out.S(i).callstat = getStat(out.S(i).call);
    out.S(i).put = exp(-param.r*N*dt) * max(putStrike - x(:,end), 0);
    out.S(i).putstat = getStat(out.S(i).put);
    
    % call on worst two stocks
    
    out.S(i).f_gc_surface = f_gc_surface;
    out.S(i).x = x;
  end
  
  
  
  
  %{
  
    RAINBOW CALL
  
    s1_returns = out.S(i).x ./ param.S(i).x_0;
  
    s1_ret = s1_returns(:, N/2);
    s2_ret = s2_returns(:, N/2);
    s3_ret = s3_returns(:, N/2);
        
    % find the 'best' return for each possible series
    maxes = max([s1_ret s2_ret s3_ret], [], 2);
    
    % compute whether it is the best or not.  If it is not, give a 1 index
    i1 = (s1_ret == maxes);
    i2 = (s2_ret == maxes);
    i3 = (s3_ret == maxes);
    
    s1_ret = s1_returns(:, end);
    s2_ret = s2_returns(:, end);
    s3_ret = s3_returns(:, end);
            
    v = max((i1.*s1_ret + i2.*s2_ret + i3.*s3_ret)-1, 0);
    prices = exp(-r*N*dt)*v;
  %}


end