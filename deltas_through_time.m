function delta_path = deltas_through_time(num_paths)
    h = 0.0001;
    
    params = getParam();
    data.original = runRainbowSuspenders();
    N = params.T * 252;

    for i = 1:3
        data.S(i).params.up = getParam();
        data.S(i).params.up.S(i).x_0 = data.S(i).params.up.S(i).x_0 + h;
        data.S(i).up = runRainbowSuspenders(data.S(i).params.up, data.original);

        data.S(i).params.down = getParam();
        data.S(i).params.down.S(i).x_0 = data.S(i).params.down.S(i).x_0 - h;
        data.S(i).down = runRainbowSuspenders(data.S(i).params.down, data.original);
       
        data.S(i).params.up.dt = params.T / N;
        data.S(i).params.down.dt = params.T / N;
    end
    
    delta_path = zeros(num_paths, N, 3);
    for i = 1:num_paths
        
        fprintf('%i\n', i);
        
        current_path.S(1).x = data.original.S(1).x(i, :);
        current_path.S(2).x = data.original.S(2).x(i, :);
        current_path.S(3).x = data.original.S(3).x(i, :);
        
            
        for today = 1:(N-1)
            for j = 1:3
                current_data.up = data.S(j).up;
                current_data.down = data.S(j).down;
                
                for k = 1:3
                    [up.S(j).s(k).x, up.S(j).s(k).strike] = cut_and_normalize(current_data.up.S(k).x(:, today), ...
                                        current_data.up.S(k).x(:, N/2), ...
                                        current_path.S(k).x(today), ...
                                        current_data.up.S(k).x(:, end));
                    
                    [down.S(j).s(k).x, down.S(j).s(k).strike] = cut_and_normalize(current_data.down.S(k).x(:, today), ...
                                        current_data.down.S(k).x(:, N/2), ...
                                        current_path.S(k).x(today), ...
                                        current_data.down.S(k).x(:, end));
                end
                
                % compute the value of the product in the up-state
                % and the down state
                value.up = product_value(up.S(j), data.S(j).params.up, N-today);
                value.down = product_value(down.S(j), data.S(j).params.down, N-today);
                
                delta = (value.up - value.down) / (2*h);
                delta_path(i, today, j) = delta;
                
                if abs(delta) > 1e7
                    5;
                end
            end
        end
    end
end

% normalize the simulation path -- reset the strike value and final value
function [m, ok] = cut_and_normalize(A, k, x, f)
    v = x ./ A; %normalizing factor
    m = v .* f;  
    ok = v .* k;     
end

function v = product_value(d, param, tau)
  % call on worst two stocks
  s1_ret = d.s(1).strike ./ param.S(1).x_0;
  s2_ret = d.s(2).strike ./ param.S(2).x_0;
  s3_ret = d.s(3).strike ./ param.S(3).x_0;
  
  % find the 'best' return for each possible series
  
  maxes = max([s1_ret s2_ret s3_ret], [], 2);
  
  % compute whether it is the best or not.  If it is not, give a 1 index
  i1 = (s1_ret ~= maxes);
  i2 = (s2_ret ~= maxes);
  i3 = (s3_ret ~= maxes);
  
  rainbow_v = max((i1.*d.s(1).x + i2.*d.s(2).x + ...
              i3.*d.s(3).x) - (i1.*d.s(1).strike + ...
              i2.*d.s(2).strike + i3.*d.s(3).strike), 0);

  
  K1 = param.K_put * param.S(1).x_0;
  K2 = param.K_put * param.S(2).x_0;
  K3 = param.K_put * param.S(3).x_0;
  
  ks = [K1 K2 K3];
  
  p1 = max(K1 - d.s(1).x, 0);
  p2 = max(K2 - d.s(2).x, 0);
  p3 = max(K3 - d.s(3).x, 0);
  
  out.rainbow = exp(-param.r*tau*param.dt)*rainbow_v;
  
  holdFraction = param.notional / sum([param.S(:).sharesHeld] .* ks);
  putQty = holdFraction .* [param.S(:).sharesHeld];
  qt = min(putQty);
  
  v = mean(exp(-param.r*tau*param.dt)*(-qt*p1 + -qt*p2 + -qt*p3 + qt*rainbow_v));
end