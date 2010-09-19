function out = deltas_through_time(num_paths)
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
   
    for i = 1:num_paths
        
        fprintf('%i\n', i);
        
        current_path.S(1).x = data.original.S(1).x(i, :);
        current_path.S(2).x = data.original.S(2).x(i, :);
        current_path.S(3).x = data.original.S(3).x(i, :);
        
        out.S(1).path(i).stock_path = current_path.S(1).x;
        out.S(2).path(i).stock_path = current_path.S(2).x;
        out.S(3).path(i).stock_path = current_path.S(3).x;
        
        out.S(1).path(i).delta = zeros(N-1, 1);
        out.S(2).path(i).delta = zeros(N-1, 1);
        out.S(3).path(i).delta = zeros(N-1, 1);
        
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
                
                % make sure we don't get some weird artifacts from the
                % +h putting the product in our out of the money
                if today > 1
                    difference = abs(delta - out.S(j).path(i).delta(today-1));
                    if difference > 3e5;
                        out.S(j).path(i).delta(today) = out.S(j).path(i).delta(today-1);
                    else
                        out.S(j).path(i).delta(today) = delta;
                    end
                else
                    out.S(j).path(i).delta(today) = delta;
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
  
  % Stock values at T(2)
  e1 = d.s(1).x;
  e2 = d.s(2).x;
  e3 = d.s(3).x;
  
  % Stock values at T(1)
  h1 = d.s(1).strike;
  h2 = d.s(2).strike;
  h3 = d.s(3).strike;
  
  % Determine the proportion that each stock gives to the basket
  n1 = i1.*(1.0 - h1 ./ (i1.*h1 + i2.*h2 + i3.*h3));
  n2 = i2.*(1.0 - h2 ./ (i1.*h1 + i2.*h2 + i3.*h3));
  n3 = i3.*(1.0 - h3 ./ (i1.*h1 + i2.*h2 + i3.*h3));
  
  % Compute the non-discounted value
  rainbow_v = max((n1.*e1 + n2.*e2 + n3.*e3) - (n1.*h1 + n2.*h2 + n3.*h3), 0);
  
  K1 = param.K_put * param.S(1).x_0;
  K2 = param.K_put * param.S(2).x_0;
  K3 = param.K_put * param.S(3).x_0;
  
  p1 = max(K1 - e1, 0);
  p2 = max(K2 - e2, 0);
  p3 = max(K3 - e3, 0);
    
  qt1 = 1.6287e+006;
  qt2 = 3.8002e+006;
  qt3 = 2.1716e+006;
  
  qt4 = 2.7620e+006;
  
  v = mean(exp(-param.r*tau*param.dt)*(-qt1*p1 + -qt2*p2 + -qt3*p3 + qt4*rainbow_v));
end