function hr = hedgeReport(data)

  p = getParam();
  % total cost of puts
  % qty purchased per day
  % delta of each put
  % 2yr ivol
  
  % number of rainbows sold
  % cost per rainbow
  % delta per rainbow
  
  getMcStat = @(x) [mean(x) std(x)/sqrt(numel(x))];
  
  
  hr.rainbow.v = getMcStat(data.rainbow);
  
  hr.S(3).call = [];
  for i=1:numel(data.S)
    K_put = p.K_put * p.S(i).x_0;
    K_call = p.K_call * p.S(i).x_0;
    
    hr.S(i).put.v = getMcStat(data.S(i).put);
    hr.S(i).put.K = K_put;
    hr.S(i).put.ivol = bsm.ivol(-1, hr.S(i).put.v(1), p.S(i).x_0, K_put, p.T, 0, p.r);
    hr.S(i).put.delta_0 = bsm.forwardDelta(-1, p.S(i).x_0, K_put, hr.S(i).put.ivol, p.T, 0, p.r);
    hr.S(i).call.v = getMcStat(data.S(i).call);
    hr.S(i).call.K = K_call;
    hr.S(i).call.ivol = bsm.ivol(1, hr.S(i).call.v(1), p.S(i).x_0, K_call, p.T, 0, p.r);
    hr.S(i).call.delta_0 = bsm.forwardDelta(1, p.S(i).x_0, K_call, hr.S(i).call.ivol, p.T, 0, p.r); 
  end
  
  holdFraction = p.notional / sum([p.S(:).sharesHeld] .* arrayfun(@(x) x.put.K, hr.S));
  
  putQty = holdFraction .* [p.S(:).sharesHeld];
  putPrice = arrayfun(@(x) x.put.v(1), hr.S);
  
  putDelta = arrayfun(@(x) x.put.delta_0, hr.S);
  callDelta = arrayfun(@(x) x.call.delta_0, hr.S);
  maxTrade = [p.S.adv] * .20;
  
  requiredDelta = callDelta - putDelta;
  hr.daysTrade = putQty .* requiredDelta ./ maxTrade;
  
  
  hr.putLeg = sum(putPrice .* putQty);
  hr.rainbow.qty = min(putQty);
  hr.clientPays = hr.putLeg - (hr.rainbow.qty * hr.rainbow.v(1));
  
  
  

  
end