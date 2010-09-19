function p = getParam()

  p.nSim = 10000;
  p.profit = 0.10;
  p.notional = 100000000;
  p.K_put = 0.60;
  p.K_call = 1.1;
  p.rho = 0.3;
  p.r = 0.045;
  p.T = 2;
  
  p.S(1).name = 'Green Co.';
  p.S(1).sharesHeld = 3e6;
  p.S(1).x_0 = 34;
  p.S(1).adv = 217000;
  p.S(1).T = [1/12; 3/12; 6/12];
  p.S(1).K = [30; 35; 40];
  p.S(1).c = [4.99 6.09 7.44; ...
              1.58 2.79 4.05; ...
              0.22 0.85 1.72]';
  p.S(1).fixedSigma = [];        % change to value to use fixed local value         
            
  p.S(2).name = 'SynerOptics';
  p.S(2).sharesHeld = 7e6;
  p.S(2).x_0 = 15;
  p.S(2).adv = 610000;
  p.S(2).T = [1/12; 3/12; 6/12];
  p.S(2).K = [10; 15; 20;];
  p.S(2).c = [5.10 5.31 5.66; ...
              1.15 1.87 2.47; ...
              0.06 0.39 0.91]';
  p.S(2).fixedSigma = []; 
            
  p.S(3).name = 'SW Industries';
  p.S(3).sharesHeld = 4e6;
  p.S(3).x_0 = 25;
  p.S(3).adv = 802000;
  p.S(3).T = [1/12; 3/12; 6/12];
  p.S(3).K = [20; 25; 30];
  p.S(3).c = [5.35 6.03 6.89; ...
              1.60 2.53 3.44; ...
              0.16 0.63 1.31]';
  p.S(3).fixedSigma = []; 
end