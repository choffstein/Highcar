function Local_vs_Dupire=Local_vs_Dupire()
%  Function Local_vs_Dupire 
%  This function plot the implied volatility and the local volatility (with different kind of methods) 
%  for a set of call price from the real market data with different strikes and maturities.
%
% Inputs = Underlying asset price
%          Exercise price
%          Call options prices for different Maturities
%          Riskfree and Dividend rate
% Outputs= Implied volatility
%          Local Volatility 
%               using Dupire Equation in terms of Strikes and Maturities
%               A) = using Dupire formula
%               B) = using Dupire equation
%               using Dupire equation in terms of imply volatility
%               A) = using Dupire equation
%
% NOTE1: The implied volatility is calculated using Newton-Raphson method and
%        Black-Schooles equation.
% NOTE2: For solving the PDE in the Dupire equation I use the theta method
%        If Theta = 1   is called the explicit method.
%        If Theta = 0.5 is called the Crank Nicholson method.
%
% If you want to see the data type in matlab worksapce "type vol.m"

% Data from the newspapers.
Underlying_asset_price = 34;

Exercise_price = [30 35 40];

Call_price =     [4.99 6.09 7.44;
                  1.58 2.79 4.05;
                  0.22 0.85 1.71]';
                
                    % Put 0 when there is not available the call price
              
Month_of_today = 0;       % month of the year, i.e.  January = 1/12
Maturity_in_years = [1/12, 3/12, 6/12];   % months of the year, i.e. June = 6/12 
Risk_free_interest_rate = 0.045;
Dividend_rate = 0;


% Newton Rapson method variables
Maximum_number_of_iterations = 1000;
Tolerance_of_convergence = 1e-6;

% Theta method variables.
Theta_variable = .5;         % 0 <= Theta <= 1

% NOTES IMPORTANTS
%
% A) Imply volatility
% This function solves for the implied volatility, using Newton's method and Black-Schooles equation.
% MAXITER is the maximum number of iterations used in solving for V.  By
% default, MAXITER = 50.
% By TOL is the tolerance when the result of the Newton-Raphson method 
% is considered to have converged; default is 1e-6.
% Note: This function uses normcdf and normpdf, the normal cumulative distribution 
% and normal probability density functions in the Statistics Toolbox. 
%
% B) Local volatility
% This function solves for the local volatility, using dupire equation.
% For solving the PDE I use the theta method. 
% If Theta = 1   is called the explicit method.
% If Theta = 0.5 is called the Crank Nicholson method.

%VARIABLES IN OUR WORKSPACE

% Define variable names in our workspace
S0 = Underlying_asset_price;
K = Exercise_price;
C = Call_price;
T0 = Month_of_today;
T = Maturity_in_years;
R = Risk_free_interest_rate;
Q = Dividend_rate;
NK = size(K,2);  % = Num_of_strikes 
NT = size(C,1);  % = Num_of_maturities 

% Newton Rapson method variables in our workspace
MAXITER = Maximum_number_of_iterations;
TOL = Tolerance_of_convergence;

% Theta method variables in our workspace 
Theta = Theta_variable;

%CALCULATIONS

% Calculate Implied volatility with MATLAB
for i=1:NT
for j=1:NK
    if C(i,j)==0
        Imp_vol(i,j) = 0;
    else
        
        Imp_vol(i,j)=bsm.ivol(C(i,j), S0, K(1,j), T(i), T0, R, Q);
    end
end
end

Imp_vol

% Calculate de Local Volatility using dupire equation in terms of Strikes and Maturity

% define delta for maturities
for i=1:NT-1
    delta_T(i)=T(i+1)-T(i);
end
% define delta for strikes
for j=1:NK-1
    delta_K(j)=K(j+1)-K(j);
end
% parcial derivative of Call / respect MATURITY
for i=1:NT-1
for j=1:NK
    if C(i+1,j)==0      dp_C_r_T(i,j)=0;
    else    dp_C_r_T(i,j)=(C(i+1,j)-C(i,j))/delta_T(i);
            dp_C_r_T(i+1,j)=0;
    end
end
end
% parcial derivative of Call / respect STRIKES
for i=1:NT-1
for j=2:NK-1   
    if C(i,j-1)==0 | C(i,j+1)==0 | C(i+1,j-1)==0 | C(i+1,j+1)==0
        dp_C_r_K(i,j)=0; 
    else   
dp_C_r_K(i,j)=[Theta*(-C(i,j-1)+C(i,j+1))+(1-Theta)*(-C(i+1,j-1)+C(i+1,j+1))]/(2*delta_K(i)); 
    end
    dp_C_r_K(i+1,j)=0;
    dp_C_r_K(i,j+1)=0;
end
end
% second parcial derivative of Call / respect STRIKES
for i=1:NT-1
for j=2:NK-1 
    if C(i,j)==0 | C(i,j-1)==0 | C(i,j+1)==0 | C(i+1,j)==0 | C(i+1,j-1)==0 | C(i+1,j+1)==0        
        sec_dp_C_r_K(i,j)=0; 
    else
sec_dp_C_r_K(i,j)=[Theta*(C(i,j-1)-2*C(i,j)+C(i,j+1))+(1-Theta)*(C(i+1,j-1)-2*C(i+1,j)+C(i+1,j+1))]/(delta_K(i)*delta_K(i));
    end
    sec_dp_C_r_K(i+1,j)=0;
    sec_dp_C_r_K(i,j+1)=0;
end
end
%Using Dupire Equation in terms of Strikes and Maturity
% Loc_vol1 = dupire formula
% Loc_vol2 = internet equation
% Loc_vol3 = JH equation
for i=1:NT-1
for j=2:NK-1
    if sec_dp_C_r_K(i,j)==0 | dp_C_r_T(i,j)==0 | dp_C_r_K(i,j)==0  
        Loc_vol1(i,j)=0;
        Loc_vol2(i,j)=0;
        Loc_vol3(i,j)=0;
    else
Loc_vol1(i,j)=sqrt(2*dp_C_r_T(i,j)/[K(1,j)*K(1,j)*sec_dp_C_r_K(i,j)]);
Loc_vol2(i,j)=sqrt([2*(dp_C_r_T(i,j)-((R-Q)*(C(i,j)-K(1,j)*dp_C_r_K(i,j))))]/[K(1,j)*K(1,j)*sec_dp_C_r_K(i,j)]);
Loc_vol3(i,j)=sqrt([2*(dp_C_r_T(i,j)+((R-Q)*(K(1,j)*dp_C_r_K(i,j))+Q*C(i,j)))]/[K(1,j)*K(1,j)*sec_dp_C_r_K(i,j)]);
    end
    Loc_vol1(i,j+1)=0;
    Loc_vol1(i+1,j)=0;
    Loc_vol2(i,j+1)=0;
    Loc_vol2(i+1,j)=0;
    Loc_vol3(i,j+1)=0;
    Loc_vol3(i+1,j)=0;
end
end
delta_K;
delta_T;
dp_C_r_T;
dp_C_r_K;
sec_dp_C_r_K;
Loc_vol1;
Loc_vol2;
Loc_vol3;

% Calculate de Local Volatility using dupire equation in terms of imply volatility
% parcial derivative of imply volatility / respect MATURITY
for i=1:NT-1
for j=1:NK
    if Imp_vol(i+1,j)==0      dp_IV_r_T(i,j)=0;
    else    dp_IV_r_T(i,j)=(Imp_vol(i+1,j)-Imp_vol(i,j))/delta_T(i);
            dp_IV_r_T(i+1,j)=0;
    end
end
end
% parcial derivative of imply volatility / respect STRIKES
for i=1:NT-1
for j=2:NK-1   
    if Imp_vol(i,j-1)==0 | Imp_vol(i,j+1)==0 | Imp_vol(i+1,j-1)==0 | Imp_vol(i+1,j+1)==0          
       dp_IV_r_K(i,j)=0;  
    else   
dp_IV_r_K(i,j)=[Theta*(-Imp_vol(i,j-1)+Imp_vol(i,j+1))+(1-Theta)*(-Imp_vol(i+1,j-1)+Imp_vol(i+1,j+1))]/(2*delta_K(i)); 
    end
    dp_IV_r_K(i+1,j)=0;
    dp_IV_r_K(i,j+1)=0;
end
end
% second parcial derivative of imply volatility / respect STRIKES
for i=1:NT-1
for j=2:NK-1 
    if Imp_vol(i,j)==0 | Imp_vol(i,j-1)==0 | Imp_vol(i,j+1)==0 | Imp_vol(i+1,j)==0 | Imp_vol(i+1,j-1)==0 | Imp_vol(i+1,j+1)==0        
        sec_dp_IV_r_K(i,j)=0; 
    else
sec_dp_IV_r_K(i,j)=[Theta*(Imp_vol(i,j-1)-2*Imp_vol(i,j)+Imp_vol(i,j+1))+(1-Theta)*(Imp_vol(i+1,j-1)-2*Imp_vol(i+1,j)+Imp_vol(i+1,j+1))]/(delta_K(i)*delta_K(i));
    end
    sec_dp_IV_r_K(i+1,j)=0;
    sec_dp_IV_r_K(i,j+1)=0;
end
end
%Using Dupire Equation in terms of implied volatility
% Loc_vol_IV using less brackets ();
% Loc_vol_IV_1 test using more brackets ();
for i=1:NT-1
for j=2:NK-1
    if sec_dp_C_r_K(i,j)==0 | dp_C_r_T(i,j)==0 | dp_C_r_K(i,j)==0  
        Loc_vol_IV(i,j)=0;
        Loc_vol_IV_1(i,j)=0;
    else
d1(i,j)=[log(S0/(K(1,j)*1.000001))+(R-Q+0.5*(Imp_vol(i,j)^2)*(T(i)-T0))]/(Imp_vol(i,j)*(sqrt(T(i)-T0)));   
d2(i,j)=[1+K(1,j)*d1(i,j)*dp_IV_r_K(i,j)*sqrt(T(i)-T0)]^2+[Imp_vol(i,j)^2*K(1,j)^2*(T(i)-T0)*(sec_dp_IV_r_K(i,j)-d1(i,j)*dp_IV_r_K(i,j)^2*(sqrt(T(i)-T0)))];
Loc_vol_IV(i,j)=sqrt([Imp_vol(i,j)^2+2*Imp_vol(i,j)*(T(i)-T0)*(dp_IV_r_T(i,j)+(R-Q)*K(1,j)*dp_IV_r_K(i,j))]/d2(i,j));
%test using more brackets ()
d3(i,j)=[(log(S0/K(1,j)))+(R-Q+0.5*(Imp_vol(i,j)^2)*(T(i)-T0))]/(Imp_vol(i,j)*(sqrt(T(i)-T0))); 
d4(i,j)=[(Imp_vol(i,j)^2)*(K(1,j)^2)*(T(i)-T0)*(sec_dp_IV_r_K(i,j)-d3(i,j)*(dp_IV_r_K(i,j)^2)*(sqrt(T(i)-T0)))];
d5(i,j)=([1+K(1,j)*d3(i,j)*dp_IV_r_K(i,j)*(sqrt(T(i)-T0))]^2)+d4(i,j);
Loc_vol_IV_1(i,j)=sqrt([(Imp_vol(i,j)^2)+(2*Imp_vol(i,j)*(T(i)-T0)*(dp_IV_r_T(i,j)+(R-Q)*K(1,j)*dp_IV_r_K(i,j)))]/d5(i,j));
    end
    Loc_vol_IV(i,j+1)=0;
    Loc_vol_IV(i+1,j)=0;
    Loc_vol_IV_1(i,j+1)=0;
    Loc_vol_IV_1(i+1,j)=0;
end
end

dp_IV_r_T;
dp_IV_r_K;
sec_dp_IV_r_K;
Loc_vol_IV
Loc_vol_IV_1;

% Difference betwwen Local and Impled Volatilit
for i=1:NT
for j=2:NK
    if Loc_vol_IV(i,j)==0 
        Dif(i,j)=0;
    else
    Dif(i,j)=Imp_vol(i,j)-Loc_vol_IV(i,j);
    end
end
end

Dif


% Plots
% Local Volatility using Dupire Equation in terms of Strikes and Maturity
%               Loc_vol1 = using Dupire formula
%               Loc_vol3 = using Dupire equation
% Local Volatility using Dupire equation in terms of imply volatility
%               Loc_vol_IV = using Dupire equation
for i=1:NT
    subplot(NT,1,i)
plot(K,Imp_vol(i,:),K,Loc_vol1(i,:),K,Loc_vol3(i,:),K,Loc_vol_IV(i,:));
%axis([K(1,2),K(1,NK-1),0,.4]);
YLABEL(+ Maturity_in_years(i));
GRID ON

if i==1    
TITLE('Implied Volatility and Local Volatility                                               ');
legend('Implied Volatility','Dupire formula','Dupire equation (in term of call prices)','dupire equation (in term of Implied volatility')
end
if i==NT    
XLABEL('Maturity vs Strikes');
end
end


