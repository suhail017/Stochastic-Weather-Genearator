%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Subfunction ComputeCloudFraction      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N,mt] = ComputeCloudFraction(M0,mtm1,sigmam,rhom,t,t0,tb,gam,SS,acloud,bcloud,Ntm1,EN1)
%%%%%%%%%%OUTPUT 
%%%N cloudiness [0-1] 
%%%mt  correlated deviations
%%%%INPUT 
% M0 = mean fair weather cloudiness 
% sigmam = std fair weather cloudiness 
% rhom = lag-1 autocorrelation of fair weather cloudiness 
% chsi = gam = cloudiness decay rate [1/h]
% a b = parame beta distribution cond cloudin previous hour [11,1]
%%% mtm1 = mt(t-1)  
%%i  = t %% time step 1 hour 
%%%Ntm1  = Nt(t-1)  
%choice of the parameters a and b  function of N(t-1)% N(t-1)-->0:1:10 
ind=round(Ntm1*10)+1;
%%%%
chsi=gam;
a = acloud(ind); %% a parameter beta pdf
b = bcloud(ind); %% b parameter beta pdf
if  isequal(SS,1) % Intrastorm period, set N = 1
    N  = 1.0;
    mt = 0;
else %%% Interstorm period
    if isequal(t0,t) %% To avoid problem at interstorm begin
        t0 = t0 - 1;
    end
    %%%% transiction function computation
    J = (1 - exp(-chsi*(t - t0)))*(1 - exp(-gam*(t0 + tb - t)));
    if  J > 0.99
        J = 1.0;
    end
    %%%%%%%% Limitation of the random deviate from  [[[ check ok]]]
    l1   = M0 + (EN1 - M0)*(1 - J);
    l2   = 1 - l1;
    den  = sigmam*J*sqrt(1 - rhom^2);
    L1   = (-l1 - rhom*J*mtm1)/den;
    L2   = (l2 - rhom*J*mtm1)/den;
    %%%%%%%%%%%
    epst = betarnd(a,b); % random deviate beta distribution
    %epst = normrnd(0,1);
    epst = L1 + (L2 - L1)*epst; % random deviate to constrain N between 0-1
    mt   = rhom*mtm1 + epst*sigmam*sqrt(1 - rhom^2); %% correlated deviations (5)
    N    = M0 + (EN1 - M0)*(1 - J) + mt*J; %% cloud cover (1)
    %%%%%%%% 
    % Constrain in case of overshoot (Why ??) 
    if (N>1.0)
        N = 1.0;
    elseif (N<0.0)
        N = 0.0;
    end
        %%%%%%%%%%
    if isnan(N)
        N = Ntm1;
    end
    %%%%%%%%%%%%%
end
return 