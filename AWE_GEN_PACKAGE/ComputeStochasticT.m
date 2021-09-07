%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Subfunction ComputeStochasticT        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dT = ComputeStochasticT(dTbar,rhodT,dTtm1,sigmadT)
%%%OTUPUT
%%% dT = Random deviate component of temperature [°C] 
%%%%%%%%%%
%%% INPUT
%%% dTbar, mean random temperature deviate 
%%% rhodT, lag-1 autocorrelation random temperature deviate 
%%% dTtm1, mean random temperature deviate (t-1) 
%%% sigmadT %% standard deviation random temperature deviate 
epsT = normrnd(0,1);   %% normal distributed
dT = dTbar + rhodT*(dTtm1 - dTbar) + epsT*sigmadT*sqrt(1 - rhodT^2);
%%%%%%%%%%%%%%%%%%%%%
return 
