%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction   Compute atmospheric pressure        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function[Pre]=ComputeAtmPressure(EPre,rhopre,Prem1,sigmapre)
%%%OTUPUT
%%% Pre Pressure (t)
%%%%%%%%%%
%%% INPUT
%%% EPre = mean pressure
%%% rhopre =lag-1 autocorrelation pressure
%%% Prem1, =pressure (t-1)
% sigmapre = standard deviation pressure [m/s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epst   =  normrnd(0,1); %% standard normal deviate %% random deviate wind
Pre    = EPre + rhopre*(Prem1 - EPre) + epst*sigmapre*sqrt(1 - rhopre^2);
%%%
end