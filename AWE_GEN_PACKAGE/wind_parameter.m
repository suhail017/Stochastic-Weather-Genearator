%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Ws_parameter       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[c,dWsm,rhodWs,sigmadWs,skedWs] = wind_parameter(Ws,Rsw)
%%% WIND SPEED ESTIMATION
%%%Monthly basis
%%%  INPUT %%%
%%% Ta = Temperature
%%%....
%%% Datam %% [Yr, MO, DA, HR]
%%% DeltaGMT [°]
%%% Lon [°]
%%% Lat [°]
% %% OUTPUT
%%% a regression coefficients of Deterministic component of wind speed 
%%% dWsm, mean random humidity
%%% rhoWs, lag-1 autocorrelation random 
%%% sigmaWs %% standard deviation random
%%% skedWs %% skewness deviation random
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Correction NaN
%%% adj dimension
%%% adj dimension
n=length(Ws);
Rsw=reshape(Rsw,1,n);
Ws=reshape(Ws,1,n);
Y=Ws(4:end); %%% [] Wind speed 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Coefficients
%%%% Computation coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XX=[ones(1,n-3);Rsw(4:end) ; Rsw(3:end-1); Rsw(2:end-2); Rsw(1:end-3)];
c=regress(Y',XX');
%%%%%%%%%%%5
%Wsi=nanmean(Ws);
Wss=c'*XX;
Wss(Wss<0.0)=0.0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dWs= Y-Wss;
dWs=dWs(not(isnan(dWs)));
dWsm=mean(dWs); %mean random
%R=autocorr(dWs,10);
%rhodWs=R(2); %lag-1 autocorrelation random 
R=xcov(dWs,dWs,10,'coeff');
rhodWs=R(12); %lag-1 autocorrelation random 
sigmadWs =std(dWs);%%  standard deviation random 
skedWs=skewness(dWs); %%
skedWs(skedWs>1.5)=1.5; 
clear R
end