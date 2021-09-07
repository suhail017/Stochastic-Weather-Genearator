%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction U_parameter       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[a,dDem,rhodDe,sigmadDe] = vap_pre_parameter(ea,esat,Ta,Rsw,Datam)
%%% TEMPERATURE PARAMETER ESTIMATION
%%%Monthly basis
%%%  INPUT %%%
%%% Ta = Temperature
%%%....
%%% Datam %% [Yr, MO, DA, HR]
% %% OUTPUT
%%% a regression coefficients of Deterministic component vapor pressure deficit 
%%% dDem, mean random humidity
%%% rhoDe, lag-1 autocorrelation random humidity
%%% sigmaDe %% standard deviation random humidity
%%% skedDe %% skewness deviation random humidity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Correzione NaN
%%% adj dimension
n=length(ea);
ea=reshape(ea,1,n);
esat=reshape(esat,1,n);
Rsw=reshape(Rsw,1,n);
Ta=reshape(Ta,1,n);
Ta3=sign(Ta).*abs(Ta).^3;  
%Y=esat(Rsw>0)-ea(Rsw>0); %%% [] Vapor pressure deficit 
Y=esat(3:end)-ea(3:end);
%%% Coefficients
%%%% Computation coefficient
t=Datam(:,4); %% actual hour
dt=1;
XX=[ones(1,length(Y)); Ta3(3:end); Rsw(2:end-1); Rsw(1:end-2);];
a=regress(Y',XX');
Des=a'*XX;
Des(Des<0)=0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dDe= Y-Des;
dDe=dDe(not(isnan(dDe)));
dDem=mean(dDe); %mean
%R=autocorr(dDe,10);
%rhodDe=R(2); %lag-1 autocorrelation 
R=xcov(dDe,dDe,10,'coeff');
rhodDe=R(12); %lag-1 autocorrelation 
sigmadDe =std(dDe);%%  standard deviation random
skedDe=skewness(dDe); %%
clear R
end