%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Subfunction ComputeAirTemp          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Ta,T_tilde,dT,qt,I2,I3,I4]=ComputeAirTemperature(Datam,DeltaGMT,...
    Lon,Lat,N,b,qtm1,dTbar,rhodT,sigmadT,dTtm1,T_tildetm1,I2tm1,I3tm1,I4tm1,Tave)
%%%OTUPUT
%%% Ta  %% temperature (t) (°C)
%%% qt   = e_cs*sigma*Kc*(Ta(i,1) + 273.15)^4; %%  q(t) estimate incoming long wave radiation
%%% dT; %%% mean random temperature deviate (t)
%%% T_tilde  %% Deterministic component of temperature [°C]
%%% I2,I3,I4  integral components
%%%%%%%%%%
%%% INPUT
%%%%% Datam %% [Yr, MO, DA, HR]
%%% DeltaGMT [°]
%%% Lon [°]
%%% Lat [°]
%%% N cloudiness
%%% b regression coefficients of Deterministic component of temperature first order differential equation
%%% qtm1, q(t-1) estimate incoming long wave radiation
%%% dTbar, mean random temperature deviate
%%% rhodT, lag-1 autocorrelation random temperature deviate
%%% sigmadT %% standard deviation random temperature deviate
%%% skedT %% skewness deviation random temperature deviate
%%% dTtm1, mean random temperature deviate (t-1)
%%% T_tildetm1  %% Deterministic component of temperature [°C] (t-1)
%%% I2tm1,I3tm1,I4tm1  integral components (t-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[delta_S,h_S,zeta_S,T_sunrise,T_sunset,L_day,E0,jDay,Delta_TSL] = SetSunVariables(Datam,DeltaGMT,Lon,Lat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction SetSunVariables   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nowHR =Datam(4);  
%%%%% Compute Deterministic component of air temperature
if nowHR == 0
    I2tm1 = 0;
    I3tm1 = 0;
    I4tm1 = 0;
end
Ti = T_tildetm1;
[T_tilde,I2,I3,I4] = ComputeDeterministicT(Ti,delta_S,Lat,nowHR,T_sunrise,T_sunset,N,b,I2tm1,I3tm1,I4tm1,qtm1,Delta_TSL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Subfunction ComputeDeterministicT      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Compute random deviate of air temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Subfunction ComputeStochasticT        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ta=Inf; kjx=0; 
while ((Ta-Tave) > 25) && (kjx < 100)
    dT = ComputeStochasticT(dTbar,rhodT,dTtm1,sigmadT);%
    Ta = T_tilde + dT; %% Ta  %% temperature (t) (°C)
    kjx=kjx+1; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kc      = 1 + 0.17*N^2; %% cloud effect
sigmaSB = 5.6704e-8; % Stefan-Boltzmann constant [W/m^2 K] 
qt    = sigmaSB*Kc*(Ta + 273.15)^4; %%  q(t) estimate incoming long wave radiation [W/m^2] 
return

