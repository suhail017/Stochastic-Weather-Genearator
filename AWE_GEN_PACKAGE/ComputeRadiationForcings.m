%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute radiation forcings %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[EB,ED,EB1,EB2,ED1,ED2,PARB,PARD]=ComputeRadiationForcings(Datam,DeltaGMT,Lon,Lat,...
    LWP0,N,Zbas,Tdew,beta_A,alpha_A,omega_A1,omega_A2,uo,un,rho_g)
%%%OTUPUT
%%%  EB %% sky beam irradiance without terrain effects
%%%% ED %% sky total diffuse irradiance  without terrain effects
%%EB1 [W/m^2] sky beam irradiance VIS band[0.29 um - 0.70um ]
%%EB2  [W/m^2]   sky beam irradiance NIR band [0.70 um - 4.0 um ]
%%ED1 [W/m^2]    sky total diffuse flux at the ground  VIS band [0.29 um -0.70um ]
%%%ED2   [W/m^2]  sky  total diffuse flux at the ground NIR band  [0.70 um-4.0 um ]
%%%%%%%%%%
%%% INPUT
%%% Datam %% [Yr, MO, DA, HR]
%%% DeltaGMT [°]
%%% Lon [°]
%%% Lat [°]
%%% LWP0 [g/m^2]
%%% N [0-1] cloudiness
%Zbas [m a.s.l.]  watershed elevation
%Tdew [°C] dew point temperature
%beta_A   0.05 [0-1] Angstrom turbidity parameters
%alpha_A  1.3 [0.5-2.5]Angstrom turbidity parameters 
%omega_A1  %% 0.92 [0.74-0.94]aerosol single-scattering albedo band 1 
%omega_A2 %% 0.84 [0.74-0.94] aerosol single-scattering albedo band 2
%%% uo 0.35 [0.22-0.35]  %% [cm] ozone amount in vertical column
%%% un 0.0002 [0.0001-0.046] %% [cm] total nitrogen dioxide amount 
%%%  rho_g 0.15 [0.05-0.3] %% [] spatial average regional albedo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
So    = 1366.1;  % [W/m^2] Solar constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[delta_S,h_S,zeta_S,T_sunrise,T_sunset,L_day,E0,jDay,Delta_TSL] = SetSunVariables(Datam,DeltaGMT,Lon,Lat); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction SetSunVariables   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sop = So*E0; %% Actual Solar constant [W/m^2]
%%%
%%% Partition Energy two bands Gueymard (2004; 2008)
So1 = Sop*0.4651;  %%  [W/m^2]Extraterrestrial  Radiation VIS band [0.29 um - 0.70 um ]
So2 = Sop*0.5195; %% [W/m^2] %%Extraterrestrial Radiation NIR band [0.70 um - 4.0 um ]
%LWP = exp(N*log(LWP0)) - 1; %%% Liquid Water Path whit cloudiness  [g/m^2]
LWP = LWP0*N; %% Liquid Water Path whit cloudiness  [g/m^2]
if(LWP < 1.1)
    Ntmp = 0; %% Cloudiness
else
    Ntmp = N; % Cloudiness
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(Ntmp,0) % Completely clear sky
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     Subfunction SetClearSkyRadiation      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [EB1,EB2,ED1,ED2,Edp1,Edp2,rho_s1,rho_s2,Mb,Mg] = SetClearSkyRadiation(h_S,Zbas,Tdew,So1,So2,beta_A,alpha_A,omega_A1,omega_A2,uo,un,rho_g);
    %%%%%%%%%%%%%%%%%%%%%%
elseif Ntmp > 0  %% Overcast Condition
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     Subfunction SetClearSkyRadiation      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Eb1,Eb2,Ed1,Ed2,Edp1,Edp2,rho_s1,rho_s2,Mb,Mg] = SetClearSkyRadiation(h_S,Zbas,Tdew,So1,So2,beta_A,alpha_A,omega_A1,omega_A2,uo,un,rho_g);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    Subfunction SetCloudySkyRadiation      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [EB1,EB2,ED1,ED2] = SetCloudySkyRadiation(LWP,h_S,Eb1,Eb2,Edp1,Edp2,N,rho_s1,rho_s2,rho_g); %
    %%%%%%%%%%%%%%%%%%%%%%%
end
    %%%%%%%%%%Correction Flat Surface 
    EB1=EB1*sin(h_S); 
    EB2=EB2*sin(h_S); 
    EB = EB1 + EB2;%% [W/m^2] % beam irradiance
    ED = ED1 + ED2; %% [W/m^2] %total diffuse flux at the ground
    %%%% PAR Radiation estimation 
    PARB = EB1*Mb; 
    PARD = Mg*(EB1+ED1)- PARB; 
    %%%%%%%%%%%%%%%%%%%%%%%%%
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55555
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%