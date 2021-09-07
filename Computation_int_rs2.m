%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunction Computation r and s  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Is,Ir]= Computation_int_rs2(t1,delta_S,Lat,T_sunrise,T_sunset,DeltaTSL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT
%%% t1 e t2 integral boundary
%%% dt [h] time step
%%delta_S,  Solar declination
phi = Lat*pi/180; %% Latitude [rad]
%%T_sunrise, [h]  sunrise time,
%%T_sunset,  [h]  sunset time,
%%%%%%%%%%%%%%%%%%%%%%%%
%%OUTPUT
%%% Ir integral r t1-t2
%%% Is  integral s t1-t2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T0 = 0 - DeltaTSL;
T23 = 23.0 - DeltaTSL;
%%% First evaluation
TL  = t1 - DeltaTSL;  %%% time [h]
%%%%%%%%%%%%%%%%%%%%%%%%%%%
SunSetHrLoc = T_sunset; 
SunRisHrLoc = T_sunrise;
aa = 0.0005;
del=delta_S;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (DeltaTSL >= 0.0)
    Rho = floor(SunRisHrLoc+1.0) - DeltaTSL;
    if (Rho < SunRisHrLoc)
        Rho=Rho+1;
    end
    Sigma = floor(SunSetHrLoc+1.0) - DeltaTSL;
    if (Sigma < SunSetHrLoc)
        Sigma=Sigma+1;
    end
    T12 = 13 - DeltaTSL;
elseif (DeltaTSL < 0.0)
    Rho = floor(SunRisHrLoc) -  DeltaTSL;
    if (Rho < SunRisHrLoc)
        Rho=Rho+1;
    end
    Sigma = floor(SunSetHrLoc) -  DeltaTSL;
    if (Sigma < SunSetHrLoc)
        Sigma=Sigma+1;
    end
    T12 = 12.0 - DeltaTSL;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%  // Determine the appropriate range for x2 and x3
%%//... before sunrise ...
if (T0 <= TL && TL < SunRisHrLoc)
    x2 = 0.0;
    x3 = 0.0;
    range = 1;
    %% //... sunrise ...
elseif (Rho-aa <= TL && Rho+aa >= TL)
    A = pi*SunRisHrLoc/12.0;
    B = pi*Rho/12.0;
    x2 = (Rho-SunRisHrLoc)*sin(phi)*sin(del);
    x2 = x2- ((12.0/pi)*cos(del)*cos(phi)*(sin(B)-sin(A)));
    x3 = cos(del)*cos(phi)*(cos(A)-cos(B));
    range = 2;
    %% //... morning hours ...
elseif (Rho+aa <= TL && TL <= 12)
    A = pi*TL/12.0;
    B = pi*(TL-1)/12.0;
    x2 = sin(phi)*sin(del);
    x2 = x2-((12.0/pi)*cos(del)*cos(phi)*(sin(A)-sin(B)));
    x3 = cos(del)*cos(phi)*(cos(B)-cos(A));
    range = 3;
    %% ... noon ...
elseif (T12-aa <= TL && T12+aa >= TL)
    A = pi*TL/12.0;
    B = pi*(TL-1)/12.0;
    C = pi*(T12-1)/12.0;
    x2 = sin(phi)*sin(del);
    x2 = x2- ((12.0/pi)*cos(del)*cos(phi)*(sin(A)-sin(B)));
    x3 = cos(del)*cos(phi)*(cos(C)+1);
    range = 4;
    %%//... afternoon hours ...
elseif (T12+aa <= TL && TL < SunSetHrLoc)
    A = pi*TL/12.0;
    B = pi*(TL-1)/12.0;
    x2 = sin(phi)*sin(del);
    x2 = x2- ((12.0/pi)*cos(del)*cos(phi)*(sin(A)-sin(B)));
    x3 = 0.0;
    range = 5;
    %%%  //... sunset ...
elseif (Sigma-aa <= TL && Sigma+aa >= TL)
    A = pi*SunSetHrLoc/12.0;
    B = pi*(Sigma-1)/12.0;
    x2 = (SunSetHrLoc-Sigma+1.0)*sin(phi)*sin(del);
    x2 = x2+ ((12.0/pi)*cos(del)*cos(phi)*(sin(B)-sin(A)));
    x3 = 0.0;
    range = 6;
    %%% //... evening hours ...
elseif (Sigma+aa <= TL && TL <= T23)  %% correction T23 instead 23 
    x2 = 0.0;
    x3 = 0.0;
    range = 7;
end
Is=x2;
Ir=x3; 
return 
  %%%%%%%%%%%
