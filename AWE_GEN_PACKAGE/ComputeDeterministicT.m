%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Subfunction ComputeDeterministicT      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T_tilde,I2,I3,I4] = ComputeDeterministicT(Ti,delta_S,Lat,t,T_sunrise,T_sunset,N,b,I2tm1,I3tm1,I4tm1,qtm1,DeltaTSL)
%%%OTUPUT
%%% T_tilde = Deterministic component of temperature [°C] 
%I2
%I3
%I4
%%%%%%%%%%
%%% INPUT
%Ti,  %% initial temperature [°C]
%%delta_S,  Solar declination
phi=Lat; %% Latitude [rad]
%t  [h] time 
%%T_sunrise, [h]  sunrise time, 
%%T_sunset,  [h]  sunset time, 
%N  cloudiness [0-1] 
%b  regression coefficients of Deterministic component of temperature first order differential equation
%I2tm1, I2(t-1) 
%I3tm1, I3(t-1) 
%I4tm1, I4(t-1)
% qtm1, q(t-1) estimate incoming long wave radiation 
%%DeltaTSL [h] Time difference between standard and local meridian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Solution of the first order differential equation, the solution to which can be found if
%%the initial condition, i.e., the initial temperature, Ti, is given. Curtis and Eagleson (1982)
%%%%%%%%% Definition of coefficients
b0    = b(1);
b1    = b(2);
b2    = b(3);
b3    = b(4);
b4    = b(5);
%K     = 1 - 0.65*N^2;
K=  1 - 0.75*N.^(3.4);
phi = phi*pi/180; %% Latitude [rad]
TL  = t - DeltaTSL;  %%% time [h]
TP  = -DeltaTSL - 1;  %% initial time [h]
T12 = 12 - DeltaTSL;
%%% Coefficients
p  = pi/12;
K1 = b0/b1;
K2 = (b2/b1)*sin(delta_S)*sin(phi);
K3 = ((b1*b2)/(b1^2 + p^2))*cos(delta_S)*cos(phi);
K4 = ((p*b2)/(b1^2 + p^2))*cos(delta_S)*cos(phi);
K5 = ((p^2*b3)/(b1^2 + p^2))*cos(delta_S)*cos(phi);
K6 = ((p*b1*b3)/(b1^2 + p^2))*cos(delta_S)*cos(phi);
%%%%%%%%%%%%%%%%%%%%
I1 = K1*(exp(b1*TL) - exp(b1*TP)); %%% Integral First Expression 
I4 = (b4/b1)*(qtm1/1000)*(1 - exp(-b1))*exp(b1*TL) + I4tm1; % Integral Fourth Expression 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Integral second Expression
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (TL >= T_sunrise) && (TL <= T_sunset)
       [I2] = integral2(TL,TL-1,b1,K,K2,K3,K4,I2tm1);
    if TL-1 <= T_sunrise
        [I2] = integral2(TL,T_sunrise,b1,K,K2,K3,K4,I2tm1);
    end
    if TL+1 >= T_sunset
        [I2] = integral2(T_sunset,TL-1,b1,K,K2,K3,K4,I2tm1);
    end
    if (TL-1 <= T_sunrise) && (TL+1 >= T_sunset)
        [I2] = integral2(T_sunset,T_sunrise,b1,K,K2,K3,K4,I2tm1);
    end
else
    I2=I2tm1;
end
%%%%%%%%%%%%%%%%%%%%%
%%% Integral third Expression
if (TL >= T_sunrise) && (TL <= T12)
        [I3] = integral3(TL,TL-1,b1,K,K5,K6,I3tm1);
    if TL-1 <= T_sunrise
        [I3] = integral3(TL,T_sunrise,b1,K,K5,K6,I3tm1);
    end
    if TL+1 >= T12
        [I3] = integral3(T12,TL-1,b1,K,K5,K6,I3tm1);
    end
    if (TL-1 <= T_sunrise) && (TL+1 >= T12)
        [I3] = integral3(T12,T_sunrise,b1,K,K5,K6,I3tm1);
    end
else
    I3=I3tm1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = I1 + I2 + I3 + I4;  %%% Auxiliary variable for the solution
T_tilde = Ti*exp(-b1*(TL - TP)) + G*exp(-b1*TL); %%%   Deterministic component of temperature [°C]
end
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Subfunctions 
function [I2] = integral2(t,t2,b1,K,K2,K3,K4,I2tm1)
I2 = K*(K2*(exp(b1*t)-exp(b1*(t2)))-K3*exp(b1*t)*cos(pi*t/12)-K4*exp(b1*t)*sin(pi*t/12)+...
    K3*exp(b1*(t2))*cos(pi*(t2)/12)+K4*exp(b1*(t2))*sin(pi*(t2)/12)) + I2tm1; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I3] = integral3(t,t2,b1,K,K5,K6,I3tm1)
I3 = K*(K6*exp(b1*t)*sin(pi*t/12)-K5*exp(b1*t)*cos(pi*t/12)-...
    K6*exp(b1*(t2))*sin(pi*(t2)/12)+K5*exp(b1*(t2))*cos(pi*(t2)/12)) + I3tm1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

