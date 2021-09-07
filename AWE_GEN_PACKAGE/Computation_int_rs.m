%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunction Computation r and s  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Is,Ir]= Computation_int_rs(t1,t2,dt,delta_S,Lat,T_sunrise,T_sunset,DeltaTSL)
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
T12 = 12 - DeltaTSL;
%%% First evaluation 
TL1  = t1 - DeltaTSL;  %%% time [h]
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (TL1 >= T_sunrise) && (TL1 <= T_sunset) 
s1= sin(delta_S)*sin(phi)-cos(delta_S)*cos(phi)*cos(pi*TL1/12);
else 
    s1=0 ; 
end 
if (TL1 >= T_sunrise) && (TL1 <= T12) 
r1= (pi/12)*cos(delta_S)*cos(phi)*sin(pi*TL1/12);
else 
    r1=0 ; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Second evaluation 
TL2  = t2 - DeltaTSL;  %%% time [h]
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (TL2 >= T_sunrise) && (TL2 <= T_sunset) 
s2= sin(delta_S)*sin(phi)-cos(delta_S)*cos(phi)*cos(pi*TL2/12);
else 
    s2=0 ; 
end 
if (TL2 >= T_sunrise) && (TL2 <= T12) 
r2= (pi/12)*cos(delta_S)*cos(phi)*sin(pi*TL2/12);
else 
    r2=0 ; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Integral first Expression
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (TL1 >= T_sunrise) && (TL2 <= T_sunset)
    [Is] = dt*(s1+s2)/2;
    if TL2 <= T_sunrise
        dt= abs(T_sunrise-TL1);
        [Is] = (s1)*dt/2;
    end
    if TL1 >= T_sunset
        dt= abs(TL1-T_sunset);
        [Is] = (s2)*dt/2;
    end
    if (TL2 <= T_sunrise) && (TL1 >= T_sunset)
        dt= abs(T_sunrise-T_sunset);
        [Is] = sin(delta_S)*sin(phi)-cos(delta_S)*cos(phi)*cos(pi*((TL1+TL2)/2)/12)*(dt/2);
    end
else
    Is=0;
end
%%%%%%%%%%%%%%%%%%%%%
%%% Integral third Expression
if (TL1 >= T_sunrise) && (TL2 <= T12)
    [Ir] = (r1+r2)*dt/2;
    if TL2 <= T_sunrise
         dt= abs(T_sunrise-TL2);
         [Ir] = (r1)*dt/2;
    end
    if TL1 >= T12
         dt= abs(TL2-T12);
         [Ir] = (r2)*dt/2;
    end
    if (TL2 <= T_sunrise) && (TL1 >= T12)
         dt= abs(T_sunrise-T12);
         [Ir] = (pi/12)*cos(delta_S)*cos(phi)*sin(pi*((TL1+TL2)/2)/12)*(dt/2);
    end
else
    Ir=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



