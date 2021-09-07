%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunction Computation r and s  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[r,s]= Computation_rs(delta_S,phi,t,T_sunrise,T_sunset,DeltaTSL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TL  = t - DeltaTSL;  %%% time [h]
T12 = 12 - DeltaTSL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (TL >= T_sunrise) && (TL <= T_sunset) 
s= sin(delta_S)*sin(phi)-cos(delta_S)*cos(phi)*cos(pi*TL/12);
else 
    s=0 ; 
end 
if (TL >= T_sunrise) && (TL <= T12) 
r= (pi/12)*cos(delta_S)*cos(phi)*sin(pi*TL/12);
else 
    r=0 ; 
end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%