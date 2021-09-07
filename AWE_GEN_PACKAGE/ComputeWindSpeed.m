%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction   Compute wind speed       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Ws,dWs]=ComputeWindSpeed(Rswtm3,Rswtm2,Rswtm1,Rsw,c,EdWs,dWstm1,rhodWs,sigmadWs,skedWs)
%%%OTUPUT
%%% Ws wind speed [m/s] (t)
%%%%%%%%%%
%%% INPUT
Wsdet= c(1)+ c(2)*Rsw+c(3)*Rswtm1+c(4)*Rswtm2+c(5)*Rswtm3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psit   = normrnd(0,1); %% standard normal deviate
skfact = 0; %% skewness factor 
while(skfact ~= 1) %%% Constrain psit between [-2.8 e 2.8 ] validity of W-H transforamtion 
    if(((skedWs < 0)&&(psit <= -2.8))||((skedWs > 0)&&(psit >= 2.8)))
        psit =  normrnd(0,1);
        skfact = 0;
    else
        skfact = 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gammae = (1 - rhodWs^3)*skedWs/((1 - rhodWs^2)^(1.5));  %% skewness of  random deviate wind
epst   = (2/gammae)*(1 + (gammae*psit)/6 - gammae^2/36)^3 - (2/gammae);  %% random deviate wind
dWs     = EdWs + rhodWs*(dWstm1 - EdWs) + epst*sigmadWs*sqrt(1 - rhodWs^2);
%%%
Ws = Wsdet + dWs; %% De  %%
if (Ws < 0.0)  %% check negative wind velocity 
    Ws = 0.01;
end
return 