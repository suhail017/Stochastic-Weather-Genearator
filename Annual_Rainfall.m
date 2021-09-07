%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pyr,MEr] = Annual_Rainfall(PE,PCv,Prho,Pskw,NYR,err)
%%% PARAMETER ESTIMATION
%%%  INPUT %%%
%PE = Mean of annual precipitation process [mm]
%PCv  = Coefficient of Variation of precipitation process [-]
%Prho = Lag-1 Autocorrelation of precipitation process [-]
%Pskw = Skewness of Precipitation process
%NYR = Number of year of simulation
%err in annual precipitation measurements [-]
%%%OUTPUT
% Pyr = Time series of annual precipitation
% MEr = Acceptance of error
if nargin <6
    err = 0.025;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MEr= PE*err;
Psigma = PE*PCv; %%[mm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
rn = zeros(1,length(NYR));
for i=1:NYR
    rn(i)=normrnd(0,1);
    skfact = 0; %% skewness factor
    while(skfact ~= 1) %%% Constrain rn between [-2.8 e 2.8 ] validity of W-H transforamtion ???
        if (((Pskw < 0) && (rn(i) <= -2.8)) || ((Pskw > 0)&&(rn(i) >= 2.8)))
            rn(i) =  normrnd(0,1);
            skfact = 0;
        else
            skfact = 1;
        end
    end
end
%%%%%%%%%%%%%%%%
gam = (1 - Prho^3)*Pskw/((1 - Prho^2)^(1.5));  %% gamma variate
er  = (2/gam)*(1 + (gam*rn)/6 - (gam^2)/36).^3 - (2/gam);  %% random deviate skewed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xs=zeros(1,length(NYR));
Xs(1)=PE;
for i=2:NYR
    Xs(i)  = PE + Prho*(Xs(i-1)-PE) + er(i)*Psigma*sqrt(1 - Prho^2);
end
Xs(Xs<=0)=NaN;
Xs(isnan(Xs))=nanmin(Xs);
Pyr=Xs;
%%%%%%%%%%%%%%%%%%%%%%%%
end 
