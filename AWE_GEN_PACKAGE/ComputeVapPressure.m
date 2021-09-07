%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Subfunction ComputeHumidity           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ea,dDe] = ComputeVapPressure(esat,Ta,Rswtm1,Rswtm2,dDetm1,a,dDem,rhodDe,sigmadDe)
%%% OUTPUT 
%%% U 
%%%  INPUT %%%
%%% Ta = Temperature
%%%.Pr = precipitation 
%%% a regression coefficients of Deterministic component vapor pressure deficit 
%%% dDem, mean random humidity
%%% rhoDe, lag-1 autocorrelation random humidity
%%% sigmaDe %% standard deviation random humidity
%%% skedDe %% skewness deviation random humidity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dedet= a(1)+ a(2)*(Ta^3)+ a(3)*Rswtm1 + a(4)*Rswtm2;%%  
epsDe= normrnd(0,1);
%rn=normrnd(0,1);
%skfact = 0; %% skewness factor
%while(skfact ~= 1) %%% Constrain rn between [-2.8 e 2.8 ] validity of W-H transforamtion ???
 %   if(((skedDe < 0)&&(rn <= -2.8))||((skedDe > 0)&&(rn >= 2.8)))
  %      rn =  normrnd(0,1);
  %      skfact = 0;
  % else
  %      skfact = 1;
  %  end
%end
%gam = (1 - rhodDe^3)*skedDe/((1 - rhodDe^2)^(1.5));  %% gamma variate 
%epsDe = (2/gam)*(1 + (gam*rn)/6 - (gam^2)/36).^3 - (2/gam);  %% random deviate skewed 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dDe= dDem + rhodDe*(dDetm1 - dDem) + epsDe*sigmadDe*sqrt(1 - rhodDe^2);
De = Dedet + dDe; %% De  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ea= esat -De; 
if ea < 0
    ea=1;
end
if ea > esat
    ea=esat;
end
return
%%%%%%%%%%%%%%%%%%%%5