%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [Nfw] = Fair_weather_cloudiness(N,SS,Tr)
%%% CLOUD COVER PARAMETER ESTIMATION
%%%  INPUT %%% N cloudiness
%%% SS Rainfall detector 1 0
%%% Tr = length transiction period [1] [h]
%%% OUTPUT 
%%% Nfw -- fair weather cloudiness 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N(isnan(N))=nanmean(N);
NT=length(SS);
[tb,t0] = ConditionSS(SS,NT); % [h] [h] Computes tb's and t0's during rainstorm events
%%%%%%%%%%%%%%%%%%%%%
i=0;  j=0;
for  i=1:NT
    if t0(i) > 0
        if ((i==1)||((t0(i-1)==0)&&(i>1)))
            j=j+1;  IS{j}=[];
            r=i-1;
            while(((r+1)<=NT) && (t0(r+1)>0))
                r=r+1;
                IS{j}=[IS{j} N(r)]; %%% Interstorm period
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Computation Fair weather series
j=0; Nfw=[];
for j=1:length(IS)
    ISN=IS{j};
    if length(ISN) > 2*Tr
        Nfw=[Nfw, ISN(1+Tr:end-Tr)];
    end
    clear mk ISN
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
%%%%%%%%%%%%%%%%%%%%%%%%