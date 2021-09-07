%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Subfunction ConditionSS           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[tb,t0] = ConditionSS(SS,NT)
%%%%%%%%%%%
%%% tb time between storm [h]
%%% t0 time begin of storm [h] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%SS rainy indicator 1 rain 0 no rian 
%%%NT length period 
    tb = zeros(NT,1); %% length time between storom 
    t0 = zeros(NT,1); %% storm begin 
    for i=1:NT
        if(SS(i)==0)
            if ((i==1)||((SS(i-1)==1)&&(i>1)))
                t0tmp = i;
                tbtmp = 1;
                j=i;
                while(((j+1)<=NT)&&(SS(j+1)~=1))
                    j     = j + 1;
                    tbtmp = tbtmp + 1;
                end
            end
            t0(i) = t0tmp;
            tb(i) = tbtmp;
        else
            t0(i) = 0;
            tb(i) = 0;
        end
    end
return 