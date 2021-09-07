function[Dfr,Dwa,Nfr,Nwa,HWe,CWe,ttre,pvre]=sample_temp_extreme(T,dt)
%%dt =  time step registration [h]
% T variable [] %% additional or mean 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%% Dfr  -- Tmax < 0  Icing Days 
%%% Dwa  -- Tmax > 25  Summer Days 
%%% Nfr -- Tmin < 0  Frost Day 
%%% Nwa -- Tmin > 20 Tropical Night 
%%% Hwe -- Consec. Days  T > T90 
%%% Cwe -- Consec. Days  T < T10 
%%% ttre 
%%% pvre 
%%%%%%%%%%%%%%%%%%%%%
h=24; %% Daily 
fr=h/dt;
%%%%%%%%%%%%%%%%%%
n=length(T);
%%%%%%%%%%%%%%%%%%%%
m=floor(n/fr);
Trp=reshape(T(1:m*fr),fr,m);
if fr == 1
    Trp=Trp(find(not(isnan(Trp))));
    X = Trp;
    Xmax=X; 
    Xmin=X; 
else
    X=nanmean(Trp);
    Xmax= nanmax(Trp);
    Xmin= nanmin(Trp); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dfr= sum(Xmax<=0)/m; 
Dwa= sum(Xmax>=25)/m;
Nfr= sum(Xmin<=0)/m; 
Nwa= sum(Xmin>=20)/m;
T10=prctile(X,10);
T90=prctile(X,90);
HW= (X<=T90);
CW= (X>=T10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NT=m; 
[thw,t0] = ConditionSS(HW,NT); 
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
i=0;  j=0;
for i=1:NT
    if t0(i) > 0 
        if ((i==1)||((t0(i-1)==0)&&(i>1)))
            j=j+1;  IS(j)=0;  
            r=i-1;
            while(((r+1)<=NT) && (t0(r+1)>0))
                r=r+1; 
                IS(j)=IS(j)+1; %%%
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xs=sort(IS);
pp=(1:m)/(m+1);
ttr=(1./(1-pp))*h/8766; % ttr = tempo di ritorno legato ad ogni probabilità espresso in anni [yr] 
pvr=-log(log(1./(1-h./(ttr.*8766)))); % pvr = variabile ridotta legata ad ogni probabilità
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% estrai gli ultimi k valori 
k=25; 
HWe=Xs(end-k:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear IS Xs pp ttr pvr t0  ppe 
[tcw,t0] = ConditionSS(CW,NT); 
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
i=0;  j=0;
for i=1:NT
    if t0(i) > 0 
        if ((i==1)||((t0(i-1)==0)&&(i>1)))
            j=j+1;  IS(j)=0;  
            r=i-1;
            while(((r+1)<=NT) && (t0(r+1)>0))
                r=r+1; 
                IS(j)=IS(j)+1; %%
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xs=sort(IS);
pp=(1:m)/(m+1);
ttr=(1./(1-pp))*h/8766; % ttr = tempo di ritorno legato ad ogni probabilità espresso in anni [yr] 
pvr=-log(log(1./(1-h./(ttr.*8766)))); % pvr = variabile ridotta legata ad ogni probabilità
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CWe=Xs(end-k:end);
ttre=ttr(end-k:end); 
pvre=pvr(end-k:end); 
return 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
