function[DS,WS,Xds,Xws,ttre,pvre]=sample_dry_and_wet_extreme(P,dt,TH)
%%dt =  intervallo di registrazione piogge [h]
%%% h = intervallo sui cui voglio proprietà [h]
%%% Pr  precipitation [mm]
%%% TH threshold on the precipitation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ANSW_daily=1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ANSW_daily == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h=24; %% daily
    fr=h/dt;
    %%%%%%%%%%%%%%%%%%
    Pr=P(find(not(isnan(P))));
    n=length(Pr);
    %%%%%%%%%%%%%%%%%%%%
    m=floor(n/fr);
    Prp=reshape(Pr(1:m*fr),fr,m);
    if fr == 1
        X = Prp;
    else
        X=sum(Prp); %% [mm]
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DRY STATS
    SS = (X > TH); %% TH mm su 1 day
    NT=length(SS);
    [tb,t0] = ConditionSS(SS,NT);
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
                    IS(j)=IS(j)+1; %%% Interstorm period
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
    k=20;
    Xds=Xs(end-k:end);
    ppe=pp(end-k:end);
    ttre=ttr(end-k:end);
    pvre=pvr(end-k:end);
    %%%%%%%%%%%%
    DS=IS;
    clear IS SS NT tb t0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WET STATS
    SS = (X <= TH); %% TH mm su 1 day
    NT=length(SS);
    [tb,t0] = ConditionSS(SS,NT);
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
                    IS(j)=IS(j)+1; %%% Interstorm period
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
    k=20;
    Xws=Xs(end-k:end);
    WS=IS;
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h=1; %% hourly
    fr=h/dt;
    %%%%%%%%%%%%%%%%%%
    Pr=P(find(not(isnan(P))));
    n=length(Pr);
    %%%%%%%%%%%%%%%%%%%%
    m=floor(n/fr);
    Prp=reshape(Pr(1:m*fr),fr,m);
    if fr == 1
        X = Prp;
    else
        X=sum(Prp); %% [mm]
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DRY STATS
    SS = (X > TH); %% TH mm su 1 day
    NT=length(SS);
    [tb,t0] = ConditionSS(SS,NT);
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
                    IS(j)=IS(j)+1; %%% Interstorm period
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
    k=20;
    Xds=Xs(end-k:end);
    ppe=pp(end-k:end);
    ttre=ttr(end-k:end);
    pvre=pvr(end-k:end);
    %%%%%%%%%%%%
    DS=IS;
    clear IS SS NT tb t0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WET STATS
    SS = (X <= TH); %% TH mm su 1 day
    NT=length(SS);
    [tb,t0] = ConditionSS(SS,NT);
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
                    IS(j)=IS(j)+1; %%% Interstorm period
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
    k=20;
    Xws=Xs(end-k:end);
    WS=IS;
end

return