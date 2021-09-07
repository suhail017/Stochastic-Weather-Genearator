%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% WEATHER GENERATOR HEART     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [Psim,Pyr_AR_orig,Ns,Tas,eas,Us,Tdews,SB,...
    SD,SAB1,SAB2,SAD1,SAD2,PARB,PARD,Wss,Pres]...
    = WEATHER_GENERATOR(DS,Par,NYR,Pr_yr,phat,Ta_yr,DeltaGMT,Lon,Lat,Zbas,IAV,VERB,DT_fut)
%%%% INPUT 
%%% Par DS NYR Pr_yr Ta_yr DeltaGMT,Lon,Lat,Zbas
%%%     J  F  M  A  M  J  J  A  S  O  N  D
Days = [31 28 31 30 31 30 31 31 30 31 30 31];
DS=reshape(DS,length(DS),1);
[Yr,Mo,Da,Hr,Mi,Se]=datevec(DS);
Datam=[Yr,Mo,Da,Hr]; %% Datam %% [Yr, MO, DA, HR]
%%%%%%%%%%%%%%%%
Par.acloud=reshape(Par.acloud,11,12); %%%%[11x12 double]
Par.bcloud=reshape(Par.bcloud,11,12); %%% %%%[11x12 double]
Par.bTemp = reshape(Par.bTemp,12,5);  %%% [12x5 double]
Par.aVap= reshape(Par.aVap,12,4);  %%%[12x4 double]
%%%%%%%%%%%%%%
%%%% OUTPUT 
%%% PH Pyrs Psim Pyr_AR_orig 
%%% Ns Tas eas Us Tdews SB SD SAB1 SAB2 SAD1 SAD2 Wss Pres 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    RAINFALL GENERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[Pyr_AR_orig,MEr] = Annual_Rainfall(Pr_yr,NYR);
[Pyr_AR_orig,MEr] = Annual_Rainfall(Pr_yr.E,Pr_yr.Cv,Pr_yr.rho,Pr_yr.skw,NYR,phat); 
Pyr_AR=Pyr_AR_orig;
%%%%%%%%%%%%%% Constrain on mean annual precipitation otherwise ARMA may
%%%%%%%%%%%%%% produce sometimes quite different values by chance 
while abs(mean(Pyr_AR_orig)-Pr_yr.E)> 3 
    %[Pyr_AR_orig,MEr] = Annual_Rainfall(Pr_yr,NYR);
    [Pyr_AR_orig,MEr] = Annual_Rainfall(Pr_yr.E,Pr_yr.Cv,Pr_yr.rho,Pr_yr.skw,NYR,phat); 
    Pyr_AR=Pyr_AR_orig;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mb =  helpdlg('RAINFALL GENERATION','PARAMETER_ESTIMATION');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lap=0; PHH=cell(12,NYR);
if VERB == 1
    bau = waitbar(0,'Toward conclusion');
    bau2 = waitbar(0,'Internal Loop');
end
%%%%%%%%%%%%%%%%
while not(nansum(Pyr_AR)==0)
    adv= NYR-nansum(Pyr_AR>0);
    lap=lap+1;
    if VERB == 1
        waitbar(adv/NYR,bau);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NYRT=50; %%% number of year target 
    Ns=250; %% number of simulated storm
    if lap*NYRT >= NYR*10
        NYRT = 10; 
        Ns= 40; 
    end 
    SPD=zeros(1,NYRT); %% Total annual
    for ii=1:12
        %%%%%%%%
        if VERB == 1
            waitbar(ii/12,bau2);
        end
        [PHy]=ComputeRainfall(Ns,Par.lan(ii),Par.bet(ii),Par.muc(ii),Par.eta(ii),Par.alp(ii),Par.tet(ii));
        while length(PHy) < (Days(ii)*24*NYRT)
            Ns=100;
            [PHy2]=ComputeRainfall(Ns,Par.lan(ii),Par.bet(ii),Par.muc(ii),Par.eta(ii),Par.alp(ii),Par.tet(ii));
            PHy=[PHy PHy2];
        end
        PHy(PHy<0.1)=0.0;
        %%%%%%%%%%%%%%%%%%
        k=1;
        for j=1:NYRT
            PHt{j,ii}=[PHy(k:(k-1)+Days(ii)*24)];
            k= k+Days(ii)*24;
            SPD(j) = SPD(j) + sum(PHt{j,ii});
        end
        clear PHy PHy2 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if VERB == 1
        close(bau2);
        bau2 = waitbar(0,'Internal Loop');
    end
    %%%%%%%%%%%%%%%%%%%%%
    sinked = 0; rrr=0; 
    for j=1:NYRT
        [V,Pos]=nanmin(abs(Pyr_AR-SPD(j)));
        %%%%% Possibility of eliminate the inter-annual variability feature
        if (V < MEr) || not(strcmp(IAV,'Yes')) 
            %%%%
            if not(strcmp(IAV,'Yes'))
                rrr=rrr+1;
                Pos = rrr; 
            end
            %%%%
            for ii=1:12;
                PHH{ii,Pos}= PHt{j,ii};
            end
            Pyr_AR(Pos)=NaN;
            sinked = sinked +1;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    if (sinked == 0)
        [SPD_max,Pos1]=max(SPD);
        [SPD_min,Pos2]=min(SPD);
        PHt_max= cell(12,1); PHt_min=cell(12,1);
        for ii=1:12
            PHt_max{ii}=PHt{Pos1,ii};
            PHt_min{ii}=PHt{Pos2,ii};
        end
        [Vmax,Pos_max]=nanmin(abs(Pyr_AR-SPD_max));
        [Vmin,Pos_min]=nanmin(abs(Pyr_AR-SPD_min));
        if Vmax < Vmin
            PHt=PHt_max; Pos=Pos_max;
            rap= Pyr_AR(Pos)/SPD_max;
        else
            PHt=PHt_min; Pos=Pos_min;
            rap= Pyr_AR(Pos)/SPD_min;
        end
        for ii=1:12;
            PHH{ii,Pos}= PHt{ii}*rap;
        end
        Pyr_AR(Pos)=NaN;
        clear SPD_max SPD_min Pos1 Pos2
        clear PHt_min PHt_max Pos_max Pos_min Vmin Vmax
        clear PHt V Pos rap
    end
    clear V Pos PHt 
end
%%%%%%%%%%%%%
if VERB == 1
    close(bau);
    close(bau2);
end
%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PH=cell(12,1);  ch=0;
Psmon=zeros(12,NYR); Psim=[];
for yr=1:NYR;
    for ii=1:12
        ch= ch+Days(ii)*24;
        %if ii==2 && day(DS(ch+1))==29
        if ii==2 && Da(ch+1)==29
            ch=ch+24;
            Ns=20; %% number of simulated storm
            [PHy]=ComputeRainfall(Ns,Par.lan(ii),Par.bet(ii),Par.muc(ii),Par.eta(ii),Par.alp(ii),Par.tet(ii));
            while length(PHy) < 1200
                [PHy]=ComputeRainfall(Ns,Par.lan(ii),Par.bet(ii),Par.muc(ii),Par.eta(ii),Par.alp(ii),Par.tet(ii));
                Ns=Ns+10;
            end
            PHy(PHy<0.1)=0.0;
            PHH{ii,yr}=[PHH{ii,yr}, PHy(1001:1000+24)];
            clear PHy Ns
        end
        PH{ii}=[PH{ii} PHH{ii,yr}];
        Psmon(ii,yr)=sum(PHH{ii,yr});
        Psim=[Psim PHH{ii,yr}];
    end
end
Pyrs=sum(Psmon);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close(mb);
NT=length(DS); 
i=0;  SS=zeros(1,NT); Ns=zeros(1,NT);
Tas=zeros(1,NT); dT=zeros(1,NT);T_tilde=zeros(1,NT);
eas=zeros(1,NT);dDe=zeros(1,NT);
esat_s=zeros(1,NT);Tdews=zeros(1,NT);Us=zeros(1,NT);
SB=zeros(1,NT); SD=zeros(1,NT);SAB1=zeros(1,NT);
SAB2=zeros(1,NT);SAD1=zeros(1,NT);SAD2=zeros(1,NT);
PARD=zeros(1,NT);PARB=zeros(1,NT);
Rsws=zeros(1,NT);
Wss=zeros(1,NT);dWs=zeros(1,NT);Pres=zeros(1,NT);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bau = waitbar(0,'WEATHER GENERATOR COMPUTATION');
SS = Psim>0;    
a_a=17.27; b_b=237.7;
[tb,t0] = ConditionSS(SS,length(SS)); % [h] Computes tb's and t0's during rainstorm events
for i=1:NT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    waitbar(i/NT,bau);  
    %%%%%%%%%%%%%%%%%%%%
    %j=month(DS(i)); 
    %Hr = hour(DS(i))+1; 
    j=Mo(i); 
    Hor = Hr(i)+1; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%    CLOUD COVER  GENERATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i == 1 
         mt = 0; Ns(1)=0;
    else
        mtm1 = mt;
        [Ns(i),mt] = ComputeCloudFraction(Par.M0(j),mtm1,Par.sigmam(j),Par.rhom(j),...
            i,t0(i),tb(i),Par.gam(j),SS(i),Par.acloud(:,j),Par.bcloud(:,j),Ns(i-1),Par.EN1(j));
        %%%%%%%%%%%%%%%%%%%%
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%   TEMPERATURE  GENERATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Tas(i),T_tilde(i),dT(i),qt,I2,I3,I4]=ComputeAirTemperature(Datam(i,:),DeltaGMT,...
            Lon,Lat,Ns(i),Par.bTemp(j,:),0,Par.dTbar(j,Hor),Par.rhodT(j),Par.sigmadT(j,Hor),...
            0,Par.Ti(j),0,0,0,Par.Ti(j));
         Ti = Par.Ti(j); 
    else
        if (Datam(i,4) == 0) 
            Ti = T_tilde(i-1);
        end

        [Tas(i),T_tilde(i),dT(i),qt,I2,I3,I4]=ComputeAirTemperature(Datam(i,:),DeltaGMT,...
            Lon,Lat,Ns(i),Par.bTemp(j,:),qt,Par.dTbar(j,Hor),Par.rhodT(j),Par.sigmadT(j,Hor),...
            dT(i-1),Ti,I2,I3,I4,Par.Ti(j));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Tas(i)= Tas(i) + DT_fut(j); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%   VAPOR PRESSURE  GENERATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    esat_s(i)=611*exp(17.27*Tas(i)/(237.3+Tas(i))); %%
    if i< 3
        [eas(i),dDe(i)] = ComputeVapPressure(esat_s(i),Tas(i),0,0,0,Par.aVap(j,:),...
            Par.dDem(j),Par.rhodDe(j),Par.sigmadDe(j));
    else
        [eas(i),dDe(i)] = ComputeVapPressure(esat_s(i),Tas(i),Rsws(i-1),Rsws(i-2),dDe(i-1),Par.aVap(j,:),...
            Par.dDem(j),Par.rhodDe(j),Par.sigmadDe(j));
    end
    Us(i)=eas(i)/esat_s(i); 
    G_g= a_a*Tas(i)./(b_b+Tas(i)) + log(Us(i));
    Tdews(i) =b_b*G_g./(a_a-G_g);
    Tdews(isinf(Tdews))=NaN; clear G_g 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%   RADIATION   GENERATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [SB(i),SD(i),SAB1(i),SAB2(i),SAD1(i),SAD2(i),PARB(i),PARD(i)]=ComputeRadiationForcings(Datam(i,:),DeltaGMT,Lon,Lat,...
        Par.LWP_R(j),Ns(i),Zbas,Tdews(i),Par.beta_A(j),Par.alpha_A,Par.omega_A1,Par.omega_A2,Par.uo,Par.un,Par.rho_g);
    Rsws(i)=SB(i)+SD(i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%   WIND   GENERATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if  i < 4 
       [Wss(i),dWs(i)]=ComputeWindSpeed(0,0,0,Rsws(i),Par.cWind,Par.EdWs,0,Par.rhodWs,Par.sigmadWs,Par.skedWs);
    else
       [Wss(i),dWs(i)]=ComputeWindSpeed(Rsws(i-3),Rsws(i-2),Rsws(i-1),Rsws(i),Par.cWind,Par.EdWs,dWs(i-1),Par.rhodWs,Par.sigmadWs,Par.skedWs);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%   PRESSURE  GENERATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if i==1
        [Pres(i)]=ComputeAtmPressure(Par.EPre,Par.rhopre,Par.EPre,Par.sigmapre);
    else
        [Pres(i)]=ComputeAtmPressure(Par.EPre,Par.rhopre,Pres(i-1),Par.sigmapre); 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Tas(i) > 100) || (Tas(i) < -100)
        errordlg('Numerical instability in the air temperature estimation','Error03')
        return
    end
    %%%%%%%%%%%%%%%%%%%%%%5
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
close(bau)
end 
%%%%%%%%%%