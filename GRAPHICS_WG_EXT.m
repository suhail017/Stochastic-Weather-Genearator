%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% GRAPHICS WEATHER GENERATOR  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pr_yr_s,Ta_yr_s]= GRAPHICS_WG_EXT(DS,Par,NYR,Psim,Pyr_AR_orig,Ns,Tas,eas,Us,Tdews,Rsws,...
    Rdirs,Rdifs,Wss,Pres,D,Pr,Pr_yr,Ta,Ta_yr,N,Tdew,U,ea,Rsw,Rdir,Rdif,Ws,Pre,...
    Prm, Nm,Tam, Rdifm, Rdirm, Rswm,Um, Tdewm,eam)
%%%% INPUT
%Prm Prms PH  Psim Pr Pr_yr Pr_yr_s  Pyr_AR_orig NYR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Yro,Moo,Dao,Hro]=datevec(D);
[Yrs,Mo,Das,Hrs]=datevec(DS);
%HI = hour(DS(1)); %%% HI start hour
HI = Hrs(1); %%% HI start hour
%Yrs= year(DS); 
%Mo = month(DS); 
%HIm = hour(D(1)); %%% HI start hour
HIm = Hro(1); %%% HI start hour
while  not(isequal(HIm,HI))
    D=D(2:end);
    Ta=Ta(2:end);
    Pr=Pr(2:end);
    N=N(2:end);
    Tdew=Tdew(2:end);
    U=U(2:end);
    ea=ea(2:end);
    Rsw=Rsw(2:end); 
    Rdir=Rdir(2:end);
    Rdif=Rdif(2:end);
    Ws=Ws(2:end);
    Pre=Pre(2:end);
    %HIm = hour(D(1)); %%% HI start hour
    Hro=Hro(2:end);
    HIm = Hro(1); %%% HI start hour
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ta_n=Ta(not(isnan(Ta)));
Rsw_n=Rsw(not(isnan(Rsw)));
U_n=U(not(isnan(U)));
ea_n=ea(not(isnan(ea)));
Tdew_n=Tdew(not(isnan(Tdew)));
Ws_n=Ws(not(isnan(Ws)));
Pre_n=Pre(not(isnan(Pre)));
N_n=N(not(isnan(N)));
%%%     J  F  M  A  M  J  J  A  S  O  N  D
Days = [31 28 31 30 31 30 31 31 30 31 30 31];
dt= 1; %% % time step 1-Hour
didascalia={'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'};
mm=1:12;
%%%%%%%%%%%%%%%%%
%%% Monthly Subdvision  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=0;
for i=1:12
    Prsm{i}= Psim(Mo==i) ;
    Nsm{i}= Ns(Mo==i) ;
    Rswsm{i}= Rsws(Mo==i);
    Rdirsm{i}=Rdirs(Mo==i);
    Rdifsm{i}=Rdifs(Mo==i);
    Tasm{i}=Tas(Mo==i);
    Usm{i}=Us(Mo==i);
    Tdewsm{i}=Tdews(Mo==i);
    easm{i}=eas(Mo==i);
end
%%%%% Year Computation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_TR = 8520; %% Threshold Year Validity %%% 355 Days
r1=0; r2=0;
for i=min(Yrs):max(Yrs)
    if length(Psim(Yrs==i)) >= H_TR
        r1=r1+1;
        Pr_yr_s(r1)= nansum(Psim(Yrs==i));
    end
    if length(Tas(Yrs==i)) >= H_TR
        r2=r2+1;
        Ta_yr_s(r2)= nanmean(Tas(Yrs==i));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
clear r1 r2 i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    RAINFALL GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:12
    fr=24*Days(i);
    %%%%%%%%%%%%%%%%%%
    [EPM(i),VARPM(i),CVPM(i),RlPM(i),SKPM(i),FFh,Fddh,Fwwh]=sample_properties(Prm{i},dt,fr,1);
    [EPS(i),VARPS(i),CVPS(i),RlPS(i),SKPS(i),FFh,Fddh,Fwwh]=sample_properties(Prsm{i},dt,fr,1);
end
%%%%%%%%%%%%%%%%%%%%%
figure(101)
set(gca,'FontSize',12);
errorbar(mm,EPM,sqrt(VARPM),'--or','LineWidth', 2.0); hold on;
errorbar(mm,EPS,sqrt(VARPS),'og','LineWidth', 2.0);
legend('OBS.','SIM.')
xlabel('Month'); ylabel('[mm]')
title('Monthly Precipitation')
grid on
axis([0.5 12.5 min(EPM-max(sqrt(VARPM)))-10 max(EPM+max(sqrt(VARPM)))+10])
%%%%%%%%%%%%%%%%
%%%%%%%%%%
for ii=1:12
    [EP24(ii),VARP24(ii),CVP24(ii),RlP24(ii),SKP24(ii),FFP24(ii),FddP24(ii),FwwP24(ii)]=sample_properties(Prm{ii},dt,24,1);
    [EP1(ii),VARP1(ii),CVP1(ii),RlP1(ii),SKP1(ii),FFP1(ii),FddP1(ii),FwwP1(ii)]=sample_properties(Prm{ii},dt,1,1);
    [EP48(ii),VARP48(ii),CVP48(ii),RlP48(ii),SKP48(ii),FFP48(ii),FddP48(ii),FwwP48(ii)]=sample_properties(Prm{ii},dt,48,1);
    %%%
    [E24(ii),VAR24(ii),CV24(ii),Rl24(ii),SK24(ii),FF24(ii),Fdd24(ii),Fww24(ii)]=sample_properties(Prsm{ii},1,24,1);
    [E1(ii),VAR1(ii),CV1(ii),Rl1(ii),SK1(ii),FF1(ii),Fdd1(ii),Fww1(ii)]=sample_properties(Prsm{ii},1,1,1);
    [E48(ii),VAR48(ii),CV48(ii),Rl48(ii),SK48(ii),FF48(ii),Fdd48(ii),Fww48(ii)]=sample_properties(Prsm{ii},1,48,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mm=1:12;
figure(102)
subplot(3,2,1); set(gca,'FontSize',8.5);
plot(mm,E1,'g',mm,EP1,'r','LineWidth', 2);
grid on ; title('a) Mean')
ylabel('[mm]')
%text(15,0.45,['1 Hour'],'hor','left','vert','top','FontSize',12);
subplot(3,2,2); set(gca,'FontSize',8.5);
plot(mm,VAR1,'g',mm,VARP1,'r','LineWidth', 2);
grid on ; title('b) Variance')
ylabel('[mm^2]')
subplot(3,2,3); set(gca,'FontSize',8.5);
plot(mm,Rl1,'g',mm,RlP1,'r','LineWidth', 2);
grid on ; title('c) Lag-1 autocorrelation')
ylabel('[-]')
subplot(3,2,4); set(gca,'FontSize',8.5);
plot(mm,SK1,'g',mm,SKP1,'r','LineWidth', 2);
grid on ; title('d) Skewness')
ylabel('[-]')
subplot(3,2,5); set(gca,'FontSize',8.5);
plot(mm,FF1,'g',mm,FFP1,'r','LineWidth', 2);
grid on ; title('e) Frequency of non-precipitation')
xlabel('Month')
ylabel('[-]')
subplot(3,2,6); set(gca,'FontSize',8.5);
plot(mm,Fww1,'g',mm,FwwP1,'r','LineWidth', 2);
grid on ; title('f) Transition probability wet-wet')
legend('SIM.','OBS.')
xlabel('Month')
ylabel('[-]')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(103)
subplot(3,2,1); set(gca,'FontSize',8.5);
plot(mm,E24,'g',mm,EP24,'r','LineWidth', 2);
grid on ;  title('a) Mean')
ylabel('[mm]')
subplot(3,2,2); set(gca,'FontSize',8.5);
plot(mm,VAR24,'g',mm,VARP24,'r','LineWidth', 2);
grid on ; title('b) Variance')
ylabel('[mm^2]')
subplot(3,2,3); set(gca,'FontSize',8.5);
plot(mm,Rl24,'g',mm,RlP24,'r','LineWidth', 2);
grid on ; title('c) Lag-1 autocorrelation')
ylabel('[-]')
subplot(3,2,4); set(gca,'FontSize',8.5);
plot(mm,SK24,'g',mm,SKP24,'r','LineWidth', 2);
grid on ;title('d) Skewness')
ylabel('[-]')
subplot(3,2,5); set(gca,'FontSize',8.5);
plot(mm,FF24,'g',mm,FFP24,'r','LineWidth', 2);
grid on ;title('e) Frequency of non-precipitation')
xlabel('Month')
ylabel('[-]')
subplot(3,2,6); set(gca,'FontSize',8.5);
plot(mm,Fww24,'g',mm,FwwP24,'r','LineWidth', 2);
grid on ; title('f) Transition probability wet-wet')
xlabel('Month')
legend('SIM.','OBS.')
ylabel('[-]')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(104)
subplot(3,2,1); set(gca,'FontSize',8.5);
plot(mm,E48,'g',mm,EP48,'r','LineWidth', 2);
grid on ; title('a) Mean')
ylabel('[mm]')
subplot(3,2,2); set(gca,'FontSize',8.5);
plot(mm,VAR48,'g',mm,VARP48,'r','LineWidth', 2);
grid on ; title('b) Variance')
ylabel('[mm^2]')
subplot(3,2,3); set(gca,'FontSize',8.5);
plot(mm,Rl48,'g',mm,RlP48,'r','LineWidth', 2);
grid on ; title('c) Lag-1 autocorrelation')
ylabel('[-]')
subplot(3,2,4); set(gca,'FontSize',8.5);
plot(mm,SK48,'g',mm,SKP48,'r','LineWidth', 2);
grid on ; title('d) Skewness')
ylabel('[-]')
subplot(3,2,5); set(gca,'FontSize',8.5);
plot(mm,FF48,'g',mm,FFP48,'r','LineWidth', 2);
grid on ; title('e) Frequency of non-precipitation')
xlabel('Month')
ylabel('[-]')
subplot(3,2,6); set(gca,'FontSize',8.5);
plot(mm,Fww48,'g',mm,FwwP48,'r','LineWidth', 2);
grid on ; title('f) Transition probability wet-wet')
xlabel('Month')
ylabel('[-]')
legend('SIM.','OBS.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%
dur=[ 1 3 6 9 12 18 24 48 72 96];
i=0; jk=1; k=7;
for i=1:length(dur)
    [N0(i),N1(i),N10(i),N20(i),Xse{i},ttre{i},pvre{i}]=sample_extreme_prop(Psim,1,dur(i));
    [N0P(i),N1P(i),N10P(i),N20P(i),XseP{i},ttreP{i},pvreP{i}]=sample_extreme_prop(Pr,dt,dur(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(105)
subplot(1,3,1); set(gca,'FontSize',10);
plot(dur,N1,'g',dur,N1P,'r','LineWidth', 2);
text(dur(8),N1(8)+0.1,['> 1 [mm]'],'hor','left','vert','top','FontSize',11);
hold on ;
plot(dur,N10,'g',dur,N10P,'r','LineWidth', 2);
text(dur(8),N10(8)+0.1,['> 10 [mm]'],'hor','left','vert','top','FontSize',11);
plot(dur,N20,'g',dur,N20P,'r','LineWidth', 2);
text(dur(8),N20(8)+0.1,['> 20 [mm]'],'hor','left','vert','top','FontSize',11);
%plot(dur,N0,'r',dur,N0P,'g','LineWidth', 2);
%text(dur(8),N0(8)+0.1,['> 0 [mm]'],'hor','left','vert','top','FontSize',11);
xlabel('Agg. Period [h]'); ylabel('Fraction')
grid on ; title('a) Fraction of time precipitation larger than # [mm]')
%%%%%%%%%%%%%%%%
figure(106)
subplot(4,1,1); set(gca,'FontSize',9);
semilogx(ttre{jk},Xse{jk},'xg',ttreP{jk},XseP{jk},'xr','LineWidth', 1.5);
grid on ; title('a) Extremes of precipitation 1 hour')
ylabel('[mm]')
figure(106)
subplot(4,1,2);set(gca,'FontSize',9);
semilogx(ttre{k},Xse{k},'xg',ttreP{k},XseP{k},'xr','LineWidth', 1.5);
grid on ; title('b) Extremes of precipitation 24 hours')
ylabel('[mm]')
legend('SIM.','OBS.')
%%%Prsm 1h  %% PPP dt
%[DSs,Drys,ttres,pvres]=sample_dry_extreme(Psim,1);
%[DSm,Drym,ttrem,pvrem]=sample_dry_extreme(PPP,dt);
TH=1; %%[mm] Threshold
[DSs,WSs,Drys,Wets,ttres,pvres]=sample_dry_and_wet_extreme(Psim,1,TH);
[DSm,WSm,Drym,Wetm,ttrem,pvrem]=sample_dry_and_wet_extreme(Pr,dt,TH);
%%%%%%%%%%%%%%%%%%
figure(105)
edges=[1:2:15]; 
subplot(1,3,2); set(gca,'FontSize',10);
[hx]=histc(DSm,edges);
[hs]=histc(DSs,edges);
bar(edges,hx/length(DSm),'r','LineWidth', 1);
hold on ; grid on;
plot(edges,hs/length(DSs),'+-g','LineWidth', 2);
title('b) Dry spell length distribution ')
xlabel('Consecutive days'); ylabel('Frequency')
legend('OBS.','SIM')
text(7,0.31,['E_{obs}=',num2str(mean(DSm))],'hor','left','vert','top','FontSize',8.5);
text(7,0.26,['E_{sim}=',num2str(mean(DSs))],'hor','left','vert','top','FontSize',8.5);
text(7,0.19,['\sigma_{obs}=',num2str(std(DSm))],'hor','left','vert','top','FontSize',8.5);
text(7,0.16,['\sigma_{sim}=',num2str(std(DSs))],'hor','left','vert','top','FontSize',8.5);
figure(106)
subplot(4,1,3); set(gca,'FontSize',9);
semilogx(ttres,Drys,'xg',ttrem,Drym,'xr','LineWidth', 1.5);
grid on ; title('c) Extreme dry spell, consecutive days P_r < 1 [mm]')
xlabel('Return period'); ylabel('Days')
%%%%%%%%%%%%%%%%%%
figure(105)
edges=[1:2:10]; 
subplot(1,3,3); set(gca,'FontSize',10);
[hx]=histc(WSm,edges);
[hs]=histc(WSs,edges);
bar(edges,hx/length(WSm),'r','LineWidth', 1);
hold on ; grid on;
plot(edges,hs/length(WSs),'+-g','LineWidth', 2);
title('c) Wet spell length distribution')
xlabel('Consecutive days'); ylabel('Frequency')
text(4,0.5,['E_{obs}=',num2str(mean(WSm))],'hor','left','vert','top','FontSize',8.5);
text(4,0.45,['E_{sim}=',num2str(mean(WSs))],'hor','left','vert','top','FontSize',8.5);
text(4,0.35,['\sigma_{obs}=',num2str(std(WSm))],'hor','left','vert','top','FontSize',8.5);
text(4,0.3,['\sigma_{sim}=',num2str(std(WSs))],'hor','left','vert','top','FontSize',8.5);
figure(106)
subplot(4,1,4); set(gca,'FontSize',9);
semilogx(ttres,Wets,'xg',ttrem,Wetm,'xr','LineWidth', 1.5);
grid on ; title('d) Extreme wet spell, consecutive days P_r >= 1 [mm]')
xlabel('Return period'); ylabel('Days')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MEr=15;
figure(108)
set(gca,'FontSize',12);
%plot(Pr_yr,'g','LineWidth', 2.5);
errorbar(1:NYR,Pyr_AR_orig,ones(1,NYR)*MEr,'--om','LineWidth', 2.5);
hold on; grid on;
plot(Pr_yr_s,'r','LineWidth', 2.5);
xlabel('Year'); title('Total annual precipitation'); ylabel('[mm]') 
legend('Simulated AR(1)','Simulated NSRP')
axis([0 NYR+1 min(Pyr_AR_orig)-MEr-25 max(Pyr_AR_orig)+MEr+25 ])
%text(NYR/2,max(Pyr_AR_orig)+MErorig+5,['Ss  ',num2str(mean(Pyr)],'hor','left','vert','top','FontSize',8.5);
%text(0.2,0.6,['Em
%',num2str(nanmean(Pr_yr_s))],'hor','left','vert','top','FontSize',8.5);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   CLOUD COVER GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:12
    SS=Prm{j}>0;
    SSs = Prsm{j}>0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Nfw_s] = Fair_weather_cloudiness(Nsm{j},SSs,Par.Tr(j));
    [Nfw] = Fair_weather_cloudiness(Nm{j},SS,Par.Tr(j));
    figure(201)
    %%%%%title('Total Cloudiness')
    subplot(6,2,j); set(gca,'FontSize',9);
    [nn, hh] = hist(Nm{j},10);  nn = nn/length(Nm{j});
    [nns, hhs] = hist(Nsm{j},10);  nns = nns/length(Nsm{j});
    bar(hh, nn, 0.75); grid on; hold on; colormap cool
    plot(hhs,nns,'-+m','LineWidth', 2)
    text(0.2,0.6,['E_{obs}=',num2str(nanmean(Nm{j}))],'hor','left','vert','top','FontSize',6.5);
    text(0.6,0.6,['E_{sim}=',num2str(nanmean(Nsm{j}))],'hor','left','vert','top','FontSize',6.5);
    text(0.2,0.35,['\sigma_{obs}=',num2str(nanstd(Nm{j}))],'hor','left','vert','top','FontSize',6.5);
    text(0.6,0.35,['\sigma_{sim}=',num2str(nanstd(Nsm{j}))],'hor','left','vert','top','FontSize',6.5);
    ey=max(0.65,max(max(nn),max(nns)))+0.1;
    axis([0 1 0 ey])
    title(didascalia{j});
    if j>=11
        xlabel('Cloudiness')
        if j>=12 
        legend('OBS.','SIM.')
        end 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%
    figure(202)
    %%%%%title('Fair Weather Cloudiness')
    subplot(6,2,j); set(gca,'FontSize',9);
    [nn, hh] = hist(Nfw,10);  nn = nn/length(Nfw);
    [nns, hhs] = hist(Nfw_s,10);  nns = nns/length(Nfw_s);
    bar(hh, nn, 0.75); grid on; hold on; colormap cool
    plot(hhs,nns,'-+m','LineWidth', 2)
    text(0.2,0.6,['E_{obs}=',num2str(nanmean(Nfw))],'hor','left','vert','top','FontSize',6.5);
    text(0.6,0.6,['E_{sim}=',num2str(nanmean(Nfw_s))],'hor','left','vert','top','FontSize',6.5);
    text(0.2,0.35,['\sigma_{obs}=',num2str(nanstd(Nfw))],'hor','left','vert','top','FontSize',6.5);
    text(0.6,0.35,['\sigma_{sim}=',num2str(nanstd(Nfw_s))],'hor','left','vert','top','FontSize',6.5);
    ey=max(0.65,max(max(nn),max(nns)))+0.1;
    axis([0 1 0 ey])
    title(didascalia{j});
    if j>=11
        xlabel('Cloudiness')
        if j>=12 
        legend('OBS.','SIM.')
        end 
    end
end
[hx,nn]=hist(N_n,10); 
[hs,nns]=hist(Ns,10);
figure(203)
set(gca,'FontSize',11.5);
bar(nn,hx/length(N_n),0.75,'c','LineWidth',2)
hold on ; grid on;
plot(nns,hs/length(Ns),'-+k','LineWidth',2)
legend('OBS.','SIM.')
xlabel('N [-]'); ylabel('Frequency')
title('PDF Cloudiness')
text(0.2,0.4,['E_{obs}=',num2str(nanmean(N))],'hor','left','vert','top','FontSize',10);
text(0.6,0.4,['E_{sim}=',num2str(nanmean(Ns))],'hor','left','vert','top','FontSize',10);
text(0.2,0.32,['\sigma_{obs}=',num2str(nanstd(N))],'hor','left','vert','top','FontSize',10);
text(0.6,0.32,['\sigma_{sim}=',num2str(nanstd(Ns))],'hor','left','vert','top','FontSize',10);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  AIR TEMPERATURE  GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:12
    %%%%%%%%%%%%%%%%%%%%
    [Eh_1s(j),VARh_1s(j),CVh_1s(j),Rlh_1s(j),SKh_1s(j)]=sample_properties2( Tasm{j},1,1,1);
    [Eh_24s(j),VARh_24s(j),CVh_24s(j),Rlh_24s(j),SKh_24s(j)]=sample_properties2( Tasm{j},1,24,1);
    %[E_max_1s(j),E_min_1s(j),Std_max_1s(j),Std_min_1s(j)]=sample_properties_mm( Tasm{j},1,1);
    [E_max_24s(j),E_min_24s(j),Std_max_24s(j),Std_min_24s(j)]=sample_properties_mm( Tasm{j},1,24);
    %%%%%%%%%%%
    [Eh_1(j),VARh_1(j),CVh_1(j),Rlh_1(j),SKh_1(j)]=sample_properties2(Tam{j},dt,1,1);
    [Eh_24(j),VARh_24(j),CVh_24(j),Rlh_24(j),SKh_24(j)]=sample_properties2(Tam{j},dt,24,1);
    %[E_max_1(j),E_min_1(j),Std_max_1(j),Std_min_1(j)]=sample_properties_mm(Tam{j},dt,1);
    [E_max_24(j),E_min_24(j),Std_max_24(j),Std_min_24(j)]=sample_properties_mm(Tam{j},dt,24);
    %%%%%%%%%%%%%%%%5
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(301)
subplot(2,1,1); set(gca,'FontSize',9.5);
errorbar(mm,Eh_1,sqrt(VARh_1),'or','LineWidth', 2);
hold on;
errorbar(mm,Eh_1s,sqrt(VARh_1s),'--og','LineWidth', 2);
xlabel('Month'); ylabel('[°C]')
title('a) Monthly average air temperature, agg. period 1 hour')
legend('OBS.','SIM.')
grid on ;
figure(301)
subplot(2,1,2); set(gca,'FontSize',9.5);
errorbar(mm,Eh_24,sqrt(VARh_24),'or','LineWidth', 2);
hold on;
errorbar(mm,Eh_24s,sqrt(VARh_24s),'--og','LineWidth', 2);
xlabel('Month'); ylabel('[°C]')
grid on
title('b) Monthly average air temperature, agg. period 24 hour')
legend('OBS.','SIM.')
figure(302)
subplot(2,1,1); set(gca,'FontSize',9.5);
errorbar(mm,E_max_24,Std_max_24,'or','LineWidth', 2);
hold on;
errorbar(mm,E_max_24s,Std_max_24s,'--og','LineWidth', 2);
xlabel('Month'); ylabel('[°C]')
title('a) Monthly average daily maximum air temperature')
grid on
legend('OBS.','SIM.')
figure(302)
subplot(2,1,2); set(gca,'FontSize',9.5);
errorbar(mm,E_min_24,Std_min_24,'or','LineWidth', 2);
hold on;
errorbar(mm,E_min_24s,Std_min_24s,'--og','LineWidth', 2);
xlabel('Month'); ylabel('[°C]')
title('b) Monthly average daily minimum air temperature')
grid on
legend('OBS.','SIM.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edges=[min(min(Tas),min(Ta_n)):3:max(max(Tas),max(Ta_n))];
[hx]=histc(Ta_n,edges);
[hs]=histc(Tas,edges);
%[Acx,lags,Bounds]=autocorr(Ta_n,150);
%[Acss,lags,Bounds]=autocorr(Tas,150);
[Acx,lags]=xcov(Ta_n,Ta_n,150,'coeff'); %% 
[Acss,lags]=xcov(Tas,Tas,150,'coeff'); %% 
figure(304)
subplot(1,2,1); set(gca,'FontSize',10.5);
bar(edges,hx/length(Ta_n),0.75,'r','LineWidth',1)
hold on ; grid on;
plot(edges,hs/length(Tas),'+-g','LineWidth',2)
legend('OBS.','SIM.')
xlabel('[°C]'); ylabel('Frequency')
title('a) PDF Air temperature')
axis([min(min(Tas),min(Ta_n)) max(max(Tas),max(Ta_n)) 0  max(hx/length(Ta_n))+0.01]);
text(5,0.02,['E_{obs}=',num2str(nanmean(Ta))],'hor','left','vert','top','FontSize',8.5);
text(5,0.04,['E_{sim}=',num2str(nanmean(Tas))],'hor','left','vert','top','FontSize',8.5);
text(5,0.06,['\sigma_{obs}=',num2str(nanstd(Ta))],'hor','left','vert','top','FontSize',8.5);
text(5,0.08,['\sigma_{sim}=',num2str(nanstd(Tas))],'hor','left','vert','top','FontSize',8.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=length(Ta); fr=24;
n1=length(Tas);
%%%%%%%%%%%%%%%%%%%%
m=floor(n/fr);m1=floor(n1/fr);
Tap=reshape(Ta(1:m*fr),fr,m);
Tasp=reshape(Tas(1:m1*fr),fr,m1);
Xm=nanmean(Tap');
Xs=mean(Tasp');
stdXm=nanstd(Tap');
stdXs=std(Tasp');
t=1:fr;
figure(304)
subplot(1,2,2); set(gca,'FontSize',10.5);
plot(t,Xm,'r','LineWidth',2)
hold on ; grid on;
plot(t,Xs,'g','LineWidth',2)
plot(t,stdXm,'^r')
plot(t,stdXs,'^g')
legend('OBS.','SIM.')
ylabel('[°C]'); xlabel('Hour')
title('b) Temperature daily cycle')
clear Tap Tasp Xm Xs stdXm stdXs
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Analysis on the year time series not monthly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Tmax_m,tr_m,pvr_m]=sample_extreme_prop2(Ta_n,dt,1);
[Tmax_s,tr_s,pvr_s]=sample_extreme_prop2(Tas,1,1);
[Tmin_m,tr_m,pvr_m]=sample_extreme_prop2(-Ta_n,dt,1);
[Tmin_s,tr_s,pvr_s]=sample_extreme_prop2(-Tas,1,1);
figure(307)
subplot(4,1,1); set(gca,'FontSize',9.1);
plot(tr_m,Tmax_m,'xr',tr_s,Tmax_s,'xg','LineWidth',2)
grid on ; title('a) Extremes of maximum air temperature, 1 hour ')
ylabel('[°C]'); 
%legend('OBS.','SIM.')
figure(307)
subplot(4,1,2); set(gca,'FontSize',9.1);
plot(tr_m,-Tmin_m,'xr',tr_s,-Tmin_s,'xg','LineWidth',2)
grid on ; title('b) Extremes of minimum air temperature,  1 hour ')
ylabel('[°C]'); 
%legend('OBS.','SIM.')
[Tmax_m,tr_m,pvr_m]=sample_extreme_prop2(Ta_n,dt,24);
[Tmax_s,tr_s,pvr_s]=sample_extreme_prop2(Tas,1,24);
[Tmin_m,tr_m,pvr_m]=sample_extreme_prop2(-Ta_n,dt,24);
[Tmin_s,tr_s,pvr_s]=sample_extreme_prop2(-Tas,1,24);
figure(307)
subplot(4,1,3); set(gca,'FontSize',9.1);
plot(tr_m,Tmax_m,'xr',tr_s,Tmax_s,'xg','LineWidth',2)
grid on ; title('c) Extremes of maximum air temperature, 24 hours')
ylabel('[°C]'); 
%legend('OBS.','SIM.')
figure(307)
subplot(4,1,4); set(gca,'FontSize',9.1);
plot(tr_m,-Tmin_m,'xr',tr_s,-Tmin_s,'xg','LineWidth',2)
grid on ; title('d) Extremes of minimum air temperature, 24 hours')
xlabel('Return period'); ylabel('[°C]'); 
legend('OBS.','SIM.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Dfr_m,Dwa_m,Nfr_m,Nwa_m,HWe_m,CWe_m,tre_m,pvre_m]=sample_temp_extreme(Ta_n,dt);
[Dfr_s,Dwa_s,Nfr_s,Nwa_s,HWe_s,CWe_s,tre_s,pvre_s]=sample_temp_extreme(Tas,1);
figure(308)
subplot(2,1,1); set(gca,'FontSize',9.5);
plot(tre_m,HWe_m,'xr',tre_s,HWe_s,'xg','LineWidth',2)
grid on ; title('a) Extreme heat waves. Consecutive days with air temp. larger than 90 percentile')
xlabel('Return period'); ylabel('Days')
legend('OBS.','SIM.')
figure(308)
subplot(2,1,2); set(gca,'FontSize',9.5);
plot(tre_m,CWe_m,'xr',tre_s,CWe_s,'xg','LineWidth',2)
grid on ; title('b) Extreme cold waves. Consecutive days with air temp. lower than 10 percentile')
%xlabel('Return period'); 
ylabel('Days')
legend('OBS.','SIM.')
%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  RADIATION   GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:12
    [Edir(j),VARdir(j),CV(j),Rl(j),SK(j)]=sample_properties2(Rdirm{j},1,24,1);
    [Edif(j),VARdif(j),CV(j),Rl(j),SK(j)]=sample_properties2(Rdifm{j},1,24,1);
    [Etot(j),VARtot(j),CV(j),Rl(j),SK(j)]=sample_properties2(Rswm{j},1,24,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Edir_s(j),VARdir_s(j),CV(j),Rl(j),SK(j)]=sample_properties2(Rdirsm{j},1,24,1);
    [Edif_s(j),VARdif_s(j),CV(j),Rl(j),SK(j)]=sample_properties2(Rdifsm{j},1,24,1);
    [Etot_s(j),VARtot_s(j),CV(j),Rl(j),SK(j)]=sample_properties2(Rswsm{j},1,24,1);%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(401)
set(gca,'FontSize',10);
subplot(3,1,2);
errorbar(mm,Edir,sqrt(VARdir),'--or','LineWidth', 2);
hold on;
errorbar(mm,Edir_s,sqrt(VARdir_s),'--og','LineWidth', 2);
ylabel('[W/m^2]')
grid on
title('b) Direct shortwave radiation')
figure(401)
subplot(3,1,3);
set(gca,'FontSize',10);
errorbar(mm,Edif,sqrt(VARdif),'--or','LineWidth', 2);
hold on;
errorbar(mm,Edif_s,sqrt(VARdif_s),'--og','LineWidth', 2);
xlabel('Month'); ylabel('[W/m^2]')
grid on
title('c) Diffuse shortwave radiation')
figure(401)
subplot(3,1,1);
set(gca,'FontSize',10);
errorbar(mm,Etot,sqrt(VARtot),'--or','LineWidth', 2);
hold on;
errorbar(mm,Etot_s,sqrt(VARtot_s),'--og','LineWidth', 2);
 ylabel('[W/m^2]')
title('a) Global shortwave radiation')
grid on
legend('OBS.','SIM.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
n=length(Rsw); fr=24/dt;
n1=length(Rsws);
%%%%%%%%%%%%%%%%%%%%
m=floor(n/fr); m1=floor(n1/fr); 
Rswp=reshape(Rsw(1:m*fr),fr,m);
Rswsp=reshape(Rsws(1:m1*fr),fr,m1);
Xm=nanmean(Rswp');
Xs=mean(Rswsp');
stdXm=nanstd(Rswp');
stdXs=std(Rswsp');
t=[HI:23, 0:HI-1];
figure(404)
subplot(1,3,1);
set(gca,'FontSize',10);
plot(t,Xm,'r','LineWidth', 2);
hold on ; grid on;
plot(t,Xs,'b','LineWidth', 2);
title('a) Global shortwave radiation daily cycle')
%legend('OBS.','SIM.')
xlabel('Hour'); ylabel('[W/m^2]')
mR=max(max(Xm,Xs))+50;
axis([0 24 0 mR])
clear Tap Tasp Xm Xs stdXm stdXs
%%%%%%%%%%%%%%%%%%%%%
n=length(Rdif); fr=24/dt;
n1=length(Rdifs);
%%%%%%%%%%%%%%%%%%%%
m=floor(n/fr);m1=floor(n1/fr); 
Rdifp=reshape(Rdif(1:m*fr),fr,m);
SDp=reshape(Rdifs(1:m1*fr),fr,m1);
Xm=nanmean(Rdifp');
Xs=mean(SDp');
stdXm=nanstd(Rdifp');
stdXs=std(SDp');
figure(404)
subplot(1,3,3);
set(gca,'FontSize',10);
plot(t,Xm,'r','LineWidth', 2);
hold on ; grid on;
plot(t,Xs,'b','LineWidth', 2);
title('c) Diffuse shortwave radiation daily cycle')
legend('OBS.','SIM.')
xlabel('Hour'); ylabel('[W/m^2]')
axis([0 24 0 mR])
clear Tap Tasp Xm Xs stdXm stdXs
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
n=length(Rdif); fr=24/dt;
n1=length(Rdifs);
%%%%%%%%%%%%%%%%%%%%
m=floor(n/fr); m1=floor(n1/fr);
Rdirp=reshape(Rdir(1:m*fr),fr,m);
SBp=reshape(Rdirs(1:m1*fr),fr,m1);
Xm=nanmean(Rdirp');
Xs=mean(SBp');
stdXm=nanstd(Rdirp');
stdXs=std(SBp');
figure(404)
subplot(1,3,2);
set(gca,'FontSize',10);
plot(t,Xm,'r','LineWidth', 2);
hold on ; grid on;
plot(t,Xs,'b','LineWidth', 2);
title('b) Direct shortwave radiation daily cycle')
%legend('OBS.','SIM.')
xlabel('Hour'); ylabel('[W/m^2]')
axis([0 24 0 mR])
clear Tap Tasp Xm Xs stdXm stdXs
%%%%%%%%%%%%%%%%%%%%%
[Yro,Moo,Dao,Hro]=datevec(D);
Mom=Moo; Dam=Dao; RswPlot=Rsw;
while  not(isequal(Dam(1),Das(1))) ||  not(isequal(Mom(1),Mo(1)))
    RswPlot=RswPlot(2:end);
    Mom = Mom(2:end); %%% HI start hour
    Dam= Dam(2:end);
end
%%%%%%%%%%%%%%%%%%
n=length(RswPlot); fr=24*365/dt;
n1=length(Rsws);
%%%%%%%%%%%%%%%%%%%%
m=floor(n/fr);m1=floor(n1/fr);
Rswp=reshape(RswPlot(1:m*fr),fr,m);
Rswsp=reshape(Rsws(1:m1*fr),fr,m1);
Xm=nanmean(Rswp');
Xs=mean(Rswsp');
stdXm=nanstd(Rswp');
stdXs=std(Rswsp');
t=1:fr;
%%%%%%%%%%%%%%%%%%%%%
k=0;Xhm=zeros(365,24);
Xhs=zeros(365,24);
for k=1:24
    i=0;
    while (k+24*(i)) <= length(Xm);
        i=i+1;
        Xhm(i,k)= Xm(k+ 24*(i-1));
        Xhs(i,k)= Xs(k+ 24*(i-1));
    end
end
%%% Choice Hour
GH=mod((HI-1),24);%%  Hour Gap
k=1:24;
for i=1:length(k)
    figure(406)
    subplot(6,4,i);
    set(gca,'FontSize',7.5);
    %if SH > 12
    r= (k(i)+GH)- 24*(k(i)+GH >24);
    plot(1:365,Xhm(:,r),'r','LineWidth', 1);
    hold on ; grid on;
    plot(1:365,Xhs(:,r),'g','LineWidth', 1);
    % title(strcat('Radiation Total Daily cycle fixed Hour =',num2str(k(i))))
    title(strcat('Hour =',num2str(r)))
    axis([0 366 0 max(max(Xhm(:,r)),max(Xhs(:,r)))+100] )
    if i==length(k)
        legend('OBS.','SIM.')
    end
    if i > 20
        xlabel('Julian Day');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% VAPOR PRESSURE GRAPHICS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:12
    %%%%%%%%%%%%%5
    [Eh_1s(j),VARh_1s(j),CVh_1s(j),Rlh_1s(j),SKh_1s(j)]=sample_properties2(Tdewsm{j},1,1,1);
    [Eh_24s(j),VARh_24s(j),CVh_24s(j),Rlh_24s(j),SKh_24s(j)]=sample_properties2(Tdewsm{j},1,24,1);
    [E_max_24s(j),E_min_24s(j),Std_max_24s(j),Std_min_24s(j)]=sample_properties_mm(Tdewsm{j},1,24);
    %%%%
    [Eh_1(j),VARh_1(j),CVh_1(j),Rlh_1(j),SKh_1(j)]=sample_properties2(Tdewm{j},dt,1,1);
    [Eh_24(j),VARh_24(j),CVh_24(j),Rlh_24(j),SKh_24(j)]=sample_properties2(Tdewm{j},dt,24,1);
    [E_max_24(j),E_min_24(j),Std_max_24(j),Std_min_24(j)]=sample_properties_mm(Tdewm{j},dt,24);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [UEh_1s(j),UVARh_1s(j),CVh_1s(j),Rlh_1s(j),SKh_1s(j)]=sample_properties2(Usm{j},1,1,1);
    [UEh_24s(j),UVARh_24s(j),CVh_24s(j),Rlh_24s(j),SKh_24s(j)]=sample_properties2(Usm{j},1,24,1);
    [UEh_1(j),UVARh_1(j),CVh_1(j),Rlh_1(j),SKh_1(j)]=sample_properties2(Um{j},dt,1,1);
    [UEh_24(j),UVARh_24(j),CVh_24(j),Rlh_24(j),SKh_24(j)]=sample_properties2(Um{j},dt,24,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [eaEh_1s(j),eaVARh_1s(j),CVh_1s(j),Rlh_1s(j),SKh_1s(j)]=sample_properties2(easm{j},1,1,1);
    [eaEh_24s(j),eaVARh_24s(j),CVh_24s(j),Rlh_24s(j),SKh_24s(j)]=sample_properties2(easm{j},1,24,1);
    [eaEh_1(j),eaVARh_1(j),CVh_1(j),Rlh_1(j),SKh_1(j)]=sample_properties2(eam{j},dt,1,1);
    [eaEh_24(j),eaVARh_24(j),CVh_24(j),Rlh_24(j),SKh_24(j)]=sample_properties2(eam{j},dt,24,1);
end
figure(503)
subplot(2,1,1); set(gca,'FontSize',8.5);
errorbar(mm,UEh_1,sqrt(UVARh_1),'or','LineWidth', 2);
hold on;
errorbar(mm,UEh_1s,sqrt(UVARh_1s),'--og','LineWidth', 2);
xlabel('Month'); ylabel('[-]')
title('a) Monthly average relative humidity, agg. period 1 hour')
%legend('OBS.','SIM.')
grid on
figure(503)
subplot(2,1,2); set(gca,'FontSize',8.5);
errorbar(mm,UEh_24,sqrt(UVARh_24),'or','LineWidth', 2);
hold on;
errorbar(mm,UEh_24s,sqrt(UVARh_24s),'--og','LineWidth', 2);
xlabel('Month'); ylabel('[-]')
grid on
title('b) Monthly average relative humidity, agg. period 24 hours')
legend('OBS.','SIM.')
figure(504)
subplot(2,1,1); set(gca,'FontSize',8.5);
errorbar(mm,eaEh_1,sqrt(eaVARh_1),'or','LineWidth', 2);
hold on;
errorbar(mm,eaEh_1s,sqrt(eaVARh_1s),'--og','LineWidth', 2);
xlabel('Month'); ylabel('[Pa]')
title('a) Monthly average vapor pressure, agg. period 1 hour')
%legend('OBS.','SIM.')
grid on
figure(504)
subplot(2,1,2); set(gca,'FontSize',8.5);
errorbar(mm,eaEh_24,sqrt(eaVARh_24),'or','LineWidth', 2);
hold on;
errorbar(mm,eaEh_24s,sqrt(eaVARh_24s),'--og','LineWidth', 2);
xlabel('Month'); ylabel('[Pa]')
grid on
title('b) Monthly average vapor pressure, agg. period 24 hours')
legend('OBS.','SIM.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Global Properites
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edges=[min(Tdew_n):3:max(Tdew_n)];
[hx]=histc(Tdew_n,edges);
[hs]=histc(Tdews,edges);
%[Acx,lags,Bounds]=autocorr(Tdew_n,150);
%[Acss,lags,Bounds]=autocorr(Tdews,150);
[Acx,lags]=xcov(Tdew_n,Tdew_n,150,'coeff'); %% 
[Acss,lags]=xcov(Tdews,Tdews,150,'coeff'); %% 
figure(505)
subplot(1,2,1); set(gca,'FontSize',8.5);
bar(edges,hx/length(Tdew_n),0.75,'r','LineWidth', 1)
hold on ; grid on;
plot(edges,hs/length(Tdews),'+-g','LineWidth', 2)
legend('OBS.','SIM.')
xlabel('[°C]'); ylabel('Frequency')
title('a) PDF Dew point temperature')
text(2,0.02,['E_{obs}=',num2str(nanmean(Tdew))],'hor','left','vert','top','FontSize',8.5);
text(2,0.04,['E_{sim}=',num2str(nanmean(Tdews))],'hor','left','vert','top','FontSize',8.5);
text(2,0.06,['\sigma_{obs}=',num2str(nanstd(Tdew))],'hor','left','vert','top','FontSize',8.5);
text(2,0.08,['\sigma_{sim}=',num2str(nanstd(Tdews))],'hor','left','vert','top','FontSize',8.5);
axis([min(Tdew_n) max(Tdew_n) 0  max(hx/length(Tdew_n))+0.01]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hx,nx]=hist(U_n,10);
[hs,ns]=hist(Us,10);
%[Acx,lags,Bounds]=autocorr(U_n,150);
%[Acss,lags,Bounds]=autocorr(Us,150);
[Acx,lags]=xcov(U_n,U_n,150,'coeff'); %% 
[Acss,lags]=xcov(Us,Us,150,'coeff'); %% 
figure(505)
subplot(1,2,2); set(gca,'FontSize',8.5);
bar(nx,hx/length(U_n),0.75,'r','LineWidth', 1)
hold on ; grid on;
plot(ns,hs/length(Us),'+-g','LineWidth', 2)
legend('OBS.','SIM.')
xlabel('[-]'); ylabel('Frequency')
title('b) PDF Relative humidity')
text(0.25,0.02,['E_{obs}=',num2str(nanmean(U))],'hor','left','vert','top','FontSize',8.5);
text(0.25,0.04,['E_{sim}=',num2str(nanmean(Us))],'hor','left','vert','top','FontSize',8.5);
text(0.25,0.06,['\sigma_{obs}=',num2str(nanstd(U))],'hor','left','vert','top','FontSize',8.5);
text(0.25,0.08,['\sigma_{sim}=',num2str(nanstd(Us))],'hor','left','vert','top','FontSize',8.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
n=length(U);
n1=length(Us); fr=24;
%%%%%%%%%%%%%%%%%%%%
m=floor(n/fr); m1=floor(n1/fr);
Tap=reshape(U(1:m*fr),fr,m);
Tasp=reshape(Us(1:m1*fr),fr,m1);
Xm=nanmean(Tap');
Xs=mean(Tasp');
stdXm=nanstd(Tap');
stdXs=std(Tasp');
t=1:fr;
figure(510)
subplot(1,2,1);
set(gca,'FontSize',10.5);
plot(t,Xm,'r','LineWidth', 2.5)
hold on ; grid on;
plot(t,Xs,'g','LineWidth', 2.5)
plot(t,stdXm,'^r')
plot(t,stdXs,'^g')
clear Xm Xs stdXm stdXs Tap Tasp
%Ta=Taa';
%Tap=reshape(Ta(1:m*fr),fr,m);
%Xm=nanmean(Tap');
%plot(t,Xm/30,'y')
%clear Xm Tap
legend('OBS.','SIM.')%,'Temp/30')
xlabel('Hour'); ylabel('[-]')
title('a) Relative humidity daily Cycle ')
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edges=[0:250:max(max(eas),max(ea_n))];
[hx]=histc(ea_n,edges);
[hs]=histc(eas,edges);
%[Acx,lags,Bounds]=autocorr(ea_n,150);
%[Acss,lags,Bounds]=autocorr(eas,150);
[Acx,lags]=xcov(ea_n,ea_n,150,'coeff'); %% 
[Acss,lags]=xcov(eas,eas,150,'coeff'); %% 
figure(510)
subplot(1,2,2);
%subplot(1,2,2); set(gca,'FontSize',8.5);
bar(edges,hx/length(ea_n),0.75,'r','LineWidth', 1)
hold on ; grid on;
plot(edges,hs/length(eas),'+-g','LineWidth', 2)
legend('OBS.','SIM.')
xlabel('[Pa]'); ylabel('Frequency')
title('b) PDF Vapor pressure')
text(1800,0.02,['E_{obs}=',num2str(nanmean(ea))],'hor','left','vert','top','FontSize',8.5);
text(1800,0.04,['E_{sim}=',num2str(nanmean(eas))],'hor','left','vert','top','FontSize',8.5);
text(1800,0.06,['\sigma_{obs}=',num2str(nanstd(ea))],'hor','left','vert','top','FontSize',8.5);
text(1800,0.08,['\sigma_{sim}=',num2str(nanstd(eas))],'hor','left','vert','top','FontSize',8.5);
axis([-100  max(max(eas),max(ea_n)) 0  max(hx/length(ea_n))+0.01]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%% WIND GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edges=[0:1:15];
[hx]=histc(Ws_n,edges);
[hs]=histc(Wss,edges);
%[Acx,lags,Bounds]=autocorr(Ws_n,150);
%[Acss,lags,Bounds]=autocorr(Wss,150);
[Acx,lags]=xcov(Ws_n,Ws_n,150,'coeff'); %% 
[Acss,lags]=xcov(Wss,Wss,150,'coeff'); %% 
figure(601)
subplot(1,2,1); set(gca,'FontSize',9.5);
bar(edges,hx/length(Ws_n),0.75,'r','LineWidth', 2)
hold on ; grid on;
plot(edges,hs/length(Wss),'+-g','LineWidth', 2)
legend('OBS.','SIM.')
xlabel('[m/s]'); ylabel('Frequency')
title('a) PDF Wind speed')
text(10,0.125,['E_{obs}=',num2str(mean(Ws_n))],'hor','left','vert','top','FontSize',8.5);
text(10,0.1,['E_{sim}=',num2str(mean(Wss))],'hor','left','vert','top','FontSize',8.5);
text(10,0.075,['\sigma_{obs}=',num2str(std(Ws_n))],'hor','left','vert','top','FontSize',8.5);
text(10,0.05,['\sigma_{sim}=',num2str(std(Wss))],'hor','left','vert','top','FontSize',8.5);
axis([-1  15 0  max(hx/length(Ws_n))+0.03]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
n1=length(Ws); fr=24;
n2=length(Wss);
%%%%%%%%%%%%%%%%%%%%
m1=floor(n1/fr); Ws=Ws';
m2=floor(n2/fr);
Wsp=reshape(Ws(1:m1*fr),fr,m1);
Wssp=reshape(Wss(1:m2*fr),fr,m2);
Xm=nanmean(Wsp');
Xs=mean(Wssp');
stdXm=nanstd(Wsp');
stdXs=std(Wssp');
t=1:24;
figure(601)
subplot(1,2,2); set(gca,'FontSize',9.5);
plot(t,Xm,'r','LineWidth', 2)
hold on ; grid on;
plot(t,Xs,'g','LineWidth', 2)
plot(t,stdXm,'^r')
plot(t,stdXs,'^g')
legend('OBS.','SIM.')
ylabel('[m/s]'); xlabel('Hour')
title('b) Wind speed daily cycle')
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%% PRESSURE GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edges=[min(Pres):5:max(Pres)];
if length(edges) < 5 
    edges = 900:10:1035; 
end
[hx]=hist(Pre_n,edges);
[hs]=hist(Pres,edges);
%[Acx,lags,Bounds]=autocorr(Pre_n,150);
%[Acss,lags,Bounds]=autocorr(Pres,150);
[Acx,lags]=xcov(Pre_n,Pre_n,150,'coeff'); %% 
[Acss,lags]=xcov(Pres,Pres,150,'coeff'); %% 
figure(701)
set(gca,'FontSize',11);
bar(edges,hx/length(Pre_n),'r','LineWidth', 1)
hold on ; grid on;
plot(edges,hs/length(Pres),'+-g','LineWidth', 2)
legend('OBS.','SIM.')
xlabel('[mbar]'); ylabel('Frequency')
title('PDF Atmospheric pressure')
text(mean(Pres),0.25,['E_{obs}=',num2str(mean(Pre_n))],'hor','left','vert','top','FontSize',8.5);
text(mean(Pres),0.20,['E_{sim}=',num2str(mean(Pres))],'hor','left','vert','top','FontSize',8.5);
text(mean(Pres),0.15,['\sigma_{obs}=',num2str(std(Pre_n))],'hor','left','vert','top','FontSize',8.5);
text(mean(Pres),0.1,['\sigma_{sim}=',num2str(std(Pres))],'hor','left','vert','top','FontSize',8.5);
axis([min(Pre)-5  max(Pre)+5  0  max(hx/length(Pre_n))+0.03]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%




%%%%%% CROSS CORRELATION GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%
for j=1:12
    [N0s(j),N1s(j),N10,N20,Xse,ttre,pvre]=sample_extreme_prop(Prsm{j},1,24);
    [N0(j),N1(j),N10P,N20P,XseP,ttreP,pvreP]=sample_extreme_prop(Prm{j},dt,24);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [EN(j),VAR,CV,Rl,SK]=sample_properties2(Nm{j},1,24,1);
    [ENs(j),VAR,CV,Rl,SK]=sample_properties2(Nsm{j},1,24,1);%%%  average Cloudiness 
end
Nd1s=N1s.*Days; Nd1=N1.*Days; %%% Rainfall days average month                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(801)
subplot(2,1,1); set(gca,'FontSize',11);
bar(mm,Nd1,'c')
hold on ; grid on ; 
plot(mm,Nd1s,'+-k','LineWidth', 1.5);
title('a) Number of wet days')
xlabel('Month'); ylabel('[#]')
subplot(2,1,2); set(gca,'FontSize',11);
bar(mm,EN,'c')
hold on ; grid on ; 
plot(mm,ENs,'+-k','LineWidth', 1.5);
title('b) Clouidiness'); xlabel('Month') ;  ylabel('[-]')
legend('OBS.','SIM.')
%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%% FAIR WEATHER %%%%%%%%%%%%%%%%%%%%
N(N>0)=NaN; Ns(Ns>0)=NaN; 
N=N+1; Ns=Ns+1; 
Rsw=Rsw.*N; Rdif=Rdif.*N;
Rdir=Rdir.*N;
Rsws=Rsws.*Ns; Rdifs=Rdifs.*Ns;
Rdirs=Rdirs.*Ns;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
n=length(Rsw); fr=24/dt;
n1=length(Rsws);
%%%%%%%%%%%%%%%%%%%%
m=floor(n/fr); m1=floor(n1/fr); 
Rswp=reshape(Rsw(1:m*fr),fr,m);
Rswsp=reshape(Rsws(1:m1*fr),fr,m1);
Xm=nanmean(Rswp');
Xs=nanmean(Rswsp');
t=[HI:23, 0:HI-1];
figure(451)
subplot(1,3,1);
set(gca,'FontSize',10);
plot(t,Xm,'r','LineWidth', 2);
hold on ; grid on;
plot(t,Xs,'b','LineWidth', 2);
title('a) Global shortwave radiation daily cycle. Clear sky')
%legend('OBS.','SIM.')
xlabel('Hour'); ylabel('[W/m^2]')
mR=max(max(Xm,Xs))+50;
axis([0 24 0 mR])
clear Tap Tasp Xm Xs stdXm stdXs
%%%%%%%%%%%%%%%%%%%%%
n=length(Rdif); fr=24/dt;
n1=length(Rdifs);
%%%%%%%%%%%%%%%%%%%%
m=floor(n/fr);m1=floor(n1/fr); 
Rdifp=reshape(Rdif(1:m*fr),fr,m);
SDp=reshape(Rdifs(1:m1*fr),fr,m1);
Xm=nanmean(Rdifp');
Xs=nanmean(SDp');
figure(451)
subplot(1,3,3);
set(gca,'FontSize',10);
plot(t,Xm,'r','LineWidth', 2);
hold on ; grid on;
plot(t,Xs,'b','LineWidth', 2);
title('c) Diffuse shortwave radiation daily cycle. Clear sky')
legend('OBS.','SIM.')
xlabel('Hour'); ylabel('[W/m^2]')
axis([0 24 0 mR])
clear Tap Tasp Xm Xs stdXm stdXs
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
n=length(Rdif); fr=24/dt;
n1=length(Rdifs);
%%%%%%%%%%%%%%%%%%%%
m=floor(n/fr); m1=floor(n1/fr);
Rdirp=reshape(Rdir(1:m*fr),fr,m);
SBp=reshape(Rdirs(1:m1*fr),fr,m1);
Xm=nanmean(Rdirp');
Xs=nanmean(SBp');
figure(451)
subplot(1,3,2);
set(gca,'FontSize',10);
plot(t,Xm,'r','LineWidth', 2);
hold on ; grid on;
plot(t,Xs,'b','LineWidth', 2);
title('b) Direct shortwave radiation daily cycle. Clear sky')
%legend('OBS.','SIM.')
xlabel('Hour'); ylabel('[W/m^2]')
axis([0 24 0 mR])
clear Tap Tasp Xm Xs stdXm stdXs
%%%%%%%%%%%%%%%%%%%%%
return 



