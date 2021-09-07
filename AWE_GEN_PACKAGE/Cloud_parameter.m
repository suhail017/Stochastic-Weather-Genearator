%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [M0,sigmam,rhom,gam,a,b,Tr,EN1] = Cloud_parameter(N,SS,Datam,RISP)
%%% CLOUD COVER PARAMETER ESTIMATION
%%%  INPUT %%% N cloudiness
%%% SS Rainfall detector 1 0
%%% RISP = 1 plot 0 no plot
%%% Scale Monthly
% M0 = mean fair weather cloudiness [1]
% sigmam = std fair weather cloudiness [1]
% rhom = lag-1 autocorrelation of fair weather cloudiness [1]
% chsi = gam = cloudiness decay rate [1] [1/h]
% a b = parame beta distribution cond cloud in previous hour [11,1]
% Tr = length transiction period [1] [h]
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j=0;
Mtb=ceil(max(tb)/2)-1;
kk= 0:1:Mtb; %% kk hour removal
MN=NaN*zeros(length(IS),Mtb+1);
w=NaN*zeros(length(IS),Mtb+1);
MN2=NaN*zeros(length(IS),Mtb+1);
w2=NaN*zeros(length(IS),Mtb+1);
for j=1:length(IS)
    ISN=IS{j};
    mk=ceil(length(ISN)/2)-1;
    if mk > 1
        for k=0:1:mk
            MN(j,k+1)=mean(ISN(1+k:end-k));
            w(j,k+1)=length(ISN(1+k:end-k));
            MN2(j,k+1)=(mean(ISN(1+k:1+k+1))+mean(ISN(end-k-1:end-k)))/2;
            w2(j,k+1)=2;
        end
    else
        MN(j,1)=mean(ISN);
        w(j,1)=length(ISN);
        MN2(j,1)=mean(ISN);
        w2(j,1)=length(ISN);
        %MN(j,1)=NaN;
        %w(j,1)=NaN;
    end
    clear ISN mk k
end
EN= nansum(w.*MN)./nansum(w); %% mean N function of kk
StdN= nanstd(MN); %% std N function of kk
EN2= nansum(w2.*MN2)./nansum(w2); %% mean N function of kk
StdN2= nanstd(MN2); %% std N function of kk
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% %%%% Fixed Tr  metodo sulla derivata
windowSize = 15;
ENf= filter(ones(1,windowSize)/windowSize,1,EN);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RISP == 1
    didascalia={'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'};
    j=Datam(1,2);
    figure(100)
    subplot(6,2,j); set(gca,'FontSize',8.5);
    plot(kk,EN,'r','LineWidth', 2)
    hold on ; grid on ;
    %plot(kk,StdN,'g','LineWidth', 2)
    axis([0 120 0 0.8])
    plot(kk,ENf,':m','LineWidth', 2)
    %%%%%%%%%%%%%%%%%%%%%%
    title(didascalia{j});
    if j>=11
    xlabel('Transition Period Hour')
    legend('Mean Cloudiness','Smooth Mean Cloudiness')
    end 
end
%%%%%%% %%%% Fixed Tr  metodo sulla derivata
dENf=diff(ENf);
i=15;
while (abs(dENf(i)) > 0.0005) && (i < length(dENf)-2)
    i=i+1;
end
Tr= i+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IS length(IS) > 2*Tr
% chsi = gam = cloudiness decay rate [1/h]
gam=4.61/Tr;
%%% Computation Fair weather series
j=0; Nfw=[];
for j=1:length(IS)
    ISN=IS{j};
    if length(ISN) > 2*Tr
        Nfw=[Nfw, ISN(1+Tr:end-Tr)];
    end
    clear mk ISN
end
% M0 = mean fair weather cloudiness
% sigmam = std fair weather cloudiness
% rhom = lag-1 autocorrelation of fair weather cloudiness
M0=mean(Nfw);
sigmam=std(Nfw);
%%rr=autocorr(Nfw,10);
%%rhom=rr(2); clear rr
rr=xcov(Nfw,Nfw,10,'coeff');
rhom=rr(12); clear rr 
%%%%%%%%%%%%%%
%%%%%%%% Estimation of random deviate for fair weather
i=0;
et = zeros(length(Nfw),1);
for i=2:length(Nfw)
    et(i) = ((Nfw(i)-M0) - rhom*(Nfw(i-1)-M0))/(sigmam*sqrt(1-rhom^2));
end
%%%%%%%%%%%%%%%%%%%%%%
[a,b] = Cloud_parameterII(Nfw,et,RISP);  
%%% Analytical computation Jsim
chsi=gam;
t=0:1:Tr;
Jsim = (1 - exp(-chsi*(t))).*(1 - exp(-gam*(2*Tr-t)));
%%%%%%%%%%%%% Compute observed J(t)
EN1=EN(1); 
Jo= 1 - (EN -M0)/(EN1 - M0);
EN21=EN2(1); 
Jo2= 1 - (EN2 -M0)/(EN21 - M0);
if RISP == 1
    j=Datam(1,2);
    figure(200)
    subplot(6,2,j); set(gca,'FontSize',8.5);
    %plot(kk,Jo,'ob','LineWidth', 2)
    hold on ; grid on ;
    plot(t,Jsim,'r','LineWidth', 1.8)
    axis([ 0  Tr 0 1])
    plot(kk,Jo2,'ok','LineWidth', 1)
    hold on ; grid on ;
    title(didascalia{j});
    if j>=11
    xlabel('Transition Period Hour')
    end 
    set(gca,'FontSize',5.5);
    legend(strcat('J(t)  ','T_R= ',num2str(Tr)),'Observations')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% observed mt durinf fair weather
m(1)=0; t=0;
for t=2:length(Nfw)
    m(t)= rhom*m(t-1) + et(t)*sigmam*sqrt(1 - rhom^2); %%
end
EN1=EN21;
return
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



