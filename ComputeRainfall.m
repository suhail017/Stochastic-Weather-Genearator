%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction   Compute Rainfall  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[PH,PP,PD]=ComputeRainfall(Ns,lan,bet,muc,eta,alp,tet)
%%%OTUPUT
%%% PH [mm/h]PP [mm]  PD [mm]
%%%%%%%%%%
%%% INPUT
N=Ns; %%% Storm number
%%lan= %% [1/h]
%%bet=% [1/h]
%%muc= %% per storm
%%eta= %% [1/h]
%%%%%%%%%%%%%%%%%%%%%%
%%% It depends from the distribution
%%% exponential tet -- gamm alp tet -- weibull alp tet
%%alp= 
%tet= ; %%%% [mm/h]
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% For the weibull they need a modification
Bw = 1/alp; %%- shape parameter weibull
Aw = tet^-Bw; %% scale parameter  weibull
%%%%%%%%%%%%%%%%%%
j=0;
Lpp=zeros(1,N); Tbs=zeros(1,N);
P=cell(N,1);
for j=1:N
    %%%%%%%%%%%%%%%%%%
    T =  poissrnd(1/lan); %% exprnd(1/lan);%%  %% %% [h]  %% Rate arrival storm
    C=   geornd(1/(1+muc)); %% poissrnd(muc); %% %%exprnd(muc);%%  %% [#]  numbers of cell  cell
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i=0;
    if C == 0
        R=0;
        L=0;
        S_T=0;
        X=0;
    else
        for i=1:C
            %R(i)=exprnd(1/fi); %% [Km] Disc Radius
            L(i)=  exprnd(1/eta); %% [h] cell lifetime
            S_T(i)= exprnd(1/bet); %% [h] Rate arrival cell
            %X(i) = wblrnd(Aw,Bw);   %% [mm/h] rainfall intensity
            X(i) = gamrnd(alp,tet); %% % [mm/h] rainfall intensity
            % X(i)=exprnd(1/tet); %% % [mm/h] rainfall intensity
        end
    end
    %%%%%%%%%%%%%%%%
    %%% Vector Precipitation Minute In Storm
    tbs=ceil(60*T); %%%% [min]  Rate arrival storm
    lpp=max(ceil(60*(S_T+L))); %%%% [min] Duration storm
    Pe=zeros(1,lpp); %% [mm/h]
    for i =1:C
        p1=zeros(1,lpp); %%
        %p1(ceil(60*S_T(i)):ceil(60*(S_T(i)+L(i))))=X(i);
        L1=length(ceil(60*S_T(i)):round(60*(S_T(i)+L(i))));
        p1(ceil(60*S_T(i)):round(60*(S_T(i)+L(i))))=X(i)*(60*L(i)/L1);
        Pe=Pe+p1; %% [mm/h]
    end
    P{j}=Pe; %% [mm/h]
    Lpp(j)=lpp/60;
    Tbs(j)=tbs/60;
    Cr(j)=C;
    contr1(j)=sum(X.*L);
    contr2(j)= sum(Pe/60);
    clear Pe p1 p2 C T R L S_T X lpp tbs DD
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CTbs=cumsum(Tbs);
Lsimu =  max(ceil(60*(CTbs+Lpp)));  %% [min] Duration simulation
PP=zeros(1,Lsimu);
j=0;
for j=1:N
    p2=zeros(1,Lsimu); %%
    in_st=ceil(60*CTbs(j));
    le_st=round(60*Lpp(j))-1;
    p2(in_st:in_st+le_st)=P{j};
    PP=PP+p2; %% [mm/h]
end
clear p2
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%errorbar
%%%%%%%%%%%%%%%%%%
n=length(PP); fr=60;
%%%%%%%%%%%%%%%%%%%% vector in hour 
m=floor(n/fr); 
PH=reshape(PP(1:m*fr),fr,m);
PH=sum(PH/60); %% [mm]  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fr=60*24;
%%%%%%%%%%%%%%%%%%%% vector in days 
m=floor(n/fr); 
PD=reshape(PP(1:m*fr),fr,m);
PD=sum(PD/60); %% [mm] 
%%% Vector in Hours
%h=1:60:length(PP);
%PH=diff(interp1(1:length(PP),cumsum(PP)/60,h)); %% [mm/h]
%%% Vectors in Days
%d=1:24:length(PH);
%PD=diff(interp1(1:length(PH),cumsum(PH),d)); %% [mm/h]
%%%%%%%%
%%%% [mm] rainfall in a minute
%length(find(PH==0)
PP=PP/60; %% [mm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear  Lpp Tbs Cr contr1 contr2
return