function[N0,N1,N10,N20,Xse,ttre,pvre]=sample_extreme_prop(P,dt,h)
%%dt =  intervallo di registrazione piogge [h]
%%% h = intervallo sui cui voglio proprietà [h]
%%% Pr  precipitation [mm]
% l lag per l'autocorrelation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%
Xp=X>0;
%Fph =(sum(Xp)/m); 
Fph=1; 
N10= length(X(find(X>10)))/m/Fph;  %%% cell > 10 mm 
N20= length(X(find(X>20)))/m/Fph; %% cell >20 mm 
N1= length(X(find(X>1)))/m/Fph; %% cell > 1 mm 
N0= length(X(find(X>0)))/m/Fph; %% cell > 0 mm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xs=sort(X);
pp=(1:m)/(m+1);
ttr=(1./(1-pp))*h/8766; % ttr = tempo di ritorno legato ad ogni probabilità espresso in anni [yr] 
pvr=-log(log(1./(1-h./(ttr.*8766)))); % pvr = variabile ridotta legata ad ogni probabilità
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% estrai gli ultimi k valori 
k=29; 
Xse=Xs(end-k:end);
ppe=pp(end-k:end);
ttre=ttr(end-k:end); 
pvre=pvr(end-k:end); 
%%%%%%%%%%%%
%[k,alph,e]=adafg_lmom2(Xse);
%Tr=(8766/h):(8766/h):(500*8766/h); % Definizione dei Tempi di ritorno [2:1:500 anni] 
%Tr=ttre; 
%w=Tr./(Tr-1); q=1./w; 
%y=-log(-log(q));% Calcolo variabile ridotta 
%%%%%%%%%%%%%%%%%%
%Hgev= e + (alph/k)*(1-(log(w)).^k);% Altezza di pioggia per un dato tempo di ritorno con la GEV 
return