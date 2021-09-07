%%%%%%%%%%% Model Moments Estimation 
%%% Ref. Cowpertwait 1996 1998 2002
function[Eh,VARh,CVh,Rh,SKh,FFh,Fddh,Fwwh]=TP_OBJ(h,lan,bet,muc,eta,alp,EP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%h = %% time interval of properties [h]
%lan= %%  [1/h] %% Rate arrival storm
%bet=; %% [1/h] %%  Rate arrival cell
%muc= %% [#]  %%   per storm
%eta= %% [1/h]  %%  cell lifetime
%alp= 
%tet=  %% [mm/h] Rainfall intensity 
%%% EP1 [mm/h]  mean precipitation 1 hour aggregation 
%  tet=(EP1*eta/(lan*###*muc))^(-alp);%%  for weibull distribution   ????
tet= eta*EP/(alp*lan*muc*h); %% for gamma distribution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%DEFINITION OF MOMENTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%WEIBULL 
%Bw = 1/alp; %%- shape parameter weibull
%Aw = tet^-Bw; %% scale parameter  weibull
%mux=  Aw*gamma(1+Bw^-1); %mux=(tet^alp)*(gamma(1+alp)); 
%EX2= (Aw^2)*gamma(1+2*Bw^-1);%EX2=(tet^(2*alp))*gamma(1+2*alp); 
%EX3= (Aw^3)*gamma(1+3*Bw^-1);%EX3=(tet^(3*alp))*gamma(1+3*alp);
%%% GAMMA 
mux= alp*tet; %%%  tet*gamma(alp)/gamma(alp); 
EX2=(tet^2)*gamma(alp+2)/gamma(alp); %% alp*tet^2;%% 
EX3= (tet^3)*gamma(alp+3)/gamma(alp); %% 2/sqrt(alp); %%
%%%EXPONENTIAL 
%mux=tet; 
%EX2=tet^2; %% tet/2; 
%EX3=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% POISSON C  
%ECC_1=muc^2; 
%ECC_12=muc^3; 
%%%%%%% POISSON C-1
%muc=muc; 
%ECC_1=muc^2-1; 
%%% GEOMETRIC 
ECC_1= 2*muc*(muc-1); 
ECC_12= 6*muc*(muc-1)^2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%
M3h=6*lan*muc*EX3*(eta*h-2+eta*h*exp(-eta*h)+2*exp(-eta*h))/(eta^4)...
    +3*lan*mux*EX2*(ECC_1)*f_2(eta,h,bet)/(2*eta^4*bet*((bet^2-eta^2)^2))...
    +lan*(mux^3)*(ECC_12)*g_1(eta,h,bet)/(2*eta^4*bet*(eta^2-bet^2)*(eta-bet)*(2*bet+eta)*(bet+2*eta)); 
%%%%%%%%%%%%%%%%%%%
Eh= lan*muc*mux*h/eta;
VARh=Gxxh(0,lan,eta,h,muc,EX2,mux,bet,ECC_1);
CVh=sqrt(VARh)/Eh; 
Rh= Gxxh(1,lan,eta,h,muc,EX2,mux,bet,ECC_1)/VARh; 
SKh=M3h/((sqrt(VARh))^3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FFh=exp(-lan*h+(lan*(bet^-1))*(muc^-1)*(1-exp(-muc + muc*exp(-bet*h)))-lan*intp(h,bet,eta,muc));
FF2h=exp(-lan*(2*h)+(lan*(bet^-1))*((muc)^-1)*(1-exp(-muc +muc*exp(-bet*(2*h))))-lan*intp((2*h),bet,eta,muc));
Fddh=FF2h/FFh; 
Fwwh=(1-2*FFh+FF2h)/(1-FFh);
return 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5  
%%%%%%%%%%%%%%%%%%%
function AA=A(h,l,eta)
if l == 0
    AA=h*eta +exp(-eta*h)-1;
else
    AA=0.5*((1-exp(-eta*h))^2)*exp(-eta*h*(l-1));
end
return
%%%%%%%%%%%%%%%%%
function BB=B(h,l,bet)
if l == 0
    BB=h*bet +exp(-bet*h)-1;
else
    BB=0.5*((1-exp(-bet*h))^2)*exp(-bet*h*(l-1));
end
return
%%%%%%%%%%%%%%%%%
function Gxxhl=Gxxh(l,lan,eta,h,muc,EX2,mux,bet,ECC_1)
Gxxhl= lan*(eta^-3)*A(h,l,eta)*( 2*muc*EX2 + ((mux^2)*(bet^2)*(ECC_1))/(bet^2-eta^2))...
    -(lan*(mux^2)*B(h,l,bet)*(ECC_1))/(bet*(bet^2-eta^2)); 
return 
%%%%%%%%%%%%%%%%%%%%%%5
function  fff=f_2(eta,h,bet)
fff= -2*eta^3*bet^2*exp(-eta*h) -2*eta^3*bet^2*exp(-bet*h) ... 
    + eta^2*bet^3*exp(-2*eta*h)+ 2*eta^4*bet*exp(-eta*h)...
     + 2*eta^4*bet*exp(-bet*h) + 2*eta^3*bet^2*exp(-(eta+bet)*h)...
     -2*eta^4*bet*exp(-(eta+bet)*h) - 8*eta^3*bet^3*h...
     + 11*eta^2*bet^3 -2*eta^4*bet+ 2*eta^3*bet^2 +4*eta*bet^5*h ... 
     + 4*eta^5*bet*h -7*bet^5 -4*eta^5 +8*bet^5*exp(-eta*h) ... 
     - bet^5*exp(-2*eta*h) - 2*h*eta^3*bet^3*exp(-eta*h)...
     -12*eta^2*bet^3*exp(-eta*h) +2*h*eta*bet^5*exp(-eta*h)...
     +4*eta^5*exp(-bet*h); 
return 
%%%%%%%%%%%%%%%%%%%%%%%
function  ggg=g_1(eta,h,bet)
ggg= 12*eta^5*bet*exp(-bet*h) +9*eta^4*bet^2 +12*eta*bet^5*exp(-eta*h)...
    +9*eta^2*bet^4 +12*eta^3*bet^3*exp(-(eta+bet)*h)- eta^2*bet^4*exp(-2*eta*h)...
     -12*eta^3*bet^3*exp(-bet*h) -9*eta^5*bet - 9*eta*bet^5 -3*eta*bet^5*exp(-2*eta*h)...
     -eta^4*bet^2*exp(-2*(bet)*h) - 12*eta^3*bet^3*exp(-eta*h)...
     + 6*eta^5*bet^2*h -10*eta^3*bet^4*h+ 6*eta^2*bet^5*h -10*eta^4*bet^3*h ... 
     + 4*eta*bet^6*h -8*bet^2*eta^4*exp(-bet*h) +4*bet*eta^6*h +12*bet^3*eta^3 ... 
     - 8*bet^4*eta^2*exp(-eta*h) -6*eta^6 -6*bet^6 - 2*eta^6*exp(-2*bet*h)...
     -2*bet^6*exp(-2*eta*h) +8*eta^6*exp(-bet*h)...
     +8*bet^6*exp(-eta*h)-3*bet*eta^5*exp(-2*bet*h); 
return 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ph_t=intp(h,bet,eta,muc)
t=0:0.01:90; 
pht = (exp(-bet.*(t+h))+1 -(eta*exp(-bet*t)-bet*exp(-eta*t))./(eta-bet)).*....
    (exp((-muc*bet.*(exp(-bet*t)-exp(-eta*t))./(eta-bet))-muc.*(exp(-bet.*t))...
    +muc*exp(-bet.*(t+h))));
ph_t=trapz(t,(1-pht)); 
return
%%%%%%%%%%%%%%%5