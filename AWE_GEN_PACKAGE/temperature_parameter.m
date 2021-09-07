%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction temperature_parameter       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[b,dTbar,rhodT,sigmadT,Ti,R2] = temperature_parameter(Ta,N,Datam,DeltaGMT,Lon,Lat)
%%% TEMPERATURE PARAMETER ESTIMATION
%%%Monthly basis
%%%  INPUT %%%
%%% Ta = Temperature
%%% N = cloudiness 
%%% Datam %% [Yr, MO, DA, HR]
%%% DeltaGMT [°]
%%% Lon [°]
%%% Lat [°]
% %% OUTPUT
%%% b regression coefficients of Deterministic component of temperature first order differential equation
%%% dTbar, mean random temperature deviate
%%% rhodT, lag-1 autocorrelation random temperature deviate
%%% sigmadT %% standard deviation random temperature deviate
%%% Ti %% first step temperature [°C]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Correzione NaN
%%% Correzione NaN
Ta(isnan(Ta))=nanmean(Ta);
N(isnan(N))=nanmean(N);
%%%%%%%%%%%%%%%%%%%%%%%%%
n=length(Ta);
Ta=reshape(Ta,1,n);
N=reshape(N,1,n);
%%% Coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y= diff(Ta); %%% [°C] Hourly temperature change
%%% Coefficients
t=Datam(:,4); %% actual hour
dt=1; 
%K     = 1 - 0.65*N.^2;
K=  1 - 0.75*N.^(3.4);
Kc = (1+0.17*N.^2); %%%%  [TVA, 1972] 
sigmaSB = 5.6704e-8; % Stefan-Boltzmann constant [W/m^2 K4]  %%
qt    = sigmaSB.*Kc.*(Ta + 273.15).^4; %%  [W/m^2] q(t) estimate incoming long wave radiation
qt=qt/1000; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Computation coefficient 
Is=zeros(1,n-1); Ir=zeros(1,n-1); X1=zeros(1,n-1); X2=zeros(1,n-1); X3=zeros(1,n-1); X4=zeros(1,n-1); 
for i=2:n    
    [delta_S,h_S,zeta_S,T_sunrise,T_sunset,L_day,E0,jDay,Delta_TSL] = SetSunVariables(Datam(i,:),DeltaGMT,Lon,Lat);
    %[Is(i-1),Ir(i-1)]= Computation_int_rs(t(i),t(i-1),dt,delta_S,Lat,T_sunrise,T_sunset,Delta_TSL); 
    [Is(i-1),Ir(i-1)]= Computation_int_rs2(t(i),delta_S,Lat,T_sunrise,T_sunset,Delta_TSL); 
    X1(i-1)=Ta(i-1);
    X2(i-1)=((K(i)+K(i-1))/2)*Is(i-1);
    X3(i-1)=((K(i)+K(i-1))/2)*Ir(i-1);
    X4(i-1)= qt(i-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XX=[ones(1,n-1); X1 ; X2; X3; X4];
[a,BINT,RES,RINT,STATS]=regress(Y',XX'); 
R2=STATS(1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b(2)=-log(1+a(2));
b(1)=(b(2)/(-a(2)))*a(1);
b(3:5)=(b(2)/(-a(2)))*a(3:5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tas=zeros(1,length(t)); T_tilde=zeros(1,length(t));
%%%%%%%%%%%%%%%%%%%%%%%
Tave= mean(Ta); 
%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(t)
    if i==1
        Ti=Ta(1);
        [Tas(i),T_tilde(i),dT,qtS,I2,I3,I4]=ComputeAirTemperature(Datam(i,:),DeltaGMT,...
            Lon,Lat,N(i),b,0,0,0,0,0,Ti,0,0,0,Tave);
    else
        if Datam(i,4) == 0
            Ti = T_tilde(i-1);
        end
        [Tas(i),T_tilde(i),dT,qtS,I2,I3,I4]=ComputeAirTemperature(Datam(i,:),DeltaGMT,...
            Lon,Lat,N(i),b,qt(i-1)*1000,0,0,0,0,Ti,I2,I3,I4,Tave);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Subfunction   Compute air temperature  %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %[Ta,T_tilde,dT,qt,I2,I3,I4]=ComputeAirTemperature(Datam,DeltaGMT,...
        %Lon,Lat,N,b,qtm1,dTbar,rhodT,sigmadT,dTtm1,T_tildetm1,I2tm1,I3tm1,I4tm1)
    end
end
clear I2 I3 I4 qtS T_tilde dT 
%%%%%%%%%%%%%%%%%%%%%%%
dT = Ta- Tas;  dTbar=zeros(1,24); sigmadT=zeros(1,24); 
%%%%%%%%%%%%%%%%%%%%%%
for i=0:23
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dTh=dT(t==i);
    dTh=dTh(not(isnan(dTh)));
    dTbar(i+1)=mean(dTh); %mean random temperature deviate
    sigmadT(i+1) =std(dTh);%% standard deviation random temperature deviate
    clear  dTh 
end
dT=dT(not(isnan(dT)));
%R=autocorr(dT,10);
%rhodT=R(2); %lag-1 autocorrelation random temperature deviate
R=xcov(dT,dT,10,'coeff');
rhodT=R(12); %lag-1 autocorrelation random temperature deviate
clear R  
rhodT(rhodT>0.96)=0.96; 
Ti = nanmean(Ta); 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
