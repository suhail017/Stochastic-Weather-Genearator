%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Radiation_parameter         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[LWP_R,beta_A] = Radiation_parameter(Tdew,N,Datam,DeltaGMT,Lon,Lat,Zbas,Rsw,uo,...
    alpha_A,un,omega_A1,omega_A2,rho_g,Rdir,Rdif,beta_UT,VERB)
%%% RADIATION PARAMETER ESTIMATION
%%%%%% Monthly basis
%%% INPUT %%%
%%%  Rsw %% short wave total irradiance
%%%%%%  Rdir %% sky beam irradiance without terrain effects
%%%% Rdif %% sky total diffuse irradiance  without terrain effects
%%% Tdew [°C] dew point temperature
%%% N = cloudiness
%%% Datam %% [Yr, MO, DA, HR]
%%% DeltaGMT [°]
%%% Lon [°]
%%% Lat [°]
%%% Zbas [m a.s.l.]  watershed elevation
%%% beta_UT 
%%% alpha_A omega_A1 omega_A2 uo un rho_g
% %% OUTPUT
%%% LWP0 [g/m^2]
%%% beta_A []
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nansum(Rdir) == 0) ||  (nansum(Rdif) ==  0)
    ANSW = 0;
else
    ANSW =1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Adj NaN
IN=not(isnan(N)); IN=single(IN); IN(IN==0)=NaN;  
Tdew(isnan(Tdew))=nanmean(Tdew);
N(isnan(N))=nanmean(N); 
%%% Adj dimension
n=length(Tdew);
Tdew=reshape(Tdew,1,n);
N=reshape(N,1,n);
Rsw=reshape(Rsw,1,n);
Rdir=reshape(Rdir,1,n);
Rdif=reshape(Rdif,1,n);
IN=reshape(IN,1,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isnan(beta_UT)
    LWP=0; bbeta=0.05;
    if VERB == 1
        options =optimset('Display','iter','MaxFunEvals',35,'MaxIter',30,'TolFun',1 ,'TolX',0.001);
    else
        options =optimset('MaxFunEvals',35,'MaxIter',30,'TolFun',1 ,'TolX',0.001);
    end
    Xf=fminbnd(@OF_Rad,0,0.2,options,Rsw,Rdir,Rdif,Tdew,N,Datam,DeltaGMT,Lon,Lat,Zbas,LWP,bbeta,...
        alpha_A,un,uo,omega_A1,omega_A2,rho_g,ANSW,1,IN);
    beta_A=Xf(1);
else
    beta_A=beta_UT;
    LWP=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if VERB == 1
    options =optimset('Display','iter','MaxFunEvals',35,'MaxIter',30,'TolFun',1 ,'TolX',1);
else
    options =optimset('MaxFunEvals',35,'MaxIter',30,'TolFun',1 ,'TolX',1);
end
Xf=fminbnd(@OF_Rad,0,250,options,Rsw,Rdir,Rdif,Tdew,N,Datam,DeltaGMT,Lon,Lat,Zbas,LWP,beta_A,...
    alpha_A,un,uo,omega_A1,omega_A2,rho_g,ANSW,0,IN);
LWP_R =Xf(1); %[g/m^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Er = OF_Rad(Xo,Rsw,Rdir,Rdif,Tdew,N,Datam,DeltaGMT,Lon,Lat,Zbas,LWP,beta_A,alpha_A,un,uo,omega_A1,omega_A2,rho_g,ANSW,Csky,IN)
%%%%%%%%%%%%%%%%
switch  Csky 
    case 1 
    beta_A=Xo;
    case 0 
    LWP=Xo;
end
%%%%%%%%%%%
NT=length(Datam(:,1));
i=0; SB=zeros(1,NT); SD=zeros(1,NT); %
for i=1:NT
    switch  Csky
        case 1
            if N(i) == 0
                [SB(i),SD(i)]=ComputeRadiationForcings(Datam(i,:),DeltaGMT,Lon,Lat,...
                    LWP,N(i),Zbas,Tdew(i),beta_A,alpha_A,omega_A1,omega_A2,uo,un,rho_g);
            else
                SB(i)=0; SD(i)=0;
            end
        case 0
            [SB(i),SD(i)]=ComputeRadiationForcings(Datam(i,:),DeltaGMT,Lon,Lat,...
                LWP,N(i),Zbas,Tdew(i),beta_A,alpha_A,omega_A1,omega_A2,uo,un,rho_g);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SB=SB.*IN; SD=SD.*IN; %%% Erase N-NaN effect  
%%%%%%%%%%%%%%%%%%%%%%%
if  Csky == 1 
    SB=SB.*(N==0); SD=SD.*(N==0); Rsw=Rsw.*(N==0); Rdif=Rdif.*(N==0); Rdir=Rdir.*(N==0);
end
if ANSW == 1
    %%%%%
    w=[1 1 1];
    Rsw_s= SB+SD;
    %%
    [Edir]=sample_properties2(Rdir,1,24,1);
    [Edif]=sample_properties2(Rdif,1,24,1);
    [Etot]=sample_properties2(Rsw,1,24,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Edir_s]=sample_properties2(SB,1,24,1);
    [Edif_s]=sample_properties2(SD,1,24,1);
    [Etot_s]=sample_properties2(Rsw_s,1,24,1);%
    Er = w(1)*abs(Edir-Edir_s) + w(2)*abs(Edif-Edif_s) + w(3)*abs(Etot-Etot_s);
else
    %%%%%%
    Rsw_s= SB+SD;
    %%
    [Etot]=sample_properties2(Rsw,1,24,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Etot_s]=sample_properties2(Rsw_s,1,24,1);%
    Er = abs(Etot-Etot_s);
end
end