%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% FUNCTION PARAMETER ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [Par] = PARAMETER_ESTIMATION(VERB,Dm,Prm,Nm,Tam,Tdewm,Rswm,Rdirm,Rdifm,eam,esatm,Rsw,Ws,Pre,DeltaGMT,Lon,Lat,Zbas,SPar)
%%% PARAMETER ESTIMATION  %%%%%
%%% DeltaGMT, Lon,Lat,Zbas
%%% Prm Nm Dm  Tam Tdewm Rswm Rdirm Rdifm  Um em esatm Ws Pre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt={ 'Ozone amount in vertical column [cm]','Total nitrogen dioxide amount [cm]',...
    'Angstrom turbidity parameter /alpha','Aerosol single-scattering albedo band #1',...
    'Aerosol single-scattering albedo band #2','Spatial average regional albedo [-]'};
DEF={'0.35','0.0002','1.3','0.92','0.84','0.15'};
r_info=inputdlg(prompt,'SET PARAMETERS FOR RADIATION ESTIMATION:',1,DEF,'on');
Par.uo=str2double(char(r_info(1,:)));
Par.un=str2double(char(r_info(2,:)));
Par.alpha_A=str2double(char(r_info(3,:)));
Par.omega_A1 =str2double(char(r_info(4,:)));
Par.omega_A2=str2double(char(r_info(5,:)));
Par.rho_g =str2double(char(r_info(6,:)));
%%%%%%%%%%%%%%%%%%
if SPar.check(4) == 1
    Ans1='blabla';
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ans1=questdlg('Do you want to specify monthly values of the Beta Angstrom Turbidity and Liquid Water Path ("No" implies their optimization)','Message','Yes',' No','Yes');
    %%%%%%%%%%%%%%%%%%
    if strcmp(Ans1,'Yes')
        prompt={ 'Specify monthly Beta Angstrom Turbidit [-]','Specify monthly Liquid Water Path for cloudy sky (N=1) [mm]'};
        DEF={' 0.010    0.017    0.020    0.027    0.035    0.046    0.060   0.044    0.023   0.010   0.004    0.010',...
            '64  47   45   41   39  45   49   48  47   46  62   74'};
        r_info=inputdlg(prompt,'SET PARAMETERS FOR RADIATION ESTIMATION:',1,DEF,'on');
        Par.beta_A=str2num(char(r_info(1,:)));
        Par.LWP_R= str2num(char(r_info(2,:)));
    else
        Ans2=questdlg('Do you want to specify monthly values of the Beta Angstrom Turbidity only ("No" implies its optimization)','Message','Yes',' No','Yes');
        if strcmp(Ans2,'Yes')
            prompt={ 'Monthly Beta Angstrom Turbidity [-]'};
            DEF={' 0.010    0.017    0.020    0.027    0.035    0.046    0.060   0.044    0.023   0.010   0.004    0.010'};
            r_info=inputdlg(prompt,'SET PARAMETERS FOR RADIATION ESTIMATION:',1,DEF,'on');
            beta_UT=str2num(char(r_info(1,:)));
        else
            beta_UT=NaN*ones(1,12);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% RAINFALL
%%%%%%%% Function Needed
%%%  sample_properties -TP_OBJ - Rainfall_parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mb =  helpdlg('Rainfall Parameter Estimation','PARAMETER ESTIMATION:');
dt = 1; %%% time step 1-Hour
[Par.lan,Par.bet,Par.muc,Par.eta,Par.alp,Par.tet,Fest] = Rainfall_parameter(Prm,dt,VERB);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear  Fest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% VARIABLES SEASONAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RUN SUI 12 MESI
didascalia={'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SPar.check(1) == 1
    Par.M0=SPar.M0; Par.sigmam=SPar.sigmam; Par.rhom=SPar.rhom; Par.gam=SPar.gam;
    Par.acloud = SPar.acloud;  Par.bcloud = SPar.bcloud; Par.EN1=SPar.EN1;  Par.Tr=SPar.Tr;
else
    for j=1:12
        %%%%%%%% SET MONTHLY VARIABLE
        SS=Prm{j}>0;
        [Yr,Mo,Da,Hr,Mi,Se]=datevec(Dm{j});
        Datam=[Yr,Mo,Da,Hr]; %% Datam %% [Yr, MO, DA, HR]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% CLOUD COVER
        %%%%%%%% Function Needed
        %%% Condition SS  Cloud_parameter Cloud_parameterII
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tit_msg= strcat('Cloud Cover Parameter Estimation ---- ',cellstr(didascalia{j}));
        mb =  helpdlg(tit_msg,'PARAMETER ESTIMATION:');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Par.M0(j),Par.sigmam(j),Par.rhom(j),Par.gam(j),Par.acloud(:,j),Par.bcloud(:,j),Par.Tr(j),Par.EN1(j)] = Cloud_parameter(Nm{j},SS,Datam,0);
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        clear tb t0 SS
    end
    %%%%%%% Correction for NaN when there are few values of cloudiness 
    for jj=1:11
        icloud1=find(isnan(Par.acloud(jj,:)));
        icloud2=find(isnan(Par.bcloud(jj,:)));
        Par.acloud(jj,icloud1)=nanmean(Par.acloud(jj,:)); 
        Par.bcloud(jj,icloud2)=nanmean(Par.bcloud(jj,:)); 
        clear icloud1 icloud2  jj 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SPar.check(2) == 1
    Par.bTemp=SPar.bTemp; Par.dTbar=SPar.dTbar; Par.rhodT=SPar.rhodT;
    Par.sigmadT = SPar.sigmadT; Par.Ti=SPar.Ti; Par.TR2 = SPar.TR2; 
else
    for j=1:12
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% AIR TEMPERATURE
        %%%%%%% Function Needed
        %%% temperature_parameter ComputeAirTemperature SetSunVariables
        %%% ComputeStochasticT ComputeDeterministicT
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        tit_msg= strcat('Air Temperature Parameter Estimation ---- ',cellstr(didascalia{j}));
        mb =  helpdlg(tit_msg,'PARAMETER ESTIMATION:');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Yr,Mo,Da,Hr,Mi,Se]=datevec(Dm{j});
        Datam=[Yr,Mo,Da,Hr]; %% Datam %% [Yr, MO, DA, HR]
        [Par.bTemp(j,:),Par.dTbar(j,:),Par.rhodT(j),Par.sigmadT(j,:),Par.Ti(j),Par.TR2(j)] =temperature_parameter(Tam{j},Nm{j},Datam,DeltaGMT,Lon,Lat);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SPar.check(4) == 1
    Par.LWP_R=SPar.LWP_R; Par.beta_A=SPar.beta_A;
else
    if strcmp(Ans1,' No')
        for j=1:12
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tit_msg= strcat('Solar Radiation Parameter Estimation ---- ',cellstr(didascalia{j}));
            mb =  helpdlg(tit_msg,'PARAMETER ESTIMATION:');
            %%%%%% SOLAR RADIATION
            %%% Function Needed
            %%% ComputeRadiationForcings SetClearSkyRadiation SetCloudySkyRadiation
            %%% SetSunVariables Radiation_parameter sample_properties2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [Yr,Mo,Da,Hr,Mi,Se]=datevec(Dm{j});
            Datam=[Yr,Mo,Da,Hr]; %% Datam %% [Yr, MO, DA, HR]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [Par.LWP_R(j),Par.beta_A(j)] = Radiation_parameter(Tdewm{j},Nm{j},Datam,DeltaGMT,Lon,Lat,Zbas,...
                Rswm{j},Par.uo,Par.alpha_A,Par.un,Par.omega_A1,Par.omega_A2,Par.rho_g,Rdirm{j},Rdifm{j},beta_UT(j),VERB); 
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555%%%%
if SPar.check(3) == 1
    Par.aVap=SPar.aVap; Par.dDem=SPar.dDem; Par.rhodDe=SPar.rhodDe;
    Par.sigmadDe = SPar.sigmadDe;
else
    for j=1:12
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%% VAPOR PRESSURE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Function Needed
        %%% vap_pre_parameter
        tit_msg= strcat('Vapor Pressure Parameter Estimation  ---- ',cellstr(didascalia{j}));
        mb =  helpdlg(tit_msg,'PARAMETER ESTIMATION:');
        [Yr,Mo,Da,Hr,Mi,Se]=datevec(Dm{j});
        Datam=[Yr,Mo,Da,Hr]; %% Datam %% [Yr, MO, DA, HR]
        %%%%% Ta(isnan(Ta))=nanmean(Ta);
        [Par.aVap(j,:),Par.dDem(j),Par.rhodDe(j),Par.sigmadDe(j)] = vap_pre_parameter(eam{j},esatm{j},Tam{j},Rswm{j},Datam);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%% WIND PARAMETER
%%% Function Needed
%%% wind_parameter
%%%%% !!!!%%%%% !!!!%%%%% !!!!
if SPar.check(5) == 1
    Par.cWind=SPar.cWind; Par.EdWs=SPar.EdWs; Par.rhodWs=SPar.rhodDWs;
    Par.sigmadWs = SPar.sigmadWs; Par.skedWs = SPar.skedWs;
else
    mb =  helpdlg('Wind Speed Parameter Estimation','PARAMETER ESTIMATION:');
    %%%% Ws_n=Ws(not(isnan(Ws)));
    [Par.cWind,Par.EdWs,Par.rhodWs,Par.sigmadWs,Par.skedWs] = wind_parameter(Ws,Rsw);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%% PRESSURE PARAMETER
if SPar.check(6) == 1
    Par.EPre=SPar.EPre; Par.rhopre=SPar.rhopre;
    Par.sigmapre = SPar.sigmapre;
else
    mb =  helpdlg('Atm Pressure Parameters Estimation','PARAMETER ESTIMATION:');
    Pre=Pre(not(isnan(Pre)));
    Par.EPre =nanmean(Pre);
    %R=autocorr(Pre,10);
    %Par.rhopre= R(2);  %lag-1 autocorrelation pressure
    R=xcov(Pre,Pre,10,'coeff');
    Par.rhopre= R(12);  %lag-1 autocorrelation pressure
    Par.sigmapre= nanstd(Pre);
end
clear R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close(mb);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


