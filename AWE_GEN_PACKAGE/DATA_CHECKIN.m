function [Pr,N,Ta,Rsw,U,Ws,Pre,Rdir,Rdif,SPar]...
    = DATA_CHECKIN(D,Pr,N,Ta,Rsw,U,Ws,Pre,Rdir,Rdif,DeltaGMT,Lon,Lat,Zbas,Directory)
    if nansum(Pr) == 0
        errordlg('Precipitation "Pr" cannot be found!','Error01')
        return
    end 
    SPar.check(1:6)=0; 
    NT=length(D); 
    D=reshape(D,length(D),1);
    [Yr,Mo,Da,Hr,Mi,Se]=datevec(D);
    Datam=[Yr,Mo,Da,Hr]; %% Datam %% [Yr, MO, DA, HR]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if  nansum(N) == 0 
         hd3= helpdlg('Cloud cover "N" cannot be found! Specify the parameters for clouds.','Alert !');
         waitfor(hd3);
         [File_p,Dir_p] = uigetfile('*.mat','LOAD PARAMETERS CLOUD COVER');
         cd(Dir_p)
         load(File_p)
         cd(Directory)
         %%%%%%%%%%%%%
         SPar.check(1)=1; 
         SPar.M0=Par.M0; SPar.sigmam=Par.sigmam; SPar.rhom=Par.rhom; SPar.gam=Par.gam; 
         SPar.acloud = Par.acloud;  SPar.bcloud = Par.bcloud; SPar.EN1=Par.EN1;  SPar.Tr=Par.Tr;
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%% COMPUTE COLUD COVER --- FOR PARAMETER ESTIMATION 
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         NT=length(D); 
         i=0;  SS=zeros(1,NT); N=zeros(1,NT);
         bau = waitbar(0,'N GENERATION for PARAMETER ESTIMATION');
         SS = Pr>0;   
         [tb,t0] = ConditionSS(SS,length(SS)); % [h] Computes tb's and t0's during rainstorm events
         for i=1:NT
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             waitbar(i/NT,bau);
             %j=month(D(i));
             j=Mo(i);
             if i == 1
                 mt = 0; N(1)=0;
             else
                 mtm1 = mt;
                 [N(i),mt] = ComputeCloudFraction(Par.M0(j),mtm1,Par.sigmam(j),Par.rhom(j),...
                     i,t0(i),tb(i),Par.gam(j),SS(i),Par.acloud(:,j),Par.bcloud(:,j),N(i-1),Par.EN1(j));
                 %%%%%%%%%%%%%%%%%%%%
             end
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         end
         %%%%%%%%%%%%%%%%5
         close(bau) 
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if  nansum(Ta) == 0 
         hd3= helpdlg('Air temperature "Ta" cannot be found! Specify the parameters for air temperature.','Alert !');
         waitfor(hd3);
         [File_p,Dir_p] = uigetfile('*.mat','LOAD PARAMETERS AIR TEMPERATURE');
         cd(Dir_p)
         load(File_p)
         cd(Directory)
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%% COMPUTE AIR TEMPERATURE --- FOR PARAMETER ESTIMATION 
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         SPar.check(2)=1; 
         SPar.bTemp=Par.bTemp; SPar.dTbar=Par.dTbar; SPar.rhodT=Par.rhodT;
         SPar.sigmadT = Par.sigmadT; SPar.Ti=Par.Ti; Spar.TR2 = Par.TR2;
         %%%%%%%%%%%%%%%%%
         i=0;  Ta=zeros(1,NT); dT=zeros(1,NT);
         T_tilde=zeros(1,NT);
         bau = waitbar(0,'Ta GENERATION for PARAMETER ESTIMATION');
         for i=1:NT
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             waitbar(i/NT,bau);
             %j=month(D(i));
             j=Mo(i);
             Hro=Hr(i)+1;
             if i==1
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 [Ta(i),T_tilde(i),dT(i),qt,I2,I3,I4]=ComputeAirTemperature(Datam(i,:),DeltaGMT,...
                     Lon,Lat,N(i),Par.bTemp(j,:),0,Par.dTbar(j,Hro),Par.rhodT(j),Par.sigmadT(j,Hro),...
                     0,Par.Ti(j),0,0,0);
                 Ti = Par.Ti(j);
             else
                 if (Datam(i,4) == 0)
                     Ti = T_tilde(i-1);
                 end
                 [Ta(i),T_tilde(i),dT(i),qt,I2,I3,I4]=ComputeAirTemperature(Datam(i,:),DeltaGMT,...
                     Lon,Lat,N(i),Par.bTemp(j,:),qt,Par.dTbar(j,Hro),Par.rhodT(j),Par.sigmadT(j,Hro),...
                     dT(i-1),Ti,I2,I3,I4);
             end
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         end
         %%%%%%%%%%%%%%%%5
         close(bau)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if  nansum(U)== 0
        i=0; 
        hd3= helpdlg('Relative humidty "U" cannot be found! Specify the parameters for vapor pressure.','Alert !');
        waitfor(hd3);
        [File_p,Dir_p] = uigetfile('*.mat','LOAD PARAMETERS VAPOR PRESSURE');
        cd(Dir_p)
        load(File_p)
        cd(Directory)
        %%%%%%%%%%%%%%%%
        SPar.check(3)=1; 
        SPar.aVap=Par.aVap; SPar.dDem=Par.dDem; SPar.rhodDe=Par.rhodDe;
        SPar.sigmadDe = Par.sigmadDe;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%   VAPOR PRESSURE  GENERATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ea=zeros(1,NT);dDe=zeros(1,NT);
        esat=zeros(1,NT);U=zeros(1,NT);
            bau = waitbar(0,'ea GENERATION for PARAMETER ESTIMATION');
         for i=1:NT
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             waitbar(i/NT,bau);
             %j=month(D(i));
             j=Mo(i);
            esat(i)=611*exp(17.27*Ta(i)/(237.3+Ta(i))); %%
            if i< 3
                [ea(i),dDe(i)] = ComputeVapPressure(esat(i),Ta(i),0,0,0,Par.aVap(j,:),...
                    Par.dDem(j),Par.rhodDe(j),Par.sigmadDe(j));
            else
                [ea(i),dDe(i)] = ComputeVapPressure(esat(i),Ta(i),0,0,dDe(i-1),Par.aVap(j,:),...
                    Par.dDem(j),Par.rhodDe(j),Par.sigmadDe(j));
            end
            U(i)=ea(i)/esat(i);
        end
     close(bau)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if  nansum(Rsw) == 0
        i=0; 
        hd3= helpdlg('Global solar radiation "Rsw" cannot be found! Specify the parameters for Rsw.','Alert !');
        waitfor(hd3);
        [File_p,Dir_p] = uigetfile('*.mat','LOAD PARAMETERS SOLAR RADIATION');
        cd(Dir_p)
        load(File_p)
        cd(Directory)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SPar.check(4)=1; 
        SPar.uo=Par.uo; SPar.un=Par.un; SPar.alpha_A=Par.alpha_A;
        SPar.omega_A1 = Par.omega_A1; SPar.omega_A2 = Par.omega_A2;
        SPar.rho_g=Par.rho_g; SPar.LWP_R=Par.LWP_R; SPar.beta_A=Par.beta_A;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Tdew=zeros(1,NT);
        SB=zeros(1,NT); SD=zeros(1,NT);
        Rsw=zeros(1,NT);
        a_a=17.27; b_b=237.7;
        G_g= a_a*Ta./(b_b+Ta) + log(U);
        Tdew =b_b*G_g./(a_a-G_g);
        Tdew(isinf(Tdew))=NaN; clear G_g
        bau = waitbar(0,'Rsw GENERATION for PARAMETER ESTIMATION');
        for i=1:NT
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            waitbar(i/NT,bau);
            %j=month(D(i));%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    j=Mo(i);
            [SB(i),SD(i),SAB1,SAB2,SAD1,SAD2]=ComputeRadiationForcings(Datam(i,:),DeltaGMT,Lon,Lat,...
                Par.LWP_R(j),N(i),Zbas,Tdew(i),Par.beta_A(j),Par.alpha_A,Par.omega_A1,Par.omega_A2,Par.uo,Par.un,Par.rho_g);
            Rsw(i)=SB(i)+SD(i);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        end
        close(bau)
    end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if  nansum(Ws) == 0
       i=0; 
       hd3= helpdlg('Wind speed "Ws" cannot be found! Specify the parameters for wind speed.','Alert !');
       waitfor(hd3);
       [File_p,Dir_p] = uigetfile('*.mat','LOAD PARAMETERS WIND SPEED');
       cd(Dir_p)
       load(File_p)
       cd(Directory)
       %%%%%%%%%%%%%%%%
       SPar.check(5)=1; 
       SPar.cWind=Par.cWind; SPar.EdWs=Par.EdWs; SPar.rhodWs=Par.rhodWs;
       SPar.sigmadWs = Par.sigmadWs; SPar.skedWs = Par.skedWs;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       Ws=zeros(1,NT);dWs=zeros(1,NT);
       bau = waitbar(0,'Ws GENERATION for PARAMETER ESTIMATION');
       for i=1:NT
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           waitbar(i/NT,bau);
           if  i < 4
               [Ws(i),dWs(i)]=ComputeWindSpeed(0,0,0,Rsw(i),Par.cWind,Par.EdWs,0,Par.rhodWs,Par.sigmadWs,Par.skedWs);
           else
               [Ws(i),dWs(i)]=ComputeWindSpeed(Rsw(i-3),Rsw(i-2),Rsw(i-1),Rsw(i),Par.cWind,Par.EdWs,dWs(i-1),Par.rhodWs,Par.sigmadWs,Par.skedWs);
           end
       end
       close(bau)
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if nansum(Pre) == 0
       i=0; 
       hd3= helpdlg('Atmospheric pressure "Pre" cannot be found! Specify the parameters for pressure.','Alert !');
       waitfor(hd3);
       [File_p,Dir_p] = uigetfile('*.mat','LOAD PARAMETERS ATM. PRESSURE');
       cd(Dir_p)
       load(File_p)
       cd(Directory)
       %%%%%%%%%%%%%%%%
       SPar.check(6)=1; 
       SPar.EPre=Par.EPre; SPar.rhopre=Par.rhopre;
       SPar.sigmapre = Par.sigmapre;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       Pre=zeros(1,NT);
       bau = waitbar(0,'Pre GENERATION for PARAMETER ESTIMATION');
       for i=1:NT
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           waitbar(i/NT,bau);
           if i==1
               [Pre(i)]=ComputeAtmPressure(Par.EPre,Par.rhopre,Par.EPre,Par.sigmapre);
           else
               [Pre(i)]=ComputeAtmPressure(Par.EPre,Par.rhopre,Pre(i-1),Par.sigmapre);
           end
       end
       close(bau)
   end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if  nansum(Rdir)== 0
      Rdir=NaN*ones(1,NT);
  end
  if  nansum(Rdif) == 0
      Rdif=NaN*ones(1,NT);
  end
  %%%
  %%%
  return 
  
  
  
    
    
