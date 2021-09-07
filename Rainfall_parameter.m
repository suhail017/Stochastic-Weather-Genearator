%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [lan,bet,muc,eta,alp,tet,Qest] = Rainfall_parameter(Prm,dt,VERB)
%%% PARAMETER ESTIMATION
%%% INPUT %%%
%%% Prm  cell(1,12) Precipitation month 
%%% dt time step 
%%%OUTPUT
%%lan=%% [1/h]
%%bet= ;%% [1/h]
%%muc=  %% per storm
%%eta=  %% [1/h]
%%alp=
%tet= %%%% [mm/h]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Optimization Procedure
%%%%%%%%%%%%%%%%%%%%%%%
qq=Prm;
i=[]; Xf=zeros(12,5);
Xo=[   0.0114	    0.0996	    5.5431	    1.4655	    1.5242];
for i=1:12 
    fd=[]; EP=[]; 
    [Eh,VARh,CVh,Rlh,SKh,FFh,Fddh,Fwwh]=sample_properties(qq{i},dt,1,1);
    fd=[fd CVh Rlh SKh FFh];
    EP=[EP ,Eh];
    [Eh,VARh,CVh,Rlh,SKh,FFh,Fddh,Fwwh]=sample_properties(qq{i},dt,6,1);
    fd=[fd CVh Rlh SKh FFh];
    EP=[EP ,Eh];
    [Eh,VARh,CVh,Rlh,SKh,FFh,Fddh,Fwwh]=sample_properties(qq{i},dt,24,1);
    fd=[fd CVh Rlh SKh FFh];
    EP=[EP ,Eh];
    [Eh,VARh,CVh,Rlh,SKh,FFh,Fddh,Fwwh]=sample_properties(qq{i},dt,72,1);
    fd=[fd CVh Rlh SKh FFh];
    EP=[EP ,Eh];
    %%% Xo=[lan,bet,muc,eta,alp];
    if VERB == 1
        options = optimset('Display','iter','MaxFunEvals',3000,'MaxIter',1500);
    else
        options = optimset('MaxFunEvals',3000,'MaxIter',1500);
    end
    %[Xf(i,:),Qest(i)] = fminsearch(@Object_function,Xo,options,fd,EP1);
    [Xf(i,:),Qest(i),exitflag]=fmincon(@Object_function,Xo,[],[],[],[],[0.0001 0.01 1 0.5 0.1 ], [ 0.05 0.99 80 30 20],[],options,fd,EP);
    Xo =Xf(i,:);
    lan(i)=Xf(i,1);% [1/h]
    bet(i)=Xf(i,2);% ;%% [1/h]
    muc(i)= Xf(i,3);% %% per storm
    eta(i)=  Xf(i,4);%%% [1/h]
    alp(i)= Xf(i,5);%
    tet(i)= eta(i)*EP(1)/(alp(i)*lan(i)*muc(i)*1); %% for gamma distribution
end
end
function Er = Object_function(Xo,fd,EP)
%%%%%%%%%%%%%%%%
lan=Xo(1);
bet=Xo(2);
muc=Xo(3);
eta=Xo(4);
alp=Xo(5);
%%%%%%%%%%%%%%
%%% EP1 mean Precipitation 1h
w=[ 1 1 1 1  1 1 1 1  1 1 1 1  1 1 1 1 ] ; %% weight vector
%%N= 18; % Number of properties fitting
h=[1 6 24 72];
for i =1:4
    [Eh,VARh,CVh,Rh,SKh,FFh,Fddh,Fwwh]=TP_OBJ(h(i),lan,bet,muc,eta,alp,EP(i));
    fm(1+ 4*(i-1))=CVh;
    fm(2+ 4*(i-1))=Rh;
    fm(4*(i-1)+3)= SKh;
    fm(4*(i-1)+4)= FFh;
end
if isnan(sum(fm))
    fm(isnan(fm))=1e+15; 
end 
Er = sum(w.*((1-fd./fm).^2 + (1-fm./fd).^2));
end