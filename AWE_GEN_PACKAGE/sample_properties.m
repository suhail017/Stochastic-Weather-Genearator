function[Eh,VARh,CVh,Rlh,SKh,FFh,Fddh,Fwwh]=sample_properties(P,dt,h,l)
%%dt =  intervallo di registrazione piogge [h]
%%% h = intervallo sui cui voglio proprietà [h]
% Pr  precipitation [mm]
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
    X=sum(Prp);
end
Eh=mean(X);
VARh=var(X);
CVh=std(X)/mean(X);
%R=autocorr(X,5);
%Rlh=R(l+1);
R=xcov(X,X,6,'coeff');
Rlh= R(7+l);  
SKh=skewness(X);
Xp=X>0;
dXp=diff(Xp);
wd=sum(dXp==-1);
dw=sum(dXp==1);
Xpw=(dXp+1).*(Xp(1:end-1));
ww= sum(Xpw);
dd= m-wd-dw-ww; 
FFh = 1-(sum(Xp)/m);
Fwwh = ww/(wd+ww); 
Fddh = dd/(dd+dw); 
return

