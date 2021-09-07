function[Eh,VARh,CVh,Rlh,SKh]=sample_properties2(V,dt,h,l)
%%dt =  time step registration [h]
%%% h = time step in which I want properties [h]
% V variable [] %% additional or mean 
% l lag for autocorrelation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fr=h/dt;
%%%%%%%%%%%%%%%%%%
n=length(V);
%%%%%%%%%%%%%%%%%%%%
m=floor(n/fr);
Vrp=reshape(V(1:m*fr),fr,m);
if fr == 1
    Vrp=Vrp(find(not(isnan(Vrp))));
    X = Vrp;
else
    X=nanmean(Vrp);
    X(isnan(X))=nanmean(X);
end
Eh=nanmean(X);
VARh=nanvar(X);
CVh=nanstd(X)/nanmean(X);
%R=autocorr(X,5);
%Rlh=R(l+1);
R=xcov(X,X,6,'coeff');
Rlh= R(7+l); 
SKh=skewness(X);
return

