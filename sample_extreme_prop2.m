function[Xse,ttre,pvre]=sample_extreme_prop2(V,dt,h)
%%dt =  time step registration [h]
%%% h = time step in which I want properties [h]
% V variable [] %% additional or mean 
% l lag for autocorrelation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fr=h/dt;
n=length(V);
%%%%%%%%%%%%%%%%%%%%
m=floor(n/fr);
Vrp=reshape(V(1:m*fr),fr,m);
if fr == 1
    Vrp=Vrp(find(not(isnan(Vrp))));
    X = Vrp;
else
    X=nanmean(Vrp);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xs=sort(X);
pp=(1:m)/(m+1);
ttr=(1./(1-pp))*h/8766; % ttr = tempo di ritorno legato ad ogni probabilità espresso in anni [yr] 
pvr=-log(log(1./(1-h./(ttr.*8766)))); % pvr = variabile ridotta legata ad ogni probabilità
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% estrai gli ultimi k valori 
k=19; 
Xse=Xs(end-k:end);
ppe=pp(end-k:end);
ttre=ttr(end-k:end); 
pvre=pvr(end-k:end); 
return