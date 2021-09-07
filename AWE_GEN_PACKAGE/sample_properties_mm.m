function[E_max,E_min,Std_max,Std_min]=sample_properties_mm(V,dt,h)
%%dt =  time step registration [h]
%%% h = time step in which I want properties [h]
% V variable [] %% additional or mean 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fr=h/dt;
%%%%%%%%%%%%%%%%%%
n=length(V);
%%%%%%%%%%%%%%%%%%%%
m=floor(n/fr);
Vrp=reshape(V(1:m*fr),fr,m);
if fr == 1
    %%% if fr == 1 the program work differently and give max and min of the
    %%% variable 
    Vrp=Vrp(find(not(isnan(Vrp))));
    E_max=max(Vrp);
    E_min=min(Vrp);
    Std_max=NaN;
    Std_min=NaN;
else
    Xmax=nanmax(Vrp);
    Xmin=nanmin(Vrp);
    E_max=nanmean(Xmax);
    E_min=nanmean(Xmin);
    Std_max=nanstd(Xmax);
    Std_min=nanstd(Xmin);
end

return