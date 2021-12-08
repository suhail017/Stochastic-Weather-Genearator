clc
close all


% Visualizing the data

Observation = Pr;% Ta %Pre % U % Rsw %Ws

Weather_gen= Prs;  %Tas %Pres %Us %Rsws %Wss

 %%  Stacked Plot
 
t = tiledlayout(3,1);

ax1 = nexttile;

plot(d,Observation,'r')

ax2 = nexttile;

plot(ds,Weather_gen,'g')

ax3 = nexttile;

plot(d,Observation,'r')

%% Basic Stats



hold on

plot(ds,Weather_gen,'g')


% Add shared title and axis labels

title(t,'Diurnal plot','FontSize',20)

xlabel(t,'Year','FontSize',20)

ylabel(t,'Precipitation (mm/h)','FontSize',20)

% Move plots closer together

linkaxes([ax1,ax2,ax3],'x');
%t.TileSpacing = 'compact';

%% Autocorrleation plot
figure()

subplot(2,1,1)

autocorr(Observation)

subplot(2,1,2)

autocorr(Weather_gen)


%% cdf plot

figure()


cdfplot(Observation)

hold on

cdfplot(Weather_gen)

legend('observation','weather generator')


%% Frequency distribuition

figure()
histogram(Observation)
hold on
histogram(Weather_gen)
xlim([0 5])
legend('Observation','Weather generator')

%% Basic Stats

median(Observation)
median(Weather_gen)
mode(Pr)
mode(Prs)
std(Pr)
std(Prs)
iqr(Pr)
iqr(Prs)
skew(Pr)
skewness(Pr)
skewness(Prs)
median(Pr)
median(Prs)


%% Misc

mon_sum = retime(TT,'monthly','sum');
stackedplot(mon_sum)
bar(mon_sum)
mon_sum1 = [mat2cell(mon_sum.Properties.RowTimes(1:12), ones(length(mon_sum.Properties.RowTimes(1:12)) ,1), 1) table2cell(mon_sum(1:12,1:2))];
Precip2 = mon_sum{1:12,2};
Precip1 = mon_sum{1:12,1};
bar(Precip1,Precip2)
bar(Precip1)
hold on
bar(Precip2)
Prcp = [Precip1; Precip2];
Prcp = [Precip1 Precip2];
bar(Prcp)
legend('Observation','Weather Generator')
title("Monthly total Precip. in Year 2000")
xlabel("Months")
ylabel("mm/h")
xlabel("Months",'FontSize',20)
ylabel("mm/h",'FontSize',20)
title("Monthly total Precip. in Year 2000",'FontSize',28)
ds_1 = month(ds(1:8784));
Prs1 =Prs(1:8784);
TR = timerange("2020-01-01","2021-01-01");
TR
Yr2021 = TT(TR,:);
Yr2021_tbl = timetable2table(Yr2021);
boxplot(Yr2021_tbl(:,2:3))
stackeplot(Yr2021)
stackedplot(Yr2021)
xlabel("year=2020",'FontSize',20)
xlabel('year=2020','FontSize',20)
stackedplot(Yr2021)
xlabel('year=2020','FontSize',20)
Precip11 = Yr2021{:,1};
Precip21 = Yr2021{:,2};
mon_sum = retime(Yr2021,'monthly','sum');
mon_sum1 = retime(Yr2021,'monthly','sum');
stackedplot(mon_sum1)
Yr2021_ms_obs=mon_sum{:,1}
Yr2021_ms_wg=mon_sum1{:,2}
Yr2021_ms_obs=mon_sum1{:,1}
bar([Yr2021_ms_obs Yr2021_ms_wg])
legend('Observation','Weather generator')
xlabel("Months","FontSize",20)
ylabel("Precip (mm/h)","FontSize",20)
title("Monthly total Precip","FontSize",20)
title("Monthly total Precip in 2020","FontSize",20)
stackeplot(Yr2021)
stackedplot(Yr2021)
TR_jan = timerange("2020","months")
bar(Yr2021.Precip2)
bar(TT.Precip1)
numTT = TT(:,vartype("numeric"));
numTT = retime(numTT,"yearly","sum");
bar(numTT.Time,numTT.Precip1)
xlabel("Year","FontSize",20)
ylabel("Total Precip by year",'FontSize',20)
title("Observational Precip","FontSize",20)
bar(numTT.Time,numTT.Precip2)
xlabel("Year","FontSize",20)
ylabel("Total Precip by year",'FontSize',20)
title("Weather generator Precip","FontSize",20)
bar(numTT.Time,numTT.Precip2)