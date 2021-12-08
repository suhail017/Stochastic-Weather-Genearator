%% 
% 
%% Weather Generator Visualization
% 
% 
% This is a script for visualizing the comparison of Weather generator and observational 
% data. We emphaisis on Precipitation and Temperature for this study. The way 
% we setup the code , we can visualize all the meteorolgical output as well. 
%% Code:

clc
close all


% Declaring the variables

Observation = Pr;% Ta %Pre % U % Rsw %Ws

Weather_gen= Prs;  %Tas %Pres %Us %Rsws %Wss

ds = datetime(DS,'ConvertFrom','datenum');

d = datetime(D,'ConvertFrom','datenum');

% Stacked Plot

t = tiledlayout(3,1);

ax1 = nexttile;

plot(d,Observation,'r')

ax2 = nexttile;

plot(ds,Weather_gen,'g')

ax3 = nexttile;

plot(d,Observation,'r')

hold on

plot(ds,Weather_gen,'g')


% Add shared title and axis labels

title(t,'Time series','FontSize',20)

xlabel(t,'Year','FontSize',20)

ylabel(t,'Precipitation (mm/h)','FontSize',20)

% Move plots closer together

linkaxes([ax1,ax2,ax3],'x');
% Autocorrelation plot

figure()
%Need to do this in raw observational data

subplot(2,1,1)

autocorr(Observation)

subplot(2,1,2)

autocorr(Weather_gen)

% CDF plot

figure()


cdfplot(Observation)

hold on

cdfplot(Weather_gen)

legend('observation','weather generator')

figure()

l_obs = log(Observation);

l_wg = log(Weather_gen);

cdfplot(l_obs)

hold on
cdfplot(l_wg)
title("Log scale CDF plot")
legend('observation','weather generator')

% Histogram

figure()
histogram(Observation)
xlabel('Precipitation','FontSize',20)
ylabel('Frequency','FontSize',20)
title("Observation")
xlim([-0.04 10])
ylim([-99 18136])


figure()
histogram(Weather_gen)
% xlim([0 5])

xlabel('Precipitation','FontSize',20)
ylabel('Frequency','FontSize',20)
title("Weather Generator")
xlim([-0.04 10])
ylim([-99 18136])
% QQplot

figure()
qqplot(l_obs,l_wg(1:length(l_obs)))
xlabel("observation","FontSize",20)
ylabel("weather generator","FontSize",20)


Year = 1:8784:length(ds);

Pr_1 = Pr(Year(1):Year(2));

Prs1 = Prs(1:length(Pr_1));
Prs1= Prs1(:);

plot(ds(1:8785),Pr_1)
hold on
plot(ds(1:8785),Prs1)
axis tight
legend('observation','weather generator')
xlabel("Year","FontSize",20)
ylabel("Precipitation","FontSize",20)
title("Year Comparison")
figure()
normplot(Observation)
hold on
normplot(Weather_gen)



% boxplot(Precip21,ds_1,"Labels",month(ds(1:8784),"shortname"))
% title("Weather Generator Precip. Data (Year = 2020)")
% 
% boxplot(Precip11,ds_1,"Labels",month(ds(1:8784),"shortname"))
% title("Observational Precip. Data (Year = 2020)")
% 
% 

% 
% figure()
% plot(ds(1:744),Pr(1:744))
% hold on
% plot(ds(1:744),Prs(1:744))
% legend("observation","Weather generator")
% title("January Data","FontSize",20)
% 
% figure()
% plot(ds(1:8784),Pr_1)
% hold on
% plot(ds(1:8784),Prs1)
% legend("observation","Weather generator")
% title("Yearly Data","FontSize",20)
% 
figure()
probplot([Observation aa(:)])
legend('Observation','Weather Generator')
title('Probability Plot')

% hold on
% probplot(Weather_gen)



% Basic Statistics



median(Observation)
median(Weather_gen)
mode(Observation)
mode(Weather_gen)
std(Observation)
std(Weather_gen)
iqr(Observation)
iqr(Weather_gen)
skewness(Observation)
skewness(Weather_gen)
prctile(Observation,[5 10 25 50 75 99],'all')
prctile(Weather_gen,[5 10 25 50 75 99],'all')
% Error Calculation

obs_data= Pr;
mod_data = Prs(1:length(Pr));
mod_data = mod_data(:);
Error = CalcPerf(obs_data,mod_data)
% Index of Agreement (Willmott et al,1981)

willmontt_index(mod_data,obs_data)