%Reading timetable conversion

%TT = readtable('SEG_obs.csv');

TT = table2timetable(selinsgrovets53);
stackedplot(TT)
% make it as the table consist of Temperature and the Relative humidity

TT1 = TT(:,[4,9]);

% Make it unique rows

[C , ia]=unique(TT.valid);

TR = TT(ia,:); 

% Selecting the 53 minutes mark

TR=(minute(TR.valid)==53);

TR_Final=TT(TR,:);

% Make a interpolation for the hourly values

baal = retime(TT_final,'hourly','spline');

% Visualize

stackedplot(TR_Final)