%%%%%%%%%%%%%%%%%%%clouds parameters from web-site
function [aa, bb] = Cloud_parameterII(Nfw,et,RISP)
%%% INPUT Nfw  N fair weather
%%% et random deviate
%%% POPT graphical input
%%% RISP
%%% OUTPUT 
%%% beta parameter for random deviate aa bb 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Initialization ---
NBINS = 10;
met   = [];
stdet = [];
skewet = [];
aa = [];
bb = [];
%%aa2 = [];
%%bb2 = [];
% ============================
% --- Parameter estimation ---
% ============================
for k=0:NBINS;
    tinds = []; 
    for i=2:length(Nfw)
        if (Nfw(i-1)*10 == k)
            tinds = [tinds; i];
        end
    end
    % --- Adjustments to produce normalized series ---
    ett = zeros(length(et(tinds)),1);
    ett = et(tinds);
    % --- Transforming the series ---
    % --- To make it bounded by '0' and '1' ---
    etr1 = min(ett);
    etr2 = max(ett);
    etrange = etr2 - etr1;
    ett=ett-etr1;
    ett = ett./etrange;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Compute means, stds, skewet ---
    met = [met; mean(ett)];
    stdet = [stdet; std(ett)];
    skewet = [skewet; skewness(ett)];
    th = met(k+1)*(1-met(k+1))/(stdet(k+1)^2) - 1;
    aa = [aa; met(k+1)*th];
    bb = [bb; (1-met(k+1))*th];
    %% alternative way to calculate the random deviate 
    %ett(ett==1)=0.999999999; 
    %ett(ett==0)=0.0000000001; 
    %phat = betafit(ett);
    %aa2 = [aa2; phat(1)];
    %bb2 = [bb2; phat(2)];
    %%%%%%%%%%%%%%%%%%%%%
    A  = aa(k+1);
    B  = bb(k+1);
    %MN = (etr1*B + etr2*A)/(A+B);
    %STD = sqrt(A*B*((etr2-etr1)^2)/((A+B+1)*((A+B)^2)));
    %SKW = 2*A*B*(B-A)/(((A+B)^3)*(A+B+1)*(A+B+2)*((A*B/...
    %    (((A+B+1)*((A+B)^2))))^1.5));
    % SKW = A*((A+2)*(A+1)/((A+B+2)*(A+B+1)) - 3*A*(A+1)/...
    %      ((A+B)*(A+B+1)) + 2*A*A/((A+B)^2))/((STD^3)*(A+B));
    if RISP == 1
        figure(353) 
        subplot(6,2,k+1); set(gca,'FontSize',8);
        %%%% PLOTS SECTION 
        % ================================================
        % --- Plot empirical distribution ---
        [nn, hh] = hist(ett,NBINS);  nn = nn/length(ett);
        dx = hh(2) - hh(1);
        % --- Plotting empirical histogram ---
        % --- The sum of AREAS of all BARS must be equal to '1' ---
        bar(hh, nn/dx, 1.0); grid on; hold on; colormap cool
        ylabel('pdf')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- Plotting analytical distribution ---
        hfx = 0:0.005:1;
        fx = hfx.^(A-1).*(1-hfx).^(B-1)*(gamma(A+B))/(gamma(A)*gamma(B));
        %hfx = etr1:0.005:etr2;
        %fx = ((hfx-etr1)./(etr2-etr1)).^(A-1).*(1-(hfx-etr1)./
        %(etr2-etr1)).^(B-1)*(gamma(A+B))/(gamma(A)*gamma(B));
        %fx2= betapdf(hfx,aa2(k+1),bb2(k+1)); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        plot(hfx, fx, '-', 'LineWidth', 1.5);
        %plot(hfx, fx2, '-r', 'LineWidth', 1.5);
        plotid = ['N(t-1) = ',num2str(k)];
        title(plotid)
        if (k == 10)
            xlabel('Value of random term')
        end
        hold off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tr = 0:k;
        subplot(6,2,12); set(gca,'FontSize',8);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot(tr, aa, '-', 'LineWidth', 2.0); grid on; hold on;
        % plot(tr, bb, '-o');
        % legend('a par-r', 'b par-r');  hold off
        plot(tr, met, '-', 'LineWidth', 1.3); grid on; hold on;
        plot(tr, stdet, '--v');
        plot(tr, skewet, '-o'); hold off
        title('Statistics of transformed beta-variables')
        legend('mean', 'std', 'skew');  hold off
        xlabel('N(t-1)')
        hold off
        orient tall
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end 
 return
