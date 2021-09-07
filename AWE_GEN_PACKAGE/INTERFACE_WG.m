%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% INTERFACE FUNCTION  %%%%%%%%%%%%
function INTERFACE_WG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sz = get(0,'ScreenSize');
AAA = figure(1);
%set(AAA,'Position',[sz(4)*0.30   sz(4)*0.30  sz(3)*0.5  sz(4)*0.5]);
set(AAA,'Position',[sz(3)*0.21   sz(4)*0.20  sz(3)*0.62  sz(4)*0.68]);
pan_1=  uipanel('Parent',AAA,'Title','WEATHER GENERATOR','FontSize',20,'BackgroundColor','cyan',...
    'Position',[.002 .002 .996 .994]);
St1='Matlab (*.mat) input files are required for both parameter estimation and generation of synthetic series';
St2='Input file must contain variables as vectors with units and name conventions as described below (e.g., precipitation: Pr [mm])';
St3='The following variables may be included (note that precipitation is obligatory):';
St01='Department of Civil and Environmental Engineering, University of Michigan,  Ann Arbor, Michigan, USA';
St02='Department of Civil and Environmental Engineering, Massachusetts Institute of Technology, Cambridge, Massachusetts, USA';
St03='Dipartimento di Ingegneria Civile e Ambientale, Università degli Studi di Firenze,  Firenze, Italy';
%St04='Institut für Geoökologie, TU Braunschweig, Langer Kamp, Braunschweig, Germany';
u1=  uicontrol('Style','frame','Position',[0.00625*sz(3) 0.451*sz(4) 0.610*sz(3)  0.168*sz(4)]);
u1=  uicontrol('Style','text','Position',[0.01*sz(3) 0.565*sz(4) 0.60*sz(3) 0.030*sz(4) ],'String',St01,'FontSize',11.0);
u1=  uicontrol('Style','text','Position',[0.01*sz(3) 0.52*sz(4) 0.60*sz(3) 0.040*sz(4) ],'String',St02,'FontSize',11.0);
u1=  uicontrol('Style','text','Position',[0.01*sz(3) 0.475*sz(4) 0.60*sz(3) 0.030*sz(4) ],'String',St03,'FontSize',11.0);
%u1=  uicontrol('Style','text','Position',[35 340 636 25 ],'String',St04,'FontSize',11.0);
u1=  uicontrol('Style','frame','Position',[0.00859*sz(3) 0.272*sz(4) 0.606*sz(3) 0.165*sz(4) ]);
u1=  uicontrol('Style','text','Position',[0.0968*sz(3) 0.372*sz(4) 0.45*sz(3) 0.056*sz(4)],'String',St1,'FontSize',14);
u1=  uicontrol('Style','text','Position',[0.0906*sz(3) 0.321*sz(4) 0.45*sz(3) 0.05*sz(4)],'String',St2,'FontSize',11);
u1=  uicontrol('Style','text','Position',[0.0976*sz(3) 0.281*sz(4) 0.45*sz(3) 0.0325*sz(4) ],'String',St3,'FontSize',11);
u1=  uicontrol('Style','frame','Position',[0.00703*sz(3) 0.0125*sz(4) 0.607*sz(3) 0.247*sz(4)]);
u1=  uicontrol('Style','text','Position',[0.0164*sz(3) 0.202*sz(4) 0.16*sz(3) 0.041*sz(4) ],'String','INPUT DATA definition','FontSize',14);
%%%%%%%%%%%%
hpri =   uicontrol('Style','listbox','Position',[0.0218*sz(3) 0.181*sz(4) 0.207*sz(3) 0.01875*sz(4)],'String','Date "Yr" "Mo" "Da" "Hr" "Mi" or Matlab Format "D"','HandleVisibility','off');
hpri =   uicontrol('Style','listbox','Position',[0.0218*sz(3) 0.16*sz(4) 0.207*sz(3) 0.01875*sz(4)],'String','Precipitation [mm] "Pr"','HandleVisibility','off');
hni =  uicontrol('Style','listbox','Position',[0.0218*sz(3) 0.138*sz(4) 0.207*sz(3) 0.01875*sz(4)],'String','Cloud Cover [0-1] "N"');
htai =  uicontrol('Style','listbox','Position',[0.0218*sz(3)  0.1175*sz(4) 0.207*sz(3) 0.01875*sz(4)],'String','Air Temperature [°C] "Ta"');
hrswi =  uicontrol('Style','listbox','Position',[0.268*sz(3) 0.181*sz(4) 0.160*sz(3) 0.01875*sz(4)],'String','Solar Radiation: Global [W/m^2] "Rsw"');
hrdiri=   uicontrol('Style','listbox','Position',[0.268*sz(3) 0.16*sz(4) 0.160*sz(3) 0.01875*sz(4)],'String','Solar Radiation: Direct [W/m^2] "Rdir"');
hrdifi =  uicontrol('Style','listbox','Position',[0.268*sz(3) 0.138*sz(4) 0.160*sz(3) 0.01875*sz(4)],'String','Solar Radiation: Diffuse [W/m^2] "Rdif"');
hui =  uicontrol('Style','listbox','Position',[0.446*sz(3) 0.181*sz(4) 0.15*sz(3) 0.01875*sz(4)],'String','Relative Humidity [0-1] "U"');
hwsi =  uicontrol('Style','listbox','Position',[0.446*sz(3)  0.16*sz(4) 0.15*sz(3) 0.01875*sz(4)],'String','Wind Speed [m/s] "Ws"');
hprei =  uicontrol('Style','listbox','Position',[0.446*sz(3)  0.138*sz(4) 0.15*sz(3) 0.01875*sz(4)],'String','Atmospheric Pressure [mbar] "Pre"');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pbi =  uicontrol('Style','pushbutton','Position',[0.238*sz(3) 0.0237*sz(4) 0.152*sz(3) 0.0375*sz(4) ],'String','Continue','Callback', @Xinput);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitfor(pbi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BBB = figure(2);
set(BBB,'Position',[sz(3)*0.28  sz(4)*0.2 sz(3)*0.397  sz(4)*0.52]);
pan_1=  uipanel('Parent',BBB,'Title','WEATHER GENERATOR','FontSize',20,'BackgroundColor','cyan',...
    'Position',[.028 .028 .95 .95]);
u1=  uicontrol('Style','frame','Position',[0.0539*sz(3) 0.0262*sz(4) 0.275*sz(3) 0.411*sz(4) ]);
u1=  uicontrol('Style','text','Position',[0.0812*sz(3) 0.368*sz(4) 0.210*sz(3) 0.041*sz(4) ],'String','OUTPUT DATA definition','FontSize',14);
%%%%%%%%%%%%
hpri =   uicontrol('Style','listbox','Position',[0.1*sz(3) 0.35*sz(4) 0.1875*sz(3) 0.01875*sz(4)],'String','Date Matlab Format "DS"','HandleVisibility','off');
hpro =   uicontrol('Style','listbox','Position',[0.1*sz(3) 0.32*sz(4) 0.1875*sz(3) 0.01875*sz(4)],'String','Precipitation [mm] Prs') ;
hno =  uicontrol('Style','listbox','Position',[0.1*sz(3) 0.29*sz(4) 0.1875*sz(3) 0.01875*sz(4)],'String','Cloud Cover [0-1] Ns');
htao =  uicontrol('Style','listbox','Position',[0.1*sz(3) 0.26*sz(4) 0.1875*sz(3) 0.01875*sz(4)],'String','Air Temperature [°C] Tas');
hrado =  uicontrol('Style','listbox','Position',[0.1*sz(3) 0.23*sz(4) 0.1875*sz(3) 0.01875*sz(4)],'String','Solar Radiation  [W/m^2] Rsws - Rdirs - Rdifs');
heao =  uicontrol('Style','listbox','Position',[0.1*sz(3) 0.20*sz(4) 0.1875*sz(3) 0.01875*sz(4)],'String','Vapor Pressure [Pa] eas ');
hwso =  uicontrol('Style','listbox','Position',[0.1*sz(3) 0.17*sz(4) 0.1875*sz(3) 0.01875*sz(4)],'String','Wind Speed [m/s] Wss');
hpreo =  uicontrol('Style','listbox','Position',[0.1*sz(3)  0.14*sz(4) 0.1875*sz(3) 0.01875*sz(4)],'String','Atmospheric Pressure [mbar] Pres');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pbi =  uicontrol('Style','pushbutton','Position',[0.117*sz(3) 0.042*sz(4) 0.152*sz(3) 0.0375*sz(4) ],'String','Continue','Callback', @Xoutput);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitfor(pbi);
    function Xinput(hObject,ed)
        close(AAA);
    end
    function Xoutput(hObject,ed)
        close(BBB);
    end
end

