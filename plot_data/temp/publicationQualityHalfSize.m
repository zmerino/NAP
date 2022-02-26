set(0,'DefaultFigureColor','white')
fig.InvertHardcopy = 'off';
width = 6;                                                                 % Width in inches
height = 2;                                                                % Height in inches
alw = 1.5;                                                                 % AxesLineWidth 
fsz = 14;                                                                  % Fontsize 
lw = 1.25;                                                                  % LineWidth 
msz = 8;                                                                   % MarkerSize 
set(0,'defaultAxesFontSize',fsz); 
set(0,'defaultLineLineWidth',lw);       
set(0,'defaultLineMarkerSize',msz); 
set(0,'defaultAxesLineWidth',alw);
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]); 
set(0,'defaultFigurePosition', [400, 50, width*100, height*110]); 