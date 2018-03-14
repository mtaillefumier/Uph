% make_avi_movie_example1.m
%
%  Discussion:
%
%    This program loads a series of sequentially indexed data files to 
%    create an Audio Video Interleaved (AVI) movie that can be played
%    independently of Matlab, using, for example, the XINE player
%    for Linux. We capture NUMFRAMES frames for the x-y plot of breaking
%    waves with phase changing from 0 to 2 pi.  
%
%    Before running this program close all figures and don't interfere 
%    with the generated figures until the recording process is done.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    10 January 2015
%
%  Author:
%
%    Chola Kalale
%
clear all
clc

% cd('/home/kchola/project/results')
cd ('/home/kchola/simulations/dambreak')
output = dir('out*.dat');
name = {output.name};
% sort the files

[~,index] = sort_nat(name);
output = output(index);


numfiles = length(output);

printstep = 5000;

num_frames_per_second = 10;

% load total kinetic energy file
load TotalKineticEnergy.dat;
tke = TotalKineticEnergy(:,1);
ke = TotalKineticEnergy(:,2);
pe = TotalKineticEnergy(:,3);
% load point pressure file
load PointPressure.dat;
t = PointPressure(:,1);
p = PointPressure(:,2);
% pB = PointPressure(:,5);

%experimental data
load P1data.txt;
t1 = P1data(:,1);
P1 = P1data(:,2);

% Compute Power Loss
load PowerCalc.dat;
tl = PowerCalc(:,1);
pl = PowerCalc(:,4);
plf = PowerCalc(:,5);

% tn = sqrt(0.3/9.81);
% pn = 1000*9.81*0.2;

dt = 1.0e-5;
cmap = jet(20);
fid = figure(1);
% set(fid, 'Position', [1000 1000 500 500]);
hold on
set(fid,'nextplot','replacechildren');
set(fid,'Renderer','zbuffer');
mymap = [0 0 0.5;1 0 0;0 1 0;0 0 1;1 0 1;1 1 0];  

pn = 1000.0*9.81*0.3; % pressure normalization
tn = sqrt(0.3/9.81);%time normalization
cn = sqrt(0.3*9.81);
en = 1000.0*9.81*0.3^3*(1-2.0*0.3/1.6);
for k = 1:numfiles
     f{k}= dlmread(output(k).name);
     x = f{k}(:,2);
     y = f{k}(:,3);
     z = f{k}(:,9);
     w1 =f{k}(:,4); 
     w2 = f{k}(:,5);
     z2 = f{k}(:,11);
     
     aa = sqrt(w1.*w1 + w2.*w2);
     time(k) = printstep*k*dt/tn; % update time elapsed ... include time colum in output file
                                % w = f{k}(:,8)    
    % make plot as desired
%      pause(0.0001);
     drawnow limitrate;
     h = figure(fid)
    
     a = 2; % select whether to do animation only or to plot and write each frame to avi
     b = 0.8;
     
     if(a ==1) % animation using plot
         plot(x,y,'o','MarkerFacecolor',...
         'g','MarkerEdgeColor','b','MarkerSize',0.5,...
         'LineWidth',0.5,'LineStyle','none');
     
     elseif(a==2) % animation using scatter to render the particles
         scatter(x,y,b,aa), colormap jet;%(bluewhitered(256));%colormap (mymap); 
%          set(h, 'Position', [736 278 568 410]);
         cb = colorbar('Direction','Reverse'); % colorbar('southside');
         cb.Label.String = 'Pressure [Pa]';
%          caxis([-100, 5000]);
          set(gca,'color','black');
      elseif(a==4) % animation using scatter to render the particles
         [X,Y] = meshgrid(x,y);
         quiver(x,y,w1,w2);  
     
          
      
%          s = scatter(x,y,b,z);
%          s.MarkerEdgeColor ='b'; %[0.8 0.3 0.1];
%          s.MarkerFaceColor ='g';%[0.3 0.2 0.1]; 
%          s.LineWidth = 0.6;
%          , colorbar, colormap jet;
%          c = linspace(1,10,length(x)); 
%          scatter(x,y,0.2,c,'filled');
     elseif(a ==3) % animation using plot and same image
         plot(x,y,'o','MarkerFacecolor',...
         'g','MarkerEdgeColor','b','MarkerSize',0.5,...
         'LineWidth',0.5,'LineStyle','none');
         saveas(h, sprintf('anime%4d.png', k));
     else % make plots and save frames as *.png files
         scatter(x,y,b,z), colormap jet;
         saveas(h, sprintf('anime%4d.png', k));
     end


     set(0, 'defaultTextInterpreter', 'latex');
     set(gca, 'FontName', 'Arial')
     set(gca,'fontsize',14);
     ylabel('y[m]','fontsize',26);
     xlabel('x[m]','fontsize',26);
     tm = title(['time = ', num2str(time(k)), 'seconds']);
    
%     axis tight
      axis([-0.2 2.1 -0.2 1]);
   
  end 
  


%  Tell MATLAB we have completed the movie.
  hold off
  % kinetic energy
  maxKE = max(ke);
  figure(2);
  plot(tke/tn, ke, 'r.');
%   plot(tke, (ipe-ke-pe)/(ipe/3-ipe), 'r.');
  xlabel('t[s]');
  ylabel('total K.E [J]/maxKE');
  title('Total Kinetic Energy of System');
  set(gca,'xminorgrid','on');
  set(gca,'yminorgrid','on');
%   pressure at specified point 
  figure(3);
      plot(t/tn, p/pn, 'r.', t1, P1, 'b.');
      xlabel('t');
      ylabel('Point Pressure');
      title('Pressure at a point');
      set(gca,'xminorgrid','on');
      set(gca,'yminorgrid','on');
%   axis([0 10 400 1800]);
  
% %   pressure at Upper point 
%   figure(4);
%       plot(t/tn, pB/pn, 'r.');
%       xlabel('t');
%       ylabel('Point Pressure');
%       title('Pressure at upper point');
%       set(gca,'xminorgrid','on');
%       set(gca,'yminorgrid','on');
%   pressure at point-Experiment 
  figure(5);
      plot(t1, P1, 'r.');
      xlabel('t');
      ylabel('Point Pressure');
      title('Pressure at a point ---experiment');
      set(gca,'xminorgrid','on');
      set(gca,'yminorgrid','on');
   %  boundary layer dissipation rate
  figure(6);
      semilogy(tl/tn, -pl*tn/en, 'r.');
      xlabel('t');
      ylabel('Power Loss');
      title('Energy DissipationRate in Boundary Layer');
      set(gca,'xminorgrid','on');
      set(gca,'yminorgrid','on');
  %  fluid dissipation rate
  figure(7);
      plot(tl/tn, -plf*tn/en, 'r.');
      xlabel('t');
      ylabel('Power Loss');
      title('Energy DissipationRate in Fluid Bulk');
      set(gca,'xminorgrid','on');
      set(gca,'yminorgrid','on');
      
  fprintf ( 1, '\n' );
  fprintf ( 1, 'MAKE_AVI_MOVIE_SPH\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '  The ".png files" have been created.\n' );
  