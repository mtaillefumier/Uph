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
cd ('/home/kchola/simulations/waves1')

output = dir('out*.dat');
name = {output.name};
% sort the files

[~,index] = sort_nat(name);
output = output(index);


numfiles = length(output);

printstep = 1000;

num_frames_per_second = 10;

% load total kinetic energy file
load TotalKineticEnergy.dat;
tke = TotalKineticEnergy(:,1);
ke = TotalKineticEnergy(:,2);
pe = TotalKineticEnergy(:,3);
%ipe = TotalKineticEnergy(:,4);

% load point pressure file
load PointPressure.dat;
t = PointPressure(:,1);
p = PointPressure(:,2);
% cr = PointPressure(:,4);

% 
% % load point pressure file
% load DensityError.dat;
% terr = DensityError(:,1);
% derr = DensityError(:,2);


% Compute Power Loss
load PowerCalc.dat;
tl = PowerCalc(:,1);
pl = PowerCalc(:,4);
plf = PowerCalc(:,5);
% Compute Power Loss in control volumes
load ControlVolPower.dat;
tl = ControlVolPower(:,1);
cv1 = ControlVolPower(:,2);
cv2 = ControlVolPower(:,3);
cv3 = ControlVolPower(:,4);
cv4 = ControlVolPower(:,5);

% Compute Power Loss in control volumes
load ControlVolKinetic.dat;
tk = ControlVolKinetic(:,1);
kcv1 = ControlVolKinetic(:,2);
kcv2 = ControlVolKinetic(:,3);
kcv3 = ControlVolKinetic(:,4);
kcv4 = ControlVolKinetic(:,5);

tn = sqrt(0.2/9.81);
pn = 1000*9.81*0.2;
en = 1000.0*9.81*0.2^3*(1-2.0*0.2/1.6);

dt = 1.0e-4;
cmap = jet(20);
fid = figure(1);
set(fid, 'Position', [0 500 1025 354]);
hold on
set(fid,'nextplot','replacechildren');
set(fid,'Renderer','zbuffer');
mymap = [0 0.4 0.3;1 0 0;0 1 0;0 0 1;1 0 1;1 1 0];  
for k = 1:numfiles
     f{k}= dlmread(output(k).name);
     x = f{k}(:,2);
     y = f{k}(:,3);
     z = f{k}(:,11);
     w1 =f{k}(:,4); 
     w2 = f{k}(:,5);
     
     aa = sqrt(w1.*w1 + w2.*w2);
     time(k) = printstep*k*dt; % update time elapsed ... include time colum in output file
                                % w = f{k}(:,8)    
    % make plot as desired
%      pause(0.0001);
     drawnow limitrate;
     h = figure(fid)
    
     a = 2; % select whether to do animation only or to plot and write each frame to avi
     b = 0.8;
     cn = sqrt(9.81*0.16);
%      pn = 1000.0*9.81*0.3; % pressure normalization
%      tn = sqrt(0.3/9.81);%time normalization
     
     if(a ==1) % animation using plot
         plot(x,y,'o','MarkerFacecolor',...
         'g','MarkerEdgeColor','b','MarkerSize',0.5,...
         'LineWidth',0.5,'LineStyle','none');
     
     elseif(a==2) % animation using scatter to render the particles
         scatter(x,y,b,z), colormap jet(256);%(bluewhitered(256));%colormap (mymap); 
         cb = colorbar('southoutside');
         cb.Label.String = 'fluid velocity/phase velocity u/(gh_{0})^{0.5}';
         set(gca,'color','black');

         caxis([0, 0.3]);
%           set(gcf, 'InvertHardCopy', 'off');
%           print(h, '-dpng','-r100', sprintf('anime%4d.png', k));
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
     ylabel('y[m]','fontsize',16);
     xlabel('x[m]','fontsize',16);
     tm = title(['time = ', num2str(time(k)), 'seconds']);
    
%     axis tight
      axis([-0.1 5.0 -0.1 0.4]);
   
  end 
  


%  Tell MATLAB we have completed the movie.
  hold off
  % kinetic energy
  maxKE = max(ke)
  minKE = min(ke)
  figure(2);
  plot(tke, ke, 'r.');
%   semilogy(tke, ke, 'r.');
%   plot(tke, (ipe-ke-pe)/(ipe/3-ipe), 'r.');
  xlabel('t[s]');
  ylabel('total K.E [J]/maxKE', 'fontsize',16);
  title('Total Kinetic Energy of System', 'fontsize',16);
  set(gca,'xminorgrid','on');
  set(gca,'yminorgrid','on');
%   pressure at specified point 
  tp = t/tn;
  pp = p/pn;
  figure(3);
      plot(t/tn, p/pn, 'r.');
      xlabel('t(g/H)^{1/2}', 'fontsize',8);
      ylabel('P/(\rho_{0}gH)', 'fontsize',8);
      title('Pressure at the point (1.6,0.06)', 'fontsize',16);
      set(gca,'xminorgrid','on');
      set(gca,'yminorgrid','on');
%       axis([0 8 0 0.7]);

% % uncertainty in density
%  figure(4);
%       plot(t, cr, 'r.');
%       xlabel('t[s]', 'fontsize',8);
%       ylabel('\Delta{C}_{h}(\rho)~[Kgm^{-3}]', 'fontsize',8);
%       title('mean uncertainty per particle in density due to kernel', 'fontsize',12);
%       set(gca,'xminorgrid','on');
%       set(gca,'yminorgrid','on');

 %  boundary layer dissipation rate
  figure(4);
      semilogy(tl, -pl, 'r.');
      xlabel('t');
      ylabel('Power Loss');
      title('Energy DissipationRate in Boundary Layer');
      set(gca,'xminorgrid','on');
      set(gca,'yminorgrid','on');
  %  fluid dissipation rate
  figure(5);
      plot(tl, -plf, 'r.');
      xlabel('t');
      ylabel('Power Loss');
      title('Energy DissipationRate in Fluid Bulk');
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
  %  fluid dissipation rate in control volume 1
  figure(7);
      plot(tl/tn, -cv1*tn/en, 'r.');
      xlabel('t');
      ylabel('Power Loss');
      title('Energy DissipationRate in control volume 1');
      set(gca,'xminorgrid','on');
      set(gca,'yminorgrid','on');
  %  fluid dissipation rate in control volume 2
  figure(8);
      plot(tl/tn, -cv2*tn/en, 'r.');
      xlabel('t');
      ylabel('Power Loss');
      title('Energy DissipationRate in control volume 2');
      set(gca,'xminorgrid','on');
      set(gca,'yminorgrid','on');    
  %fluid dissipation rate in control volume 3
  figure(9);
      plot(tl/tn, -cv3*tn/en, 'r.');
      xlabel('t');
      ylabel('Power Loss');
      title('Energy DissipationRate in control volume 3');
      set(gca,'xminorgrid','on');
      set(gca,'yminorgrid','on');
   %fluid dissipation rate in control volume 4
  figure(10);
      plot(tl/tn, -cv4*tn/en, 'r.');
      xlabel('t');
      ylabel('Power Loss');
      title('Energy DissipationRate in control volume 4');
      set(gca,'xminorgrid','on');
      set(gca,'yminorgrid','on'); 
  %fluid dissipation rate in control volume 4
  figure(11);
      plot(tk/tn, kcv1, 'r.',tk/tn, kcv2, 'g.',tk/tn, kcv3, 'b.',tk/tn, kcv4, 'k.');
      xlabel('t');
      ylabel('Mechanical Energy');
      title('Energy Kineticin control volume ');
      set(gca,'xminorgrid','on');
      set(gca,'yminorgrid','on');    
  
  fprintf ( 1, '\n' );
  fprintf ( 1, 'MAKE_AVI_MOVIE_SPH\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '  The ".png files" have been created.\n' );
  