function dop_test
% Funtion to simulate the channel model as demo file.
%
%
% Copyright 2016 Thomas Jost, German Aerospace Center (DLR)
%
%  This file is part of the satellite-to-indoor channel simulator.
% 
%     The satellite-to-indoor channel simulator is free software: 
%     you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 2 of the License, or
%     (at your option) any later version.
% 
%     The satellite-to-indoor channel simulator 
%     is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with the program.  If not, see <http://www.gnu.org/licenses/>.
% 

 c = 299792458;  % speed of light in [m/s]
 minPow = -30;   % minimum power for path to be within Npath count

  % define the size of the room (origin will be in the center of the room)
 par.room.dx = 10;  % x-dimension of the room [m]
 par.room.dy = 10;  % y-dimension of the room [m]
 par.room.dz = 3;   % z-dimension of the room [m] (height)
 par.room.winWidth = 1;  % window width [m]
 par.room.winUpper = 0.3; % distance to upper edge for the windows [m]
 par.room.winLower = 1.2; % distance to lower edge for the windows [m]
 par.room.winNums = 5;  % number of windows on wall room.winWall
 par.room.winWall = 1;  % wall number where the windows are placed
 par.room.roofWall = 5; % wall with the transmission coefficient room.Troof
 par.bandwidth = 100e6;   % bandwidth used in [Hz]
 par.del_spread = 30e-9;  % delay spread in [s]
  % define electromagnetic properties 
 par.freq = 1.5e9;  % carrier frequency in [Hz]
 par.room.Twin = 10^(-2/20);   % perpendicular transmission factor for windows [linear]
 par.room.Twall = 10^(-15/20); % perpendicular transmission factor for walls [linear] 
 par.room.Troof = 10^(-25/20); % perpendicular transmission factor for roof [linear]
  % transmitter position
 elev = 20;   % elevation in [deg] for the transmitter
 azi = 100;          % azimuth in [deg] for the transmitter
 par.trans_pos = 20e6*[cosd(azi)*cosd(elev); cosd(elev)*sind(azi); sind(elev)]; % transmitter positions
  % define multipath parameters
 par.Nscatt = 10000;   % number of scatterers to be used (optional parameter, default 10000)
  
  % receiver positions (one position per coloumn)
   % circular receiver position
 rec = [3*cosd(0:1:719).',3*sind(0:1:719).',(-par.room.dz/2+1)*ones(720,1)].';  % receiver position within the room [1m before window, 1m height]
%    % complete room for PDP verification
%  grid = 0.5;   % grid size for regular simulation
%  [x,y,z] = meshgrid((-par.room.dx/2+grid/2):grid:(par.room.dx/2-grid/2),...
%                     (-par.room.dy/2+grid/2):grid:(par.room.dx/2-grid/2),...
%                     (-par.room.dz/2+grid/2):grid:(par.room.dz/2-grid/2));  % grid in z-direction
%  rec = [x(:).';y(:).';z(:).'];  % all receiver points 
%  fprintf('Got x: %g : %g : %g\n    y: %g : %g : %g\n   z: %g : %g : %g\n',...
%           (-par.room.dx/2+grid/2),grid,(par.room.dx/2-grid/2),...
%           (-par.room.dy/2+grid/2),grid,(par.room.dx/2-grid/2),...
%           (-par.room.dz/2+grid/2),grid,(par.room.dz/2-grid/2));


  % initialisation
 walls = createRoom(par.room);  % use script to create the parameter walls for sim_chan
 tic; obj = sim_chan(par,walls);  % initialise model
 fprintf('Initialisation in %g s\n',toc);
%  plotGeometry(obj.walls,rec,obj.trans_pos); % plot for checking the geometry

 samples = 101; % number of samples per CIR
 pdpEst = zeros(1,samples);
 h = zeros(1,samples);   % dummy variable to calculate CIR
 Npath = zeros(1,size(rec,2));
 for i = 1:size(rec,2)   % go through for each receiver position
  fprintf('Calculating CIR %i/%i \n',i,size(rec,2));
  tic; cir = sim_chan(obj,rec(:,i));  % get channel output for RX position rec(:,i)
  fprintf('   needed %g s \n',toc);
  data(i).delay = cir(:,obj.info.delay)+norm(par.trans_pos(:)-rec(:,i))/c-norm(par.trans_pos(:)-rec(:,1))/c;  % save delays normalised to first receiver position
  data(i).amp = cir(:,obj.info.amp)./(exp(-1i*2*pi*norm(par.trans_pos(:)-rec(:,i))*obj.freq/c)./norm(par.trans_pos(:)-rec(:,i))); % with normalisation to LoS amplitude   % save amplitudes
  data(i).type = cir(:,obj.info.typeCol);
  fprintf(' max ampl = %g dB  otherwise = %g \n',max(20*log10(abs(data(i).amp))),max(20*log10(abs(cir(:,obj.info.amp)))));
  Npath(i) = sum(double(20*log10(abs(data(i).amp))>minPow));  % save number of paths with power above minPow
  h(i,:) = calcCIR(cir(:,obj.info.delay),cir(:,obj.info.amp),samples,obj.bandwidth);
  pdpEst = pdpEst+abs(h(i,:)).^2/size(rec,2); % overall PDP calculation
 end;
 
  % plot the geometry of the simulation
 plotGeometry(obj.walls,rec,obj.trans_pos); % plot for checking the geometry
 
  % plot CIRs over receiver positions including power in color-mode
 plotCIR3d(data,[-25,0]); title('CIR over receiver position'); 
 
  % plot PDP
 maxSample = ceil(par.del_spread*c^2/par.bandwidth);
 figure; plot(0:maxSample-1,10*log10(pdpEst(1:maxSample)),'.-b',0:maxSample-1,10*log10(calcPDP(obj.pdp,(0:maxSample-1)/obj.bandwidth,obj.bandwidth)),'.-r');
 title('Comparison to Exp. PDP'); xlabel('Delay in samples'); ylabel('PDP [dB]'); grid on;
 legend('Simulated Spatial PDP','Exp. Model');

  % plot number of paths
 plot(Npath); xlabel('Receiver Positions'); ylabel(sprintf('Number of paths above %i dB',minPow));  % plot number of paths above minPow 
 title('Number of Paths'); grid on;

return;



%% ----- Supplementary Functions ----- %%



function plotCIR3d(cir,pow)
% Function to plot the 3d CIR. cir should be a struct with
% cir(t).[delay,amp] for delay and complex amplitude at 
% time instant t. pow is optional for [minPow,maxPow] in dB to be 
% displayed.
 c0 = 299792458;  % speed of light in [m/s]
 f = figure;
 if ~exist('pow','var') || (length(pow)~=2) || ~isnumeric(pow)
  minPow = inf;  maxPow = -inf;
  for i = 1:length(cir)   % find the minimum power and the maximum power
   if any(minPow>min(10*log10(abs(cir(i).amp).^2)))
    minPow = min(10*log10(abs(cir(i).amp).^2));
   end;
   if any(maxPow<max(10*log10(abs(cir(i).amp).^2)))
    maxPow = max(10*log10(abs(cir(i).amp).^2));
   end;
  end;
 else
  pow = sort(pow,'ascend');
  [minPow,maxPow] = deal(pow(1),pow(2));
 end;
 caxis([minPow,maxPow]); 
%  cmap = flipud(colormap('Gray'));  % get colormap for gray plot 
 cmap = flipud(colormap('Hot'));  % get colormap
 for i = 1:length(cir)  % go through each time instant
  za = 10*log10(abs(cir(i).amp).^2);
  za(za<minPow) = minPow-0.1;
  za(za>maxPow) = maxPow;
  c = round((length(cmap)-1)/(maxPow-minPow)*za+1-minPow*(length(cmap)-1)/(maxPow-minPow));
  for l = 1:length(cir(i).delay)  % got through each path separately
   if (za(l)>minPow)
    plot(i,cir(i).delay(l)*1e9,'.','MarkerEdgeColor',cmap(c(l),:),'MarkerFaceColor',cmap(c(l),:),'MarkerSize',5);
    hold on;
   end;
  end;   
 end;
 colormap(cmap);  % invert colormap
 cb = colorbar; grid on;
 ylabel('delay [ns]');
 xlabel('Receiver Position Point');
 title('CIR over RX Positions');
 if ~exist('pow','var')
  pow = [0:-5:minPow,maxPow,minPow];  % powers to be shown on the colorbar
 else
  pow = [maxPow:-5:minPow];
 end;
 m = (length(cmap)-1)/(maxPow-minPow);
 for i = 1:length(pow)
  ytickVal(i) = pow(length(pow)-i+1)*m+1-minPow*m;
  yticklabel{i} = strcat(num2str(pow(length(pow)-i+1)),' [dB]');
 end;
 [ytickVal,idx] = sort(ytickVal/length(cmap),'ascend'); % sort values for setting YTick property
 yticklabel = yticklabel(idx);
 set(cb,'YTick',ytickVal);         % set new ticks
 set(cb,'YTickLabel',yticklabel);  % set new tick labels
return;


function fig = plotRoom(walls,fig)
% function to plot the room geometry.
 color = 'grcmyk';
 if ~exist('fig','var')
  fig = figure;  % create a new figure
 else
  figure(fig);   % set the figure 
 end;
 hold on; grid on;
 s = [walls(:).s];  % center of the room
 s = [(min(s(1,:))+max(s(1,:)))/2,...
      (min(s(2,:))+max(s(2,:)))/2,...
      (min(s(3,:))+max(s(3,:)))/2];
 plot3(s(1),s(2),s(3),'ob');  % plot a point in the center of the room
 for i = 1:length(walls)  % for each wall
  point = [walls(i).s(:).';
           walls(i).s(:).'+walls(i).v2(:).';
           walls(i).s(:).'+walls(i).v2(:).'+walls(i).v1(:).';
           walls(i).s(:).'+walls(i).v1(:).'];
  patch(point(:,1),point(:,2),point(:,3),color(i),'FaceAlpha',0.2);
  midPoint = walls(i).s+(walls(i).v1+walls(i).v2)/2; % middle point of plane
  quiver3(midPoint(1),midPoint(2),midPoint(3),...
          walls(i).n(1),walls(i).n(2),walls(i).n(3),color(i));
  if isfield(walls(i),'win')
   % place windows on wall
   for l = 1:length(walls(i).win)
    point = [walls(i).win(l).s(:).';
           walls(i).win(l).s(:).'+walls(i).win(l).v2(:).';
           walls(i).win(l).s(:).'+walls(i).win(l).v2(:).'+walls(i).win(l).v1(:).';
           walls(i).win(l).s(:).'+walls(i).win(l).v1(:).'];
    patch(point(:,1),point(:,2),point(:,3),'b','FaceAlpha',0.4);
   end;
  end;         
  disp(sprintf(' Wall %i : color is %c',i,color(i)));
 end;
 axis equal;
 xlabel('x-axis [m]'); ylabel('y-axis [m]'); zlabel('z-axis [m]');
return;


function fig = plotGeometry(walls,rx,tx,fig)
% Function to plot the complete geometry for visualisation
 if (exist('fig','var'))
  fig = plotRoom(walls,fig);
 else
  fig = plotRoom(walls);
 end;
  % figure is hold on in plotRoom !
 if (size(rx,2)==1) && (size(tx,2)>1)  % more tx points than rx?
  rx = repmat(rx,1,size(tx,2));
 end;
 if (size(rx,2)>1) && (size(tx,2)==1)  % more rx points than tx?
  tx = repmat(tx,1,size(rx,2));
 end;
 plot3(rx(1,:),rx(2,:),rx(3,:),'.r');  % receiver points
 t = 10*tx(:,1)/norm(tx(:,1));  % only for the first transmitter
 plot3([0,t(1)],[0,t(2)],[0,t(3)],'-b');  % plot transm. vector
%  for i = 1:size(rx,2)
%   t = 10*tx(:,i)/norm(tx(:,i));
%   plot3(rx(1,i)+[0,t(1)],rx(2,i)+[0,t(2)],rx(3,i)+[0,t(3)],'-b');
%  end;
return;


function x = doTransform(xi,k,inverse)
% function to perform a coordinate transformation with 
% affine transform in k.A and k.b calculated by coorTransform.
% Every coloum of xi may be a single point.
% If inverse does not exist or equals to false, the forward
% transform is done. If inverse==true, the inverse is done.
 if (size(xi,2)==3) && (size(xi,1)~=3)
  xi = xi.';  % shape such that each coloumn is a point
 end;
 if (~exist('inverse','var') || (inverse==false))
  x = k.A*xi+repmat(k.b,1,size(xi,2)); % do forward transform
 else
  x = pinv(k.A)*(xi-repmat(k.b,1,size(xi,2))); % take inverse transform
 end;
return;


function options = getOptionsNLS()
% Wrapper function to obtaine the options for lsqnonlin
 global nlsOpt;  % NLS optimiser used
 funOpt = func2str(nlsOpt);
 if (any(funOpt(1:3)=='my_'))
  func = @my_optimset;  % function to be called
 else
  func = @optimset;
 end;
 % func = @optimset;  % function to be called
 options = func(funOpt);
 options = func(options,'MaxFunEvals',1000);
 options = func(options,'MaxIter',1000); 
 options = func(options,'display','off'); 
return;


function h = calcCIR(delay,amp,sampN,B)
% Function to calculate the CIR.
% delay is the delay in [s].
% amp is the amplitude 
% sampeN the number of samples
% B the bandwidth
 w = generateFreqAxis(sampN).';
 h = zeros(1,sampN);
 for i = 1:length(delay)
  h = h+ifft(ifftshift(amp(i).*exp(-j*w*delay(i)*B)));
 end;
return;


function w = generateFreqAxis(N)
% function generates the frequency axis
% returns a coloumn vector.
 if (mod(N,2)==1)   % odd number case
  w = [-2*pi/N*(N-1)/2:+2*pi/N:-2*pi/N, 0:2*pi/N:2*pi/N*(N-1)/2].';
 else               % even number case
  w = [-2*pi/N*(N/2):+2*pi/N:-2*pi/N, 0:2*pi/N:2*pi/N*((N-1)/2)].';
  w = [-2*pi/N*(N/2):+2*pi/N:-2*pi/N, 0:2*pi/N:2*pi/N*((N-1)/2)].';
 end;
return;


function pdpVal = calcPDP(pdp,delay,B)
% Function to calculate the PDP according to the struct
% defined in calcEntryLoss(). delay is the delay axis in [s]
% and B the bandwidht in [Hz].
 pdpVal = pdp.gamma*exp(-pdp.beta*delay*B);
return;

