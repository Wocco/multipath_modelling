function varargout = sim_chan(obj,varargin)
% Funtion to simulate the satellite-to-indoor channel model.
%
% Initialisation call:
%   obj = sim_chan(param,walls)
% Call when running:
%   [cir,p] = sim_chan(obj,rx)
%
% Within the initialisation, the PDP is calculated using the
% entry loss. Additionally, the multipath parameters are drawn
% and initilised.
% At run time, the LoS (incl. diffractions etc. is calculated
% using deterministic methods as used for the entry loss model.
% For the multipath components (MPCs), path parameters are calculated
% deterministically based on the initialisation.
%
% Input: param.bandwidth = simulated bandwidth in [Hz]
%        param.del_spread = delay spread set by the user [s].
%        param.freq = carrier frequency [Hz]
%        param.Nscatt = number of scatterers to be used 
%                             (default 10000)
%        param.trans_pos = transmitter position in [x,y,z], local coordinate
%                          system in [m].
%        walls = structure created by createRoom(), see help createRoom.
%        obj.multi.Nscatt = Number of scatterers within the simulation
%        obj.del_spread = delay spread set by the user [s].
%        obj.trans_pos = transmitter position in [x,y,z], local coordinate
%                         system in [m].
%        obj.info.amp = index for the amplitude value. 
%        obj.info.delay = index for the delay value.
%        obj.freq = carrier frequency in [Hz]
%        obj.bandwidth = used bandwidth in [Hz]
%        obj.multi.Nscatt = number of scatterers to be placed
%        rx = receiver position as [x,y,z] in the local coordinate system in [m].
% Output: cir = CIR output as one row per ray. cir(i,obj.info.delay) is the
%                delay while cir(i,obj.info.amp) is the complex amplitude
%                non-normalised as electrical field strength of the path.
%                The delay is normalised to the LoS path delay.
%         p = coordinates of the interaction point for the ray (1 path per
%              row as in cir)
%
% Version 1.0
%  - as described in the PhD thesis of T. Jost.
% Version 1.5
%  - includes now different types of paths in terms of power.
% Version 1.6
%  - paths are now splitted into inside and outside single scatterers in
%    terms of power.
% Version 2.0
%  - path variable power, according to journal paper.
% Version 3.0
%  - published unter GPL V2.0
% Version 4.0
%  - Version with updates according to tests from Sebastien and submitted
%     to ITU
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
 
 if length(varargin)~=1
  error('Please call sim_chan() correctly!');
 end; 
  
  % -- Multipath parameter initialisation -- %
 if isstruct(varargin{1})  % in initialisation step, load parameters
  par = load('par.mat','par');  par = par.par;  % load parameter struct from file, par.mat is saved by estimPar.m
 end;
 par.winPattern = @(N)cosWin(N,0.1);  % pattern for smooth on and off of the amplitude
 par.scattAngRes = 0.1; % angular resolution for scatterer amplitude pattern in [deg]
 par.reflRes = 1e-3;    % reflection resolution for reflections on walls in [m]  
 par.gridDelayPDF = 0.5;  % grid size of receiver points to calculate the delay pdf in [m]
 global simRefl;
 simRefl = [true,true]; % simulate [first order, second order] reflections
 global interpMeth;
 interpMeth = '*spline';  % interpolation method for interp1(.)
  % -- Entry loss parameter initialisation -- %
 global angPat;  
 angPat = @calcPatCosSquare;  % angular pattern subroutine
 gridSize = 0.5;    % grid size [m] to calculate gamma0, take normally to 0.5m 
  
  % -- tests -- %
 testMPC = [false,false]; % for test output
 if (exist('testMPC','var') && testMPC(1))
  testRoutine(par);   % MPC test
 end;
 testEntr = [false,false];  % for test output
 if (exist('testEntr','var') &&testEntr(1))
  simLoSMilan(1.0); % Entry loss test, check against Milans penetration loss
 end;
  
  % -- real procedure -- %
 if isstruct(varargin{1})
   %------------------------------%
   % ------ INITIALISATION ------ %
   %------------------------------%
   % initialisation mode, calculate parameters for the different paths
  obj = createObj(obj,varargin{1});  % create obj structure 
    % -- initialise entry loss model -- %
  obj.pdp.gamma = calcPowerDelay(obj.trans_pos,obj.walls,obj.freq,obj.bandwidth,gridSize);  % calculate gamma of the PDP
  obj.pdp.beta = 2*log((1+sqrt(1+4*obj.bandwidth^2*obj.del_spread^2))/(2*obj.bandwidth*obj.del_spread)); % calculate decay factor
  obj.pdp.bandwidth = obj.bandwidth;    % save bandwith as well to make sure
  obj.pdp.type = 'exp';
   % exp type is defined as P(delay) = pdp.gamma*exp(-pdp.beta*delay*B)
   % with B the bandwidth
  if testEntr(2)
   figure; n = 0:0.01:5;  % time steps normalised to bandwidth, so in [samples]
   plot(n/obj.bandwidth*1e9,obj.pdp.gamma*exp(-n*obj.pdp.beta));
   xlabel('$\tau_n=\frac{n}{B}$ [ns]','Interpreter','latex'); ylabel('PDP non-normalised'); title('PDP for the Room');
  end;  
    % -- initialise MPCs -- %    
  rand('seed',prod(clock));  % take random seed
  obj.multi.scatt = initScatt(obj.multi.Nscatt,obj.maxRang*obj.del_spread,...
                              obj.walls,obj.trans_pos,par,obj.info.type);  % initialise scatterers
%   plotpdf(obj.multi.scatt.type,min(obj.multi.scatt.type):max(obj.multi.scatt.type));  % for testing                             
  obj.multi.refl = initWalls(obj.walls,par);  % initilialise pure reflectors at walls
  [obj.multi.lengPath,obj.multi.lagx,obj.multi.lagPr] = ...
           calcDelayPDF(obj.walls,obj.trans_pos(:),par.gridDelayPDF,obj.multi.scatt,obj.bandwidth);  % calculate the delay distribution over the whole room, NEEDS TO BE DONE!!
    % -- return output -- %       
  varargout = {obj};  % return the output
 else
   %------------------------%
   % ------ RUN MODE ------ %
   %------------------------%
   % run mode, calculate different rays
  rx = varargin{1};
  if ~isnumeric(rx) || (numel(rx)~=3)
   error('function sim_chan: parameter rx is wrong!');
  end;
    % -- entry loss model calculation -- %  
  [EEntr,distEntr,pEntr,entry] = calcPaths(obj.trans_pos(:),rx(:),obj.walls,obj.freq,obj.bandwidth); % calculate LoS paths
  typeEntr = obj.info.type.directPath*ones(length(EEntr),1);  % set direct path type for entry loss paths
  if testEntr(2)  % testing case for enty loss model, plot CIR
   figure; stem(1e9*distEntr/c,20*log10(abs(EEntr./(exp(-1i*2*pi*norm(obj.trans_pos(:)-rx(:))*obj.freq/c)./norm(obj.trans_pos(:)-rx(:))))),'b'); hold on; grid on;
   stem(0,20*log10(abs(entry./(exp(-1i*2*pi*norm(obj.trans_pos(:)-rx(:))*obj.freq/c)./norm(obj.trans_pos(:)-rx(:))))),'r');
   title('CIR, det. comp.'); xlabel('delay [ns]'); ylabel('$\left|E/E_{\mathrm{LoS}}\right|$ [dB]','Interpreter','latex');
   legend('Single components','Entry Loss');
   keyboard;
  end;
     % -- MPC calculation -- %  
  [Escat,dScatt,pScatt,typeScatt,avPowScatt] = calcScatt(obj.trans_pos(:),rx(:),obj.multi.scatt,obj.freq,obj.info.type,obj.walls); % calculate multipath originated by scattering
  [Ewall,dWall,pWall,typeWall,avPowRefl] = calcRefl(obj.trans_pos(:),rx(:),obj.freq,obj.multi.refl,obj.info.type);  % calculate wall reflections
  if (isempty(Escat) && isempty(Ewall))
   disp('No MPC paths!!!');
  end;
   % combine paths and normalise to LoS distance and PDP
   % -- generate output -- %
  distMPC = [dScatt; dWall]-norm(obj.trans_pos(:)-rx(:));  % path length calc.
%   EMPC = normToPDP([Escat; Ewall],distMPC,obj.multi,obj.pdp,length(EEntr)+length(Escat)+length(Ewall));    % amplitude
  EMPC = normToPDP([Escat; Ewall],distMPC,obj.multi,obj.pdp,[avPowScatt;avPowRefl]);    % amplitude
  cir(:,obj.info.amp) = [EEntr(:);EMPC(:)];  % set amplitude of paths (not normalised to LoS path)
  cir(:,obj.info.delay) = [distEntr(:)/c;distMPC(:)/c];  % delay for the paths (normalised to LoS path)
  cir(:,obj.info.typeCol) = [typeEntr(:);typeScatt(:);typeWall(:)];  % get the type of path
  p = [pEntr;pScatt; pWall];  % interaction points, one point per line
  varargout = {cir,p};   % and ... generate output
  if any(cir(:,obj.info.delay)<0) || any(abs(cir(:,obj.info.amp)).^2>100*obj.pdp.gamma)
   disp('Something strange, got some delay < 0 or amplitude is very large!');
   keyboard;
  end;
  if testMPC(2)  % testing case, plot CIR
   figure; stem(1e9*cir(:,obj.info.delay),20*log10(abs(cir(:,obj.info.amp)./(exp(-1i*2*pi*norm(obj.trans_pos(:)-rx(:))*obj.freq/c)./norm(obj.trans_pos(:)-rx(:))))),'b'); hold on; grid on;
   title('CIR, MPCs'); xlabel('delay [ns]'); ylabel('$\left|E/E_{\mathrm{LoS}}\right|$ [dB]','Interpreter','latex');
   keyboard;
  end;
 end; 
return;


%% ----- General initialisation of the struct variable ----- %%

function obj = createObj(par,walls)
% Function to create the structure obj from the parameter struct par.
 % multipath components parameters
 c = 299792458;  % speed of light in [m/s]
 %maxRang = 5*c;   % nominal value
 maxRang = 10*c;   % nominal value
 defScatt = 10000; % default value for the number of scatterers
 if isfield(par,'Nscatt')
  obj.multi.Nscatt = par.Nscatt;
 else
  obj.multi.Nscatt = defScatt;  % default value
 end;
 obj.walls = walls;
  % get additional range security for scatterer distribution
 lenMax = getMaxWallLength(walls)+maxRang/c;  % get maximum length among all walls
 if (1+lenMax/(par.del_spread*c)>maxRang/c)
  obj.maxRang = (1+lenMax/(par.del_spread*c))*c;       % additional multiplicator to del_spread
 else
  obj.maxRang = maxRang;  % additional multiplicator to del_spread
 end; 
  % take over values from the parameter struct
 obj.freq = par.freq;
 obj.trans_pos = par.trans_pos(:);  % always a coloumn vector
 obj.bandwidth = par.bandwidth;
 obj.del_spread = par.del_spread;
  % info structure for output
 obj.info.delay = 2;   % delay of the path
 obj.info.amp = 1;     % complex amplitude
 obj.info.typeCol = 3; % type of path, see obj.info.type.*
  % info structure for type
 obj.info.type.directPath = 1; 
 obj.info.type.firstReflection = 2;
 obj.info.type.secondReflection = 3;
 obj.info.type.ins_fixedScatt = 4;
 obj.info.type.ins_reflScatt = 5;
 obj.info.type.ins_scattRefl = 6;
 obj.info.type.ins_ellipScatt = 7;
 obj.info.type.out_fixedScatt = 8;
 obj.info.type.out_ellipScatt = 9;
 obj.info.type.description{obj.info.type.directPath} = 'Direct components';
 obj.info.type.description{obj.info.type.firstReflection} = 'First Order Reflection';
 obj.info.type.description{obj.info.type.secondReflection} = 'Second Order Reflection';
 obj.info.type.description{obj.info.type.ins_fixedScatt} = 'Inside Fixed Scatterer';
 obj.info.type.description{obj.info.type.ins_reflScatt} = 'Inside Reflection Scatterer';
 obj.info.type.description{obj.info.type.ins_scattRefl} = 'Inside Scatterer Reflection';
 obj.info.type.description{obj.info.type.ins_ellipScatt} = 'Inside Elliptical Scatterer';
 obj.info.type.description{obj.info.type.out_fixedScatt} = 'Outside Fixed Scatterer';
 obj.info.type.description{obj.info.type.out_ellipScatt} = 'Outside Elliptical Scatterer';
return;


function len = getMaxWallLength(walls)
% Function to get the length of the longest wall out
% of all walls defining the room for simulation. Value is
% used to spread out scatterers as maxRang.
 len = 0;
 for i = 1:length(walls)
  l = max([norm(walls(i).v1),norm(walls(i).v2)]);
  if l>len
   len = l; % get new maximum length
  end;
 end;
return;


%% ----- Entry Loss Subroutines (Initialisation & Runtime) ----- %%

function gamma0 = calcPowerDelay(tx,walls,freq,bandwidth,grid)
% Function to calculate the value at tau=0 for the PDP.
% Take care, a normalisation to free space loss is not done, so 
% gamma0 corresponds directly to \bar{\gamma} in the entry loss paper.
%  c = 299792458;  % speed of light in [m/s]
 normLoS = false; % set normalisation to LoS to false for all subroutines
 if ~exist('grid','var')
  grid = 1;   % default grid size for RX receiver points
 end;
 s = [walls.s];  % all points to describe the different walls
 [x,y,z] = meshgrid((min(s(1,:))+grid/2):grid:(max(s(1,:))-grid/2),...
                    (min(s(2,:))+grid/2):grid:(max(s(2,:))-grid/2),...
                    (min(s(3,:))+grid/2):grid:(max(s(3,:))-grid/2));
%  fprintf('Got x: %g : %g : %g\n    y: %g : %g : %g\n   z: %g : %g : %g\n',...
%           (min(s(1,:))+grid/2),grid,(max(s(1,:))-grid/2),...
%           (min(s(2,:))+grid/2),grid,(max(s(2,:))-grid/2),...
%           (min(s(3,:))+grid/2),grid,(max(s(3,:))-grid/2));
 rec = [x(:).';y(:).';z(:).'];  % all receiver points
%  plot3(rec(1,:),rec(2,:),rec(3,:),'.');  % plot for checking
 [E] = calcField(tx,rec,walls,freq,bandwidth,normLoS);  % calculate field components at each rec position
 gamma0 = mean(abs(E).^2);  % calculate mean entry loss power  
return;


function [E,dist,p,entry] = calcPaths(tx,rx,walls,freq,bandwidth)
% Function to calculate the amplitude as E, the delay*c as dist
% and the interaction point as p for each different elements 
% of the room defined by walls. entry is the entry loss model.
 normLos = false;   % no normalisation to the LoS field strength
 E = []; dist = []; p = zeros(0,3,size(rx,2)); entry = 0; % init values
 for i = 1:length(walls)   % loop over the walls
  [Esum,amp,d,xp] = calcFieldWall(tx,rx,walls(i),freq,bandwidth,normLos);  % calculate field strength for components
  entry = entry+Esum;
  if ~isempty(amp)    % check if new components are available
   E(end+1:end+size(amp,1),:) = amp;   % add new components in amplitude by adding new lines, [NrPaths x NrRx]
   dist(end+1:end+size(d,1),:) = d;    % add new components in prop. distance by adding new lines, [NrPaths x NrRx]
   p(end+1:end+size(xp,1),:,:) = xp;    % add new components in prop. distance by adding new lines, [NrPaths x 3 x NrRx]
  end;
 end;
return;



%% ----- Initialisation Multipath Components ----- %%

function scatt = initScatt(Nscatt,maxRang,walls,tx_pos,par,typeStr)
% Function to initilialise all scatterers. The origin provided by the
% room (described by walls) is assumed to be in the centre.
% Input : Nscatt = number of scatterers
%         maxRang = maximum distance from the center of the room to
%                   scattering points in [m] (not in use any more !!!)
%         walls = walls for the room definition
%         tx_pos = transmitter position
%         par = parameter struct for the multipath components
%         typeStr = information struct with type of scatterer description
%                     (see also createObj())
 epsPol = 1e-4;  % minimum distance of a pole to the unit circle
  % calculate center of room
 s = [walls(:).s];  % center of the room
 midRoom = [(min(s(1,:))+max(s(1,:)))/2,...
            (min(s(2,:))+max(s(2,:)))/2,...
            (min(s(3,:))+max(s(3,:)))/2].'; % coloumn vector [x,y,z].'
  % draw scatterer location
%  el = uniRnd(Nscatt,[0,90]);  % elevation angle in degree
%  az = uniRnd(Nscatt,[0,360]);   % azimuth angle in degree
  % scatterer position, relevant for the type and the amplitude pattern,
  % not for the delay calculation, there scatt.pDelay will be set.
    % use uniform range with angles
%  scatt.pAmp = repmat(uniRnd(Nscatt,[0,maxRang]),1,3).*[cosd(el).*cosd(az), cosd(el).*sind(az), sind(el)]...
%                +repmat(midRoom.',Nscatt,1);  % normed scatterer position, per row [Nscatt x 3], from center of room, here using uniform range
    % use range estimate from data
%   scatt.pAmp = repmat(drawRnd([Nscatt,1],par.Scatt.point.absPDF,par.Scatt.point.phat),1,3).*[cosd(el).*cosd(az), cosd(el).*sind(az), sind(el)]...
%                +repmat(midRoom.',Nscatt,1);  % normed scatterer position, per row [Nscatt x 3], from center of room, here using range estimate from data
    % use a box around the center
      % only positive z-axis
%   scatt.pAmp = [reshape(uniRnd(Nscatt*2,[-maxRang/2,maxRang/2]),Nscatt,2),uniRnd(Nscatt,[0,maxRang/2])]...
%                +repmat(midRoom.',Nscatt,1);  % normed scatterer position, per row [Nscatt x 3], from center of room, here using range estimate from data
      % positive and negative z-axis
  scatt.pAmp = reshape(uniRnd(Nscatt*3,[-maxRang/2,maxRang/2]),Nscatt,3)...
               +repmat(midRoom.',Nscatt,1);  % normed scatterer position, per row [Nscatt x 3], from center of room, here using range estimate from data
% -- for testing
%  for i = 1:Nscatt
%   dist(i) = norm(trans_pos-scatt.p(i,:).')+norm(scatt.p(i,:))-norm(trans_pos);  % additional length to direct
%  end;
%  plotRoom(walls); plot3(scatt.pAmp(:,1),scatt.pAmp(:,2),scatt.pAmp(:,3),'.'); % plot scatter points
%  figure; plot(dist);   % plot additional distance to LoS path
  % draw opening angle and direction of scattering amplitude
 scatt.direc = uniRnd(Nscatt,[0,360]);  % draw direction of visibility region in [deg]
  % draw the opening angle according to scatterer inside or outside
 [within] = checkWinthinRoom(scatt.pAmp,walls);   % check which one is within the room
 scatt.open(within) = drawRnd([sum(double(within)),1],par.Scatt.open.inside.pdf,par.Scatt.open.inside.Par);  % draw opening angle in [deg] for inside scatterers
 scatt.open(~within) = drawRnd([sum(double(~within)),1],par.Scatt.open.outside.pdf,par.Scatt.open.outside.Par);  % draw opening angle in [deg] for inside scatterers 
  % version with only one opening distribution
%  scatt.open = drawRnd([Nscatt,1],par.Scatt.openPDF,par.Scatt.openPar);  % draw opening angle in [deg]
  % now draw the normalised amplitude variation over angle
 arOrd = drawRnd([Nscatt,1],par.Scatt.ar.order.pdf,par.Scatt.ar.order.p)+1;  % draw the order for each scatterer
 scatt.angLimits = zeros(Nscatt,2);  % angular limits for the scatterer
 for i = 1:Nscatt
  pol = (1-drawRnd([1,arOrd(i)],par.Scatt.ar.poles.absPdf,par.Scatt.ar.poles.absPar)).*...
           exp(1i*drawRnd([1,arOrd(i)],par.Scatt.ar.poles.phasePDF,par.Scatt.ar.poles.phasePar)); % draw complex poles
  while any(abs(pol)>=(1-epsPol))  % check for stability
   % disp('Pole outside unit circle or close to it, taking a new one!');
   idx = find(abs(pol)>=(1-epsPol));
   pol(idx) = (1-drawRnd([1,length(idx)],par.Scatt.ar.poles.absPdf,par.Scatt.ar.poles.absPar)).*...
               exp(1i*drawRnd([1,length(idx)],par.Scatt.ar.poles.phasePDF,par.Scatt.ar.poles.phasePar));  % draw new pole for stability reasons
  end;
  arCoef = poly(pol);  % calculate AR coefficients
  [cyy,lagCyy] = calcACF(30,arCoef,1);  % calculate covariance of the AR process assuming unit driving variance
  [win,lag] = par.winPattern(ceil(scatt.open(i)/par.scattAngRes)); % calculate window function for the pattern
  scatt.pat(i).K = drawRnd([1,1],par.Scatt.rice.pdf,par.Scatt.rice.par);  % draw Rice factor for scatterer in linear domain
  scatt.pat(i).ang = lag*par.scattAngRes+scatt.direc(i)-scatt.open(i)/2;   % angular axis for pattern including direction and opening
  scatt.pat(i).amp = win(:).*(simAr(arCoef,1,length(win))+sqrt(scatt.pat(i).K*cyy(1))*exp(1i*2*pi*rand(1)))/sqrt(cyy(1)*(1+scatt.pat(i).K)); % calculate normalised amplitdue with AR process, window, constant and Rice factor
  if (any(abs(scatt.pat(i).amp)>100))  % verification, should not happen
   disp('Big amplitude, please verify!');
   keyboard;
  end;
%   figure; plot(scatt.pat(i).ang,real(scatt.pat(i).amp),'.-b',scatt.pat(i).ang,imag(scatt.pat(i).amp),'.-r',scatt.pat(i).ang,abs(scatt.pat(i).amp),'g'); legend('Real Part of norm Ampl','Imag Part of norm Ampl','Abs of Ampl'); grid on;
%   figure; plot(real(pol),imag(pol),'x',cos(0:0.1:2*pi+0.1),sin(0:0.1:2*pi+0.1),'-r'); grid on; axis equal; title('Pole distribution');
%   figure; plot(lagCyy,abs(cyy/cyy(1)),'.-'); grid on; title('Covariance function normalised');
  scatt.angLimits(i,:) = [scatt.pat(i).ang(1),scatt.pat(i).ang(end)];  % save angular limits for the scatterer in [deg]
  if (scatt.angLimits(i,1)<0)
   scatt.angLimits(i,:) = scatt.angLimits(i,:)+360;   % take it always to positive angular domain
   scatt.pat(i).ang = scatt.pat(i).ang+360;    % avoid negative angles also for the axis
  end;
 end;
  % at this point, we got the amplitude and the position of the scatterer,
  % now we need to determine the type and additional parameters belonging
  % to the type of scatterer. 
 % we take the variable "within" now again into consideration ...
  % declare types as 1,2,3 =  inside as [single scatterer, moving, scatterer, refl+scatterer]
  %                  4,5 = outside as [single scatterer, moving scatterer]
 scatt.info.inside.singleScatt = typeStr.ins_fixedScatt;
 scatt.info.inside.movScatt = typeStr.ins_ellipScatt;
 scatt.info.inside.reflScatt = typeStr.ins_reflScatt;
 scatt.info.inside.scattRefl = typeStr.ins_scattRefl;
 scatt.info.outside.singleScatt = typeStr.out_fixedScatt;
 scatt.info.outside.movScatt = typeStr.out_ellipScatt;
  % now draw the type of scatterer
 scatt.type(within) = discrRnd(sum(double(within)),1,...
       [par.prob.inside.singScatt,par.prob.inside.movScatt,par.prob.inside.reflScatt,par.prob.inside.scattRefl],...
       [scatt.info.inside.singleScatt,scatt.info.inside.movScatt,scatt.info.inside.reflScatt,scatt.info.inside.scattRefl]); % draw type for inside
 scatt.type(~within) = discrRnd(sum(double(~within)),1,...
       [par.prob.outside.singScatt,par.prob.outside.movScatt],...
       [scatt.info.outside.singleScatt,scatt.info.outside.movScatt]); % draw type for outside
 scatt.pDelay = scatt.pAmp;  % set first position of scatterer for delay and amplitude to the same point
 scatt.addDelay = zeros(Nscatt,1);  % set an additional delay variable to zero first for all scatterers
  % now the type of the scatterer is fixed, so we can draw the delay
  % independent power for each one.
 scatt.av_pow(scatt.type==4) = drawRnd([1,sum(scatt.type==4)],...
                                        par.Scatt.power.pdf,par.Scatt.power.single_inside_par);  % draw for fixed scatterer inside
 scatt.av_pow(scatt.type==8) = drawRnd([1,sum(scatt.type==8)],...
                                        par.Scatt.power.pdf,par.Scatt.power.single_outside_par);  % draw for fixed scatterer outside
 scatt.av_pow((scatt.type==5)|(scatt.type==6)) = drawRnd([1,sum((scatt.type==5)|(scatt.type==6))],...
                                                   par.Scatt.power.pdf,par.Scatt.power.multiWall_par);  % draw for scatterer plus reflection at wall
 scatt.av_pow((scatt.type==7)|(scatt.type==9)) = drawRnd([1,sum((scatt.type==7)|(scatt.type==9))],...
                                                   par.Scatt.power.pdf,par.Scatt.power.multiMov_par);  % draw for elliptical scatterer
  % now go through the special types, first the reflection+scattering
 for i = find((scatt.type==scatt.info.inside.reflScatt) | (scatt.type==scatt.info.inside.scattRefl))  % go through each type reflection + scattering
  if (scatt.type==scatt.info.inside.reflScatt)   % check for reflection->scattering type
    % assume first reflection->scattering, so check for possibilities on the walls
   [xr,wallInd,len] = findReflectors(walls,scatt.pDelay(i,:),tx_pos);  % calculate possible reflectors on the walls
   if isempty(xr)
    scatt.type = scatt.info.inside.scattRefl; % no possible reflection points on the walls found -> so set to scattering->reflection
   else
     % for this type of scatterer, we only have to add an additional delay
     % which is constant each time
    scatt.addDelay(i) = len(discrRnd(1,1,(1/length(len))*ones(1,length(len))))-norm(tx_pos-scatt.pDelay(i,:).');  % draw with equal prob. the wall and calculate the additional path length
   end;
  end; 
  if (scatt.type==scatt.info.inside.scattRefl)
    % assume now a scattering->reflection, this can be easily emulated
    % by calculating the delay not by scatt.pAmp but by using the mirrored
    % point instead and adding a delay which might be positive or negative.
   scatt.pDelay(i,:) = calcMirrorPoint(walls(discrRnd(1,1,(1/length(walls))*ones(1,length(walls)))),scatt.pAmp(i,:)); % take a uniformly drawn wall and calculate the mirrored point 
   scatt.addDelay(i) = norm(tx_pos(:)-scatt.pAmp(i,:).')-norm(tx_pos(:)-scatt.pDelay(i,:).');  % correcting factor, when scatt.pDelay(i,:) is used for delay calculation
  end;
 end; 
   % scatt.pDelay should be used for delay and AoA calculation while
   % scatt.pAmp is used for amplitude pattern calculation.
  % now we consider the type of moving scatterers
   % Actually, the ellipse/hyperbola is only redisplayed as an additional
   % delay which can be easily implemented by scatt.addDelay! This is done
   % as in the COST2100 channel model.
 ellScat = (scatt.type==scatt.info.inside.movScatt) | (scatt.type==scatt.info.outside.movScatt);  % ellipse scatterers
 if (sum(double(ellScat))>0)
  scatt.addDelay(ellScat) = drawRnd([sum(double(ellScat)),1],par.Scatt.ellip.addDelayPDF,par.Scatt.ellip.addDelayPar); % draw additional length for ellipse/hyperbola
 end;
return;


function refl = initWalls(walls,par)
% Function to initialise the amplitude response from the different walls
% as walls(i). par is the parameter struct for the multipath components.
 epsPol = 1e-4;  % minimum distance of a pole to the unit circle
 Nwalls = length(walls);  % number of walls
 for i = 1:Nwalls  % go through each wall separately
  arOrd = drawRnd([1,1],par.Wall.ar.order.pdf,par.Wall.ar.order.p)+1;  % draw the order of the AR process
  pol = (1-drawRnd([1,arOrd],par.Wall.ar.poles.absPdf,par.Wall.ar.poles.absPar)).*...
           exp(1i*drawRnd([1,arOrd],par.Wall.ar.poles.phasePDF,par.Wall.ar.poles.phasePar));  % draw complex poles
  while any(abs(pol)>=(1-epsPol))  % check for stability
   % disp('Pole outside unit circle or close to it, taking a new one!');
   idx = find(abs(pol)>=(1-epsPol));
   pol(idx) = (1-drawRnd([1,length(idx)],par.Scatt.ar.poles.absPdf,par.Scatt.ar.poles.absPar)).*...
               exp(1i*drawRnd([1,length(idx)],par.Scatt.ar.poles.phasePDF,par.Scatt.ar.poles.phasePar));  % draw new pole for stability reasons
  end;          
  if any(abs(pol)>=(1-epsPol))  % check for stability
   disp('Pole outside unit circle or close to it, please check!');
   keyboard;
  end; 
  refl.reflRes = par.reflRes;  % take over the resolution on the wall
  arCoef = poly(pol);  % calculate AR coefficients
  [cyy,lagCyy] = calcACF(30,arCoef,1);  % calculate covariance of the AR process assuming unit driving variance (we only need lag "0" here)
  [refl.wall(i).win_x,refl.wall(i).win_lag_x] = par.winPattern(ceil(norm(walls(i).v1)/par.reflRes)); % calculate window function in v1 direction (norm(v1) is the length of the wall in this direction)
  [refl.wall(i).win_y,refl.wall(i).win_lag_y] = par.winPattern(ceil(norm(walls(i).v2)/par.reflRes)); % calculate window function in v2 direction (norm(v2) is the length of the wall in this direction)  
   % now adapt the wall itself as extended wall to consider the full window
   % function.
  refl.extWall(i).s = walls(i).s+refl.wall(i).win_lag_x(1)*par.reflRes*unitVec(walls(i).v1)+...
                      refl.wall(i).win_lag_y(1)*par.reflRes*unitVec(walls(i).v2);
  refl.extWall(i).v1 = unitVec(walls(i).v1)*(refl.wall(i).win_lag_x(end)-refl.wall(i).win_lag_x(1))*par.reflRes;
  refl.extWall(i).v2 = unitVec(walls(i).v2)*(refl.wall(i).win_lag_y(end)-refl.wall(i).win_lag_y(1))*par.reflRes;
  refl.extWall(i).n = walls(i).n;   % take over the normal vector for plotting
   % now draw the amplitude process
  refl.wall(i).K = drawRnd([1,1],par.Wall.rice.pdf,par.Wall.rice.par);  % draw Rice factor for wall in linear domain
  refl.wall(i).amp_axes = par.reflRes*((refl.wall(i).win_lag_x(1)+refl.wall(i).win_lag_y(1)):(refl.wall(i).win_lag_x(end)+refl.wall(i).win_lag_y(end))); % x+y axis for the amplitude in [m] with (x in v1 direction and y in v2 direction)
  refl.wall(i).amp = (simAr(arCoef,1,length(refl.wall(i).amp_axes))+sqrt(refl.wall(i).K*cyy(1))*exp(1i*2*pi*rand(1)))/sqrt(cyy(1)*(1+refl.wall(i).K)); % calculate normalised amplitdue with AR process, constant and Rice factor
  if (any(abs(refl.wall(i).amp)>100))  % verification, should not happen
   disp('Big amplitude, please verify!');
   keyboard;
  end;
   % now draw the power for each wall, first and second order
  refl.av_pow(:,1) = drawRnd([1,Nwalls],par.Wall.power.pdf,par.Wall.power.single_par);  % draw for first order reflections
  refl.av_pow(:,2) = drawRnd([1,Nwalls],par.Wall.power.pdf,par.Wall.power.double_par);  % draw for double order reflections
   % The window function is provided as wall(i).win_x in v1 direction with lags wall(i).win_lag_x
   % and wall(i).win_y in v2 direction with lags wall(i).win_lag_y.
   % To calculate the window use win(x,y) = win_x(x)*win_y(y)!
%   figure; imagesc(refl.wall(i).win_lag_x*par.reflRes,refl.wall(i).win_lag_y*par.reflRes,repmat(refl.wall(i).win_x,length(refl.wall(i).win_y),1).*repmat(refl.wall(i).win_y(:),1,length(refl.wall(i).win_x))); colorbar; title('window function'); % plot window function
%   figure; plot(refl.wall(i).amp_axes,real(refl.wall(i).amp),'.-b',refl.wall(i).amp_axes,imag(refl.wall(i).amp),'.-r',refl.wall(i).amp_axes,abs(refl.wall(i).amp),'g'); legend('Real Part of norm Ampl','Imag Part of norm Ampl','Abs of Ampl'); grid on;
%   figure; plot(real(pol),imag(pol),'x',cos(0:0.1:2*pi+0.1),sin(0:0.1:2*pi+0.1),'-r'); grid on; axis equal; title('Pole distribution');
%   figure; plot(lagCyy,abs(cyy/cyy(1)),'.-'); grid on; title('Covariance function normalised');  
 end;
return;


function [leng,lagx,lagPr] = calcDelayPDF(walls,tx,grid,scatt,bandwidth)
% Function to calculate the delay distribution over the whole room.
 global simRefl;
 % ---- need to be done for the PDP normalisation ---- %
 c = 299792458;  % speed of light in [m/s]
 s = [walls.s];  % all points to describe the different walls
 [x,y,z] = meshgrid((min(s(1,:))+grid/2):grid:(max(s(1,:))-grid/2),...
                    (min(s(2,:))+grid/2):grid:(max(s(2,:))-grid/2),...
                    (min(s(3,:))+max(s(3,:)))/2);   % only one height
                    %(min(s(3,:))+grid/2):grid:(max(s(3,:))-grid/2));  % grid in z-direction
 rec = [x(:).';y(:).';z(:).'];  % all receiver points
%  plot3(rec(1,:),rec(2,:),rec(3,:),'.');  % plot for checking
 leng = [];  % delays in [s] normalised to the LoS path delay
 for i = 1:size(rec,2)  % go through each point
  losDist = norm(rec(:,i)-tx);   % LoS distance
   % --- go through the entry loss model for the position rec(:,i) --- %  
%   [E,distEntr] = calcPaths(tx(:),rec(:,i),walls,1e9,100e6);  % frequency and bandwidth not of interest here!
  distEntr = [];   % leave entry loss out for the PDF
   %  --- ENTRY LOSS PATHS LEFT OUT HERE ! --- %
%   plotGeometry(walls,rec(:,i),tx(:));
  leng = [leng,distEntr(:).'];  % add entry loss delays to the list (are already normalised to LoS length)
   % --- calculate reflector lengths for the position rec(:,i) --- %
  if (simRefl(1))   % check if first order shall be simulated
   [xr,wallInd,len,xrWall] = findReflectors(walls,rec(:,i),tx);  % first order reflections
   leng = [leng,(len-losDist)/c];   % ok, add first order reflections to the list (normalised to LoS)
  end;
  if (simRefl(2))  % check if second order shall be simulated
   [xrDouble,wallIndDouble,lenDouble,xrWallDouble] = findDoubleReflectors(walls,rec(:,i),tx); % double bounce reflections
   leng = [leng,(lenDouble-losDist)/c];   % add second order to the list (normalised to LoS)
  end;
   % --- now check for the scatterers --- %
  [inside,ang] = checkAngLimits(scatt.angLimits,scatt.pAmp,rec(:,i)); % check which ones are active at this receiver position
  dist = zeros(1,sum(double(inside)));  % init dist
  indScattAct = find(inside);        % get the index for the active scatterers
  NScattAct = length(indScattAct);   % number of scatterers visible
  if (NScattAct>0)    % if they exist
   dist(1:NScattAct) = norm_col(repmat(tx(:),1,NScattAct)-scatt.pDelay(indScattAct,:).')+...
                       norm_col(repmat(rec(:,i),1,NScattAct)-scatt.pDelay(indScattAct,:).')+scatt.addDelay(indScattAct).';  % path length calculation in [m]
  end;
   % now add them all to the list
  leng = [leng,dist-losDist];  % add all the scatterers to the delay list
  if (any(leng<0))  % check for weird error case
   disp('delay less than 0, setting to 0');
   leng(leng<0) = 0;
%    keyboard;
  end;
 end;
   % Probability in terms of samples 
   % leng is in [m], so leng/c*bandwidth is in [samples]
 maxBin = ceil(max(leng)/c*bandwidth+7);  % maximum Bin to be used
 bins = [0,(0.5:1:maxBin-0.5); (0.5:1:maxBin+0.5)];  
 lagPr = zeros(1,size(bins,2));
 for i = 1:length(lagPr)
  lagPr(i) = quad(@(x)kden(x,leng/c*bandwidth,[0,inf]),bins(1,i),bins(2,i));  % numerical integration from bins(1,i):bins(2,i)
 end;  
 lagx = 0:maxBin;
%  plot(lagx,lagPr,'.-'); title('Probability per Bin'); xlabel('lag l'); ylabel('Pr(l)'); 
% lagCDF = kden([0,0.5+1:(max(leng)/c*bandwidth)+7],leng/c*bandwidth,[0,inf],'cdf');  % estimate CDF using a kernel estimator (delay values in [samples])
% plot([0,0.5+1:(max(leng)/c*bandwidth)+7],lagCDF)
%  lagx = (0:0.1:max(leng))/c*bandwidth;    % lags in terms of samples
%  lagPDF = kden(lagx,leng/c*bandwidth,[0,inf]);  % estimate PDF using a kernel estimator (delay values in [samples])
%  plot(lagx,lagPDF);
   % PDF in terms of length [m]
%  lagx = (0:0.01:max(leng));  % in terms of length [m]
%  lagPDF = kden(lagx,leng,[0,inf]);  % estimate PDF using a kernel estimator (delay values in [samples])
  % leng are calculate according to the receiver positions and
  % corresponds to a realisation of a random variable. To calculate the pdf
  % of the leng, we need to use a kernel estimator.
%   plotpdf(leng); % for visualisation of the delay pdf
return;



%% ----- RunTime Functions ----- %%

function [E,dist,p,type,avPow] = calcScatt(tx,rx,scatt,freq,typeStr,walls)
% Function to calculate the wall reflections (first and second order)
% according to TX position tx and RX position rx. freq is the 
% carrier frequency for the simulation. 
% Output is: dist = path length of the path (without normalisation to LoS)
%            E = amplitude without the normalisation to the PDP
%            dist = path length [m] of the path (without normalisation to LoS)
%            p = interaction point to calculate AoA (not delay!), per row
%            avPow = average power value of the path
% Input is:  typeStr = type of path defined in typeStr 
%                       (see also createObj()) (only for debugging!)
%            walls = walls defining the room (only for debugging!)
 global interpMeth;
 c = 299792458;  % speed of light in [m/s]
  % first, we need to find the scatterer points where the amp>0, so the
  % ones which are within the amplitude window. Check via scatt.angLimits,
  % where scatt.angLimits(:,1) is the starting and scatt.angLimits(:,2) is
  % the end of the amplitude pattern. Take care, scatt.angLimits(:,1) might
  % be negative and scatt.angLimits(:,2) might be > 360 deg.
 [inside,ang] = checkAngLimits(scatt.angLimits,scatt.pAmp,rx); 
 E = zeros(sum(double(inside)),1);  % init vector of non-normalised amplitudes
 avPow = zeros(sum(double(inside)),1);  % init vector of average power for the paths
 dist = zeros(sum(double(inside)),1);  % init vector of path lengths
 p = zeros(sum(double(inside)),3);  % init vector of interaction points
%   % plot for testing
%  fig = plotRoom2d(walls);  
%  plotScatt2d(scatt.pAmp(~inside,:),scatt.angLimits(~inside,:),'b',fig);
%  plotScatt2d(scatt.pAmp(inside,:),scatt.angLimits(inside,:),'r',fig);
%  plot(rx(1),rx(2),'or');  % plot receiver position
  % now we got the ones active so save the type for verification
 type = scatt.type(inside).';   % get the type of active scatterers in the correct order
  % now we got the types of scatterers, let's go on and calculate delay and
  % amplitude as well as the point.
 indScattAct = find(inside);        % get the index for the active scatterers
 NScattAct = length(indScattAct);   % get the number of active ones
 if (NScattAct>0)    % if they exist
  dist(1:NScattAct) = norm_col(repmat(tx(:),1,NScattAct)-scatt.pDelay(indScattAct,:).')+...
                      norm_col(repmat(rx(:),1,NScattAct)-scatt.pDelay(indScattAct,:).')+scatt.addDelay(indScattAct).';  % path length calculation in [m]
  for i = 1:NScattAct % go through each one for amplitude interpolation
   if (ang(indScattAct(i))>=scatt.pat(indScattAct(i)).ang(1)) && (ang(indScattAct(i))<=scatt.pat(indScattAct(i)).ang(end)) % additional not necessary check for opening angle
    try
     avPow(i) = scatt.av_pow(indScattAct(i));  % get average path power of the scatterer
     E(i) = sqrt(avPow(i))*interp1(scatt.pat(indScattAct(i)).ang,scatt.pat(indScattAct(i)).amp,ang(indScattAct(i)),interpMeth)*exp(-1i*2*pi*freq*dist(i)/c);  % interpolate for non-normalised amplitude
    catch
     disp('interp error');
     keyboard;
    end;
   else
    disp('Error, angle is out of range, please check! Should not happen!');
    keyboard;
   end;
%    plot(scatt.pat(indScattAct(i)).ang,abs(scatt.pat(indScattAct(i)).amp)); 
  end;
  p(1:NScattAct,:) = scatt.pDelay(indScattAct,:);  % take over the interaction points (not necessarily the equivalent scatterer!)
 end;
  % at the end, we got for this function
  %   dist = path length [m] of the path (without normalisation to LoS)
  %   E = amplitude without the normalisation to the PDP
  %   p = interaction point to calculate AoA (not delay!)
  %   type = type of scatterer
return;


function [E,dist,p,type,avPow] = calcRefl(tx,rx,freq,refl,typeStr)
% Function to calculate the wall reflections (first and second order)
% according to TX position tx and RX position rx. freq is the 
% carrier frequency for the simulation.
% Output is: dist = path length of the path (without normalisation to LoS)
%            E = amplitude without the normalisation to the PDP
%            dist = path length [m] of the path (without normalisation to LoS)
%            p = interaction point to calculate AoA (not delay!), per row
%            type = type of path defined in typeStr (see also createObj())
%            avPow = average power of the path
 global simRefl;
 global interpMeth;
 c = 299792458;  % speed of light in [m/s]
 walls = refl.extWall;  % use the extended walls to include the full window function
 if (simRefl(1))  % check if first order shall be simulated
  [xr,wallInd,len,xrWall] = findReflectors(walls,rx,tx);  % first order reflections
 else
  wallInd = []; len = [];  % nothing to take into account
 end;
 if (simRefl(2))  % check if second order shall be simulated
  [xrDouble,wallIndDouble,lenDouble,xrWallDouble] = findDoubleReflectors(walls,rx,tx); % double bounce reflections
 else
  wallIndDouble = []; lenDouble = [];  % nothing to take into account
 end
 type = [typeStr.firstReflection*ones(length(len),1); typeStr.secondReflection*ones(length(lenDouble),1)]; % set type of paths for output
%   % for visual inspection
%  plotRoom(walls);  % plot room
%  plot3(rx(1),rx(2),rx(3),'xr'); ntx = unitVec(tx-rx);  % norm vector to transmitter
%  quiver3(rx(1),rx(2),rx(3),ntx(1),ntx(2),ntx(3),'r');  % plot LoS for direct
%  for i = 1:size(xr,2)   % plot single bounce reflection
%   plot3(xr(1,i),xr(2,i),xr(3,i),'ob'); quiver3(xr(1,i),xr(2,i),xr(3,i),ntx(1),ntx(2),ntx(3),'r'); % single bounce reflections
%   plot3([xr(1,i),rx(1)],[xr(2,i),rx(2)],[xr(3,i),rx(3)],'-b');
%  end;
%  for i = 1:size(xrDouble,1)   % plot single bounce reflection
%   plot3(xrDouble(i,:,1),xrDouble(i,:,2),xrDouble(i,:,3),'ob'); quiver3(xrDouble(i,1,1),xrDouble(i,1,2),xrDouble(i,1,3),ntx(1),ntx(2),ntx(3),'r'); % single bounce reflections
%   plot3([xrDouble(i,:,1),rx(1)],[xrDouble(i,:,2),rx(2)],[xrDouble(i,:,3),rx(3)],'-b');
%  end;
%   % end of plotting
 E = zeros(length(wallInd)+size(wallIndDouble,1),1); % init amplitude vector
 avPow = zeros(length(wallInd)+size(wallIndDouble,1),1); % init average power vector
 dist = zeros(length(wallInd)+size(wallIndDouble,1),1); % init distance vector
 p = zeros(length(wallInd)+size(wallIndDouble,1),3); % init interaction points, per row
 Nsingle = length(wallInd);   % number of single bounce reflections
  % calculate contributions from single bounce reflections
 for i = 1:Nsingle  % go through the single bounce reflections
  xi = xrWall(1,i)*(refl.wall(wallInd(i)).win_lag_x(end)-refl.wall(wallInd(i)).win_lag_x(1))+refl.wall(wallInd(i)).win_lag_x(1); % v1 direction
  yi = xrWall(2,i)*(refl.wall(wallInd(i)).win_lag_y(end)-refl.wall(wallInd(i)).win_lag_y(1))+refl.wall(wallInd(i)).win_lag_y(1); % v2 direction
  dist(i) = len(i);    % path length without normalisation to LoS
  avPow(i) = refl.av_pow(wallInd(i),1);  % get average power for a single reflection
  E(i) = sqrt(avPow(i))*interpWallAmp([xi,yi],refl.wall(wallInd(i)))*exp(-1i*2*pi*freq*dist(i)/c);   % interpolate amplitude and add path length with delay  
  p(i,:) = xr(:,i).';    % save interaction point
   % for testing of interpWallAmp (compares local to global interpolation) :
%   val = interpWallAmp([xi,yi],refl.wall(wallInd(i)));
%   winFunc = repmat(refl.wall(wallInd(i)).win_x(:),1,length(refl.wall(wallInd(i)).win_y)).*...
%             repmat(refl.wall(wallInd(i)).win_y(:).',length(refl.wall(wallInd(i)).win_x),1); % access winFunc as winFunc(x,y) where x is in direction of v1 and y in direction of v2!
%   [x,y] = meshgrid(refl.wall(wallInd(i)).win_lag_x,refl.wall(wallInd(i)).win_lag_y); % meshgrid for x and y coordinate matrices
%   winFunc = winFunc.'.*...
%                 refl.wall(wallInd(i)).amp(x+y-refl.wall(wallInd(i)).win_lag_x(1)-refl.wall(wallInd(i)).win_lag_y(1)+1); % window function times AR(p) process
%   val2 = interp2(x,y,winFunc,xi,yi,interpMeth)      
%   fprintf('difference is %g + %g *j\n',real(val-val2),imag(val-val2));
   % and some nice figures
%   figure; imagesc(refl.wall(wallInd(i)).win_lag_x*refl.reflRes,refl.wall(wallInd(i)).win_lag_y*refl.reflRes,abs(refl.wall(wallInd(i)).winFunc)); 
%   colorbar; title(sprintf('window function, wall %i',wallInd(i))); xlabel('direction v1'); ylabel('direction v2'); % plot window function
 end;
  % calculate contributions from double bounce reflections
 Ndouble = size(wallIndDouble,1);  % number of double bounce reflections
 for i = 1:Ndouble  % go through the double bounce reflections
  wallNr = wallIndDouble(i,2);  % wall to be considered for last reflection
  xi = xrWallDouble(i,2,1)*(refl.wall(wallNr).win_lag_x(end)-refl.wall(wallNr).win_lag_x(1))+refl.wall(wallNr).win_lag_x(1); % v1 direction
  yi = xrWallDouble(i,2,2)*(refl.wall(wallNr).win_lag_y(end)-refl.wall(wallNr).win_lag_y(1))+refl.wall(wallNr).win_lag_y(1); % v2 direction
  dist(i+Nsingle) = lenDouble(i);    % path length without normalisation to LoS
  avPow(i+Nsingle) = refl.av_pow(wallNr,2);  % get average power for a double reflection
  E(i+Nsingle) = sqrt(avPow(i+Nsingle))*interpWallAmp([xi,yi],refl.wall(wallNr))*exp(-1i*2*pi*freq*dist(i+Nsingle)/c);   % interpolate amplitude and add path length with delay
  p(i+Nsingle,:) = squeeze(xrDouble(i,2,:));    % save interaction point
%   figure; imagesc(refl.wall(wallNr).win_lag_x*refl.reflRes,refl.wall(wallNr).win_lag_y*refl.reflRes,abs(refl.wall(wallNr).winFunc)); 
%   colorbar; title(sprintf('window function, wall %i',wallNr)); xlabel('direction v1'); ylabel('direction v2'); % plot window function
 end;
  % at the end, we got for this function
  %   dist = path length of the path (without normalisation to LoS)
  %   E = amplitude without the normalisation to the PDP
  %   p = interaction point to calculate AoA (not delay!)
return;


function [amp] = normToPDP(E,dist,multi,pdp,avPow)
% Function to calculate the normalisation value for each path
% in order to simulate the correct PDP. 
% Input: E = normalised amplitude without PDP
%        dist = path length normalised to LoS path in [m]
%        multi = struct of multi.lengPath, multi.lagx, multi.lagPdf,
%                 see calcDelayPDF() for explanation
%        pdp = struct defining the PDP from calcEntryLoss().
%        avPow = Average powers for all the paths with contribution in E.
 global interpMeth;
 c = 299792458;  % speed of light in [m/s]
 if (strcmpi(pdp.type,'exp'))  % check for exponential decaying PDP 
  ndist = dist(:)/c*pdp.bandwidth;     % in terms of samples
  amp = E(:).*sqrt(pdp.gamma*exp(-pdp.beta*ndist)./(multi.lagPr(round(ndist)+1).'*sum(avPow)));
  if (any(multi.lagPr(round(ndist)+1)==0))
   disp('Something wrong, probability is 0!');
   keyboard;
  end;
 else
  error('calcMulti: PDP type is unknown!');
 end;
%  keyboard;
%  plotpdf(multi.lengPath); hold on; grid on;
%  plot(multi.lagx/pdp.bandwidth*c,multi.lagPdf/c*pdp.bandwidth,'k','LineWidth',2)
%  x = 0:0.1:100; plot(x,kden(x,multi.lengPath,[0,inf]),'c','LineWidth',2); title('Delay PDF'); xlabel('Path length [m]'); grid on;
return;



%% ----- Supplementary Functions ----- %%

function E = calcField(tx,rx,walls,freq,bandwidth,normLoS)
% Function to calculate the entry loss field for all walls given.
 if (~exist('normLoS','var'))
  normLoS = true;   % set to default
 end;
 E = zeros(1,size(rx,2));  % init output
 for r = 1:length(walls)
  E = E+calcFieldWall(tx,rx,walls(r),freq,bandwidth,normLoS); % calculate field contribution from wall r
 end;
return;


function [E,amp,dist,xp] = calcFieldWall(tx,rx,wall,freq,B,normLoS)
% Function to calculate the field contribution from the wall given with
% tx and rx positions. tx is assumed to be a single position while rx
% might have multiple positions, where each coloumn is one pos as [x,y,z].'.
% freq is the carrier frequency. B is the bandwidth.
% wall is one wall given as a struct as generated by
% createWalls() and addWindows(). If normLoS is true (default)
% then the calculated field component is normalised with
% respect to the LoS field component.
 global angPat;  % subroutine to calculate the angular pattern
 c = 299792458;  % speed of light in [m/s]
 meth = 6;       % method for window_integ (6=ALBANI 5=RUBI, 1=2d-Kirchhoff)
 tol = 1e-8;     % tolerance for integration
 sumUp = false;      % summing parameter for window_integ
 if (unitVec(wall.n).'*unitVec(tx)<=0)  % check for kappa(psi)
  E = zeros(1,size(rx,2));
  amp = []; dist = []; xp = [];   % null all the output
  return;
 end;
 Nrx = size(rx,2);   % number of receiver positions
 if isfield(wall,'win')
  Nwin = length(wall.win);    % number of windows 
 else
  Nwin = 0;
 end;
 % Routine does not return the distance Tx-Element-Rx so far !!! %
 Twall = wall.Twall;  % get transmission factor for wall
 dTxRx = norm_col(repmat(tx,1,Nrx)-rx); % TX-RX distance for each rx point
%  plotGeometry(wall,rx,tx); hold on;  % for testing
 tx_loc = doTransform(tx,wall.k,false);   % transform to local coordinate system
 rx_loc = doTransform(rx,wall.k,false);
 pl_wall = getPlane(wall);  % get the plane of the wall
 Dwall = window_integ(pl_wall,tx_loc,rx_loc,c/freq,sumUp,meth,tol).';  % perform integration for plane, Dwall is [1 x Nrx]
 [losWall,xpWall] = checkForLos(pl_wall,rx_loc,tx_loc);  % calculate point on the plane, xpWall is [3 x Nrx]
 distWall = norm_col(repmat(tx_loc(:),1,Nrx)-xpWall)+norm_col(rx_loc-xpWall)-dTxRx;  % calculate tx->point->rx distance normalised to LoS!
 if (isfield(wall,'win') && ~isempty(wall.win))  % check for windows
   % windows use the same affine transformation as the wall itself
  Twin = [wall.win(:).Twin];  % get transmission factor for windows, Twin is [1 x Nwin]
  pl_win = getPlane(wall.win);   % get planes of the window
  [losWin,xpWin] = checkForLos(pl_win,rx_loc,tx_loc); % cal. points on the plane for each window as [nWin x 3 x Nrx], take care if Nwin==1, then [3 x Nrx]  
  % plotGeometry(wall,rx(:,183),tx); % for graphically checking of special points
  distWin = zeros(Nwin,Nrx);  % init distance calc. 
  if (Nwin>1)  % case of several windows
   for i = 1:Nwin  % for each window
    distWin(i,:) = norm_col(repmat(tx_loc(:),1,Nrx)-reshape(xpWin(i,:,:),size(xpWin,2),Nrx))...
                    +norm_col(rx_loc-reshape(xpWin(i,:,:),size(xpWin,2),Nrx))-dTxRx; % cal. distance tx->point->rx normalised to LoS!
   end;
  else
    % single window to be considered
   distWin = norm_col(repmat(tx_loc(:),1,Nrx)-xpWin)+norm_col(rx_loc-xpWin)-dTxRx;
  end;
   % now calculate field contribution
  Dwin = window_integ(pl_win,tx_loc,rx_loc,c/freq,sumUp,meth,tol);  % calculate window contribution, Dwin is [Nrx x 1 x Nwin]
  if (Nrx>1)
   Dwin = squeeze(Dwin).'; % [Nwin x Nrx], squeeze out the middle singleton
  else
   Dwin = squeeze(Dwin); % [Nwin x Nrx], squeeze out the middle singleton
  end;
   % Dwin is now a vector of size [Nwin x Nrx]
 else
  Dwin = zeros(1,Nrx);  % no windows do exist, so take to zero
  Twin = 0;
  distWin = 0;
 end;
  % calculate the field strength, not normalised to LoS
 E = angPat(wall.n,tx)*Twall*(Dwall-sum(Dwin,1)).*sinc(distWall/c,B)... % wall contribution itself, distWall is already normalised to LoS
      +angPat(wall.n,tx)*Twin*(Dwin.*sinc(distWin/c,B));                % window contributions, distWin is already normalised to LoS 
 if (~exist('normLoS','var') || normLoS)
  E = E./(exp(-1i*2*pi*dTxRx*freq/c)./dTxRx);  % normalise to LoS
 end;
 if (nargout>1)  % check for additional output
   % has to be done here
  amp = angPat(wall.n,tx)*[Twall*(Dwall-sum(Dwin,1));...  % wall path
                           repmat(Twin(:),1,Nrx).*Dwin];  % window contributions
  if (~exist('normLoS','var') || normLoS)  % normalisation to LoS path amplitude
   amp = amp./repmat(exp(-1i*2*pi*dTxRx*freq/c)./dTxRx,size(amp,1),1);  % normalise each component to the LoS path amplitude
  end;
   % dist is normalised to the LoS path
  dist = [distWall;...        % wall path
          distWin];           % windows
  if (Nwin==0)
   amp = amp(1:end-1,:);  
   dist = dist(1:end-1,:);   % take off the last value, will be 0 as no windows do exist
  end;
  if (any(dist/dTxRx<1e-13))  % to prevent rounding errors, set very small dist values to zero
   dist(dist/dTxRx<1e-13) = 0;
  end;
   % interaction point for each contribution as [NrPath x 3 x Nrx]
  xp(1,:,:) = doTransform(xpWall,wall.k,true);  % inverse transformation from local wall to room coordinates
  if (Nwin>1)
    % several windows
   for i = 1:Nwin
    xp(i+1,:,:) = doTransform(reshape(xpWin(i,:,:),size(xpWin,2),Nrx),wall.k,true);  % inverse transformation from local wall to room coordinates
   end;
  elseif (Nwin==1)
    % single window
   xp(2,:,:) = doTransform(xpWin,wall.k,true);  % inverse transformation from local wall to room coordinates
  end;
 end;
return;


function pat = calcPatCosSquare(n,tx)
% Function to calculate the angular pattern as a function of the 
% wall normal given by n and the TX position tx.
 cos_phi = unitVec(n).'*unitVec(tx);  
 if (cos_phi>0)  
  pat = (cos_phi)^2;  % cos(phi)*cos(phi)
 else
  pat = 0;   % TX is not visible for the wall
 end;
return;


function h = sinc(del,B)
% Function to calculate the value for n=0 if the
% sinc function with the two-sided bandwidth B.
% B is in circular frequency. del should be given in [s].
 h = zeros(size(del,1),size(del,2));  % init
 h(del~=0) = sin(pi*B*(0-del(del~=0)))./(pi*B*(0-del(del~=0))); % sinc function
 h(del==0) = 1;
return;


function siz = getPlane(wall)
% Function to extract the size of the plane as [ymin, ymax, zmin, zmax].'
% as used for window_integ(). Output matrix is oriented in coloumns per
% plane.
 siz = zeros(4,length(wall));
 for i = 1:length(wall)
  plane = doTransform(repmat(wall(i).s,1,3)+[zeros(3,1),wall(i).v1,wall(i).v2],wall(i).k,false);
  siz(:,i) = [min(plane(2,:)),max(plane(2,:)),min(plane(3,:)),max(plane(3,:))].'; 
 end;
return;


function [los,xp] = checkForLos(window,r,t)
% Function to check for a direct path through the window
% Input:
%  window is defined as [ymin, ymax, zmin, zmax].' at x=0 (each coloumn a window)
%  r is the receiver as [x,y,z].' (one point per coloumn)
%  t is the transmitter as [x,y,z]
% Output: 
%  los = (true,false) if point is on the plane defined by window
%  xp = point at x=0 between t and r if los=true
%       otherwise xp will be on the edge closest to original xp
%        xp(:,i) is associated to window(:,i)!
 los = zeros(size(window,2),size(r,2));
 xp = zeros(size(window,2),3,size(r,2)); % init matrix of points
 for i = 1:size(window,2)  % go through each coloumn
  lam = -r(1,:)./(t(1)-r(1,:));  % calculate lambda for line from receiver to transmitter at x=0
  if (t(1)-r(1,:)==0)
   keyboard;
  end;
  v = r+repmat(lam,3,1).*(repmat(t(:),1,size(r,2))-r);  % position at x=0 at the sceen
  if any(abs(v(1,:))>1e3*eps)
   disp('wrong calculation, please check!');
   keyboard;
  end;
   % create output values
%   los = zeros(1,size(v,2));  % init los vector
  los(i,(v(2,:)>=window(1,i)) & (v(2,:)<=window(2,i)) & (v(3,:)>=window(3,i)) & (v(3,:)<=window(4,i))) = 1;  % LoS condition fulfilled
  xp(i,:,:) = v;  % initialise output
  xp(i,3,(los(i,:)==0).' & (squeeze(xp(i,3,:))>window(4,i))) = window(4,i); % fix nLoS points to border of plane
  xp(i,3,(los(i,:)==0).' & (squeeze(xp(i,3,:))<window(3,i))) = window(3,i);
  xp(i,2,(los(i,:)==0).' & (squeeze(xp(i,2,:))>window(2,i))) = window(2,i);
  xp(i,2,(los(i,:)==0).' & (squeeze(xp(i,2,:))<window(1,i))) = window(1,i);
%   % for testing
%  checkGeometry(t,r.',window)
%  plot3(v(1,:),v(2,:),v(3,:),'or');  % original point
%  plot3(xp(1,:),xp(2,:),xp(3,:),'xc'); % new point either inside or on the edge
%  for i = 1:size(r,2)
%   plot3([t(1),xp(1,i),r(1,i)],[t(2),xp(2,i),r(2,i)],[t(3),xp(3,i),r(3,i)],'-c');
%  end;
%  keyboard;
 end;
 if (size(window,2)==1)
  xp = reshape(xp,size(xp,2),size(xp,3));  % squeeze out first dimension if not needed
 end;
return;


function val = interpWallAmp(p,wall)
% Function to interpolate the window function including the amplitude
% of the given wall wall. p is provided as [x,y] using the axes 
% wall.win_lag_x and wall.win_lag_y.
 global interpMeth;
 rang = [4,4];   % number of points before and after  
 xInd = find(round(p(1))==wall.win_lag_x); % get index for x
 yInd = find(round(p(2))==wall.win_lag_y); % get index for y
 if (isempty(xInd) || isempty(yInd))
  fprintf(' Point %g, %g could not be found in index!\n',p(1),p(2));
  error('Something wrong! Should not happen!');
 end;
 xIndRang = max(1,xInd-rang(1)):min(length(wall.win_lag_x),xInd+rang(1)); % range including min and max
 yIndRang = max(1,yInd-rang(2)):min(length(wall.win_lag_y),yInd+rang(2));
 win = repmat(wall.win_x(xIndRang),length(yIndRang),1).*repmat(wall.win_y(yIndRang).',1,length(xIndRang)); % build up 2d around p x-> direction are the coloumns, y->direction are the rows
 [x,y] = meshgrid(wall.win_lag_x(xIndRang),wall.win_lag_y(yIndRang)); % meshgrid for x and y coordinate matrices for matrix win
 winFuncAR = win.*wall.amp(x+y-wall.win_lag_x(1)-wall.win_lag_y(1)+1); % window function times AR(p) process
  % axes of winFunc are refl.win_lag_x in terms of coloumns
  %                     refl.win_lag_y in terms of rows, lags are in [par.reflRes], so not in [m]!
 val = interp2(x,y,winFuncAR,p(1),p(2),interpMeth);  % perform the interpolation
return;


function [x,A] = calcEllipsePoint(f1,f2,p,t,alpha,test)
% Function to generate sample points of the 2d ellipse.
% Input: f1, f2 = foci in [x,y] 
%        p = semimajor axis-norm(center-f1) 
%        t = angles to be computed in degree.
% Ouput x = matrix with length(t) x 2 coordinates as [x,y]
 f1 = f1(:);  f2 = f2(:);  % just to make sure
 [x0,y0] = deal(0.5*(f2(1)-f1(1))+f1(1),0.5*(f2(2)-f1(2))+f1(2)); % calculate center point
 alp = atan2(f2(2)-f1(2),f2(1)-f1(1));   % calculate rotation angle
 a = p+norm([x0;y0]-f1(1:2));
 b = sqrt(a^2-norm([x0,y0]-f1(1:2).')^2);    % semiminor axis calc.
 x = [cos(alp),-sin(alp); sin(alp),cos(alp)]*[a*cos(deg2rad(t(:)).'-alp); b*sin(deg2rad(t(:)).'-alp)];
 x = x+repmat([x0;y0],1,length(t));
 x = x.';    % one point per row
 xorig = x;   % save for test plot
 if (any(imag(x)))
  disp('Somehow imaginary !! Please verify!');
  keyboard;
 end; 
 if exist('alpha','var') && (length(f1)==3) && (length(f2)==3)
   % Ok, everything given, so extend to a plane in the 3d-space
   % --  calculate 3d space coordinate system
  xax3 = (f2-f1);  % x-axis in 3d space
  yax3 = cross(unitVec(xax3),[0;0;1]);  % orthogonal vector between [0;0;1] and (f2-f1) projected to x/y plane
  n = unitVec(xax3);  % unit vector of xax3
   % turning matrix from http://de.wikipedia.org/wiki/Drehmatrix#Drehmatrizen_des_Raumes_R.C2.B3
  R = [n(1)^2*(1-cosd(alpha))+cosd(alpha), n(1)*n(2)*(1-cosd(alpha))-n(3)*sind(alpha), n(1)*n(3)*(1-cosd(alpha))+n(2)*sind(alpha);...
       n(2)*n(1)*(1-cosd(alpha))+n(3)*sind(alpha), n(2)^2*(1-cosd(alpha))+cosd(alpha), n(2)*n(3)*(1-cosd(alpha))-n(1)*sind(alpha);...
       n(3)*n(1)*(1-cosd(alpha))-n(2)*sind(alpha), n(3)*n(2)*(1-cosd(alpha))+n(1)*sind(alpha), n(3)^2*(1-cosd(alpha))+cosd(alpha)];
  yax3 = R*yax3;  % now calculate real y-axis for 3d space
  zax3 = cross(unitVec(xax3),yax3);  % z-axis for the coordinate system
   % --  calculate 2d space coordinate system
  xax2 = [f2(1:2);0]-[f1(1:2);0];  % x-axis along the major axis placed on x/y plane
  yax2 = cross(unitVec(xax2),[0;0;1]);  % orthogonal vector between [0;0;1] and (f2-f1) projected to x/y plane
  zax2 = cross(unitVec(xax2),yax2);  % z-axis for the coordinate system
   % now calculate projection matrix from 2d->3d space
  A = [xax3,yax3,zax3]/[xax2,yax2,zax2];  % projection matrix from 2d->3d
  x = A*[(xorig-repmat([x0,y0],size(xorig,1),1)).';zeros(1,size(xorig,1))]+repmat((f1+f2)/2,1,size(xorig,1));
  x = x.';  % now, one point per line
  if exist('test','var') && test
   figure; plot3([f1(1),f2(1)],[f1(2),f2(2)],[f1(3),f2(3)],'xr'); grid on; hold on;
   plot3([f1(1),f2(1)],[f1(2),f2(2)],[0,0],'xg'); 
   quiver3(x0,y0,0,xax2(1),xax2(2),xax2(3),'g');  % 2d, x-axis
   quiver3(x0,y0,0,yax2(1),yax2(2),yax2(3),'g');  % 2d, x-axis
   quiver3(x0,y0,0,zax2(1),zax2(2),zax2(3),'g');  % 2d, x-axis
   x0 = (f1+f2)/2;  % in 3d space
   quiver3(x0(1),x0(2),x0(3),xax3(1),xax3(2),xax2(3),'r');  % 2d, x-axis
   quiver3(x0(1),x0(2),x0(3),yax3(1),yax3(2),yax3(3),'r');  % 2d, x-axis
   quiver3(x0(1),x0(2),x0(3),zax3(1),zax3(2),zax3(3),'r');  % 2d, x-axis
   plot3(xorig(:,1),xorig(:,2),zeros(size(x,1),1),'.-g');
   plot3(x(:,1),x(:,2),x(:,3),'.-r');
   xlabel('x'); ylabel('y'); zlabel('z');
  end;
 end;
return;


function [inside,ang] = checkAngLimits(limits,pScatt,pRx)
% Function to check if the angle between pScatt(i,:) and pRx and checking, if
% the angle is within the limits [limits(:,1),limits(:,2]. For pScatt, many
% points may be given, where each point is a row while for pRx only one
% point is allowed. The check will be done based on 2d, so only for pScatt(i,1:2) and
% pRx(1:2). inside will be a boolean with true if the angle ang is within
% the provided limits. ang output will be provided in [deg]. 
 ang = atan2(pRx(2)-pScatt(:,2),pRx(1)-pScatt(:,1))*180/pi;  % angle in degree
 ang(limits(:,2)>360) = ang(limits(:,2)>360)+360;   % angular wrapping if upper bound is above 360 deg
 inside = (ang>=limits(:,1)) & (ang<=limits(:,2));  % boolean for the angle to be inside the limits
return;


function [w,lag] = cosWin(N,ramp)
% Function to provide a window for smooth on/off of the amplitude
% pattern. The window here uses a cos-function to provide the ramping.
% N is the number of points for condition on and the 
% before and after will be of length ramp*N. 
%  lag = (-ceil(ramp*N):(N+ceil(ramp*N)-1));
%  w = [0.5*cos((pi/(ramp*N))*(-ceil(ramp*N):-1))+0.5,ones(1,N),0.5*cos((pi/(ramp*N))*(1:ceil(ramp*N)))+0.5];
  % half rising within window
 w = [0.5*cos((pi/(ramp*N))*(-ceil(ramp*N):-1))+0.5,ones(1,N-ceil(ramp*N)),0.5*cos((pi/(ramp*N))*(1:ceil(ramp*N)))+0.5];
 if (length(w)<3)
  w = [0.25,0.5,0.25];  % take a kind of minimum window for interpolation by interp1
 end;
 lag = (-ceil(ramp*N*0.5):(length(w)-ceil(ramp*N*0.5))-1);
 w = w/sqrt(mean(w.^2));   % normalise power per sample
 if (length(lag)~=length(w))
  disp('length is unequal');
  keyboard;
 end;
 %  plot(lag,w,'.-'); grid on;
return;


function [within] = checkWinthinRoom(p,walls)
% Function verifies if the point p is within the room given by walls.
% The routine uses the Jordan curve theorem to determine if p is inside.
% p must be one point per row, so [Npoints x 3].
% walls is the typical room structure created by createRoom().
% Output will be a [1 x Npoints] with true for point is inside.
 Npoints = size(p,1);   % number of points to be checked
 v = [1,0,0].';  % direction from point
 count = zeros(1,Npoints);  % counter for each point
 within = false(1,Npoints);  % set everything to false first
 for i = 1:length(walls)   % go through each wall
  if ~isinf(cond([v,-walls(i).v1,-walls(i).v2]))  % check for rank deficiancy
   lam = [v,-walls(i).v1,-walls(i).v2]\(repmat(walls(i).s,1,Npoints)-p.');  % calc line intersection point [3 x Npoints]
     % check if ray with s=p+lam(1)*v is hitting the plane
   count((lam(1,:)>0) & ((lam(2,:)>=0) & (lam(2,:)<=1)) & ((lam(3,:)>=0) & (lam(3,:)<=1))) = ...
     count((lam(1,:)>0) & ((lam(2,:)>=0) & (lam(2,:)<=1)) & ((lam(3,:)>=0) & (lam(3,:)<=1)))+1;  % line hitted, so add 1
  end;
 end;
 within(mod(count,2)==1) = true; % all counts which are odd are inside the room
  % check graphically 
%  plotRoom(walls); hold on; plot3(p(within,1),p(within,2),p(within,3),'xr'); title('points within');  % points within
%  plotRoom(walls); hold on; plot3(p(~within,1),p(~within,2),p(~within,3),'xb'); title('points outside');  % points outside
return;


function [xr,wallInd,len,xrWall] = findReflectors(walls,rx,tx)
% Function to calculate the reflection point at walls "walls" for a transmission
% from tx to rx. xr(:,i) is the reflection point on the wall wallInd(i)
% referencing to walls(i). len(i) refers to the path length
% len(i)=norm(tx-xr(:,i))+norm(xr(:,i)-rx). If isempty(xr) then no
% reflections on walls(:) could be found.
% xrWall is the same as xr, just in terms of wall coordinates with
% xrWall(1,i) in direction walls(wallInd(i)).v1 with [0,1] and 
% xrWall(2,i) in direction walls(wallInd(i)).v2 with [0,1].
 xr = zeros(3,0);
 wallInd = zeros(0,1);
 len = zeros(0,1);
 xrWall = zeros(2,0);
 for i = 1:length(walls)
  [xmir] = calcMirrorPoint(walls(i),rx(:));  % calculate the mirror point of rx by plane walls(i)
   % now get the reflection point on the wall
     % inversion may throw a warning if a reflection is not possible
  r = [walls(i).v1(:),walls(i).v2(:),-(tx(:)-xmir)]\(xmir-walls(i).s(:)); % r = [m,n,lambda] with g: x = xr + lambda * (xt-xr)
  if ((r(1)<=1) && (r(1)>=0)) && ((r(2)<=1) && (r(2)>=0)) && ((r(3)>=0) && (r(3)<=1))  % check if point is still on the wall
   xr(:,end+1) = xmir+r(3)*(tx(:)-xmir);  % calculate reflection point
   xrWall(:,end+1) = [r(1);r(2)];  % coordinates on the wall in [0,1] for vector v1 and v2
   wallInd(end+1) = i;    % add wall to list
   len(end+1) = norm(xr(:,end)-tx(:))+norm(xr(:,end)-rx(:));  % path length including reflection
  end;
 end;
return;


function [xr,wallInd,len,xrWall] = findDoubleReflectors(walls,rx,tx)
% Function to calculate the reflection points on the walls 
% for double reflection from transmitter tx to receiver point rx.
% wallInd are the indices for the walls as wallInd(i,:) for the i-th
% double bounce reflection. Please note wallInd(i,1) is the first bounce
% and wallInd(i,2) is the second bounce seen from the transmitter.
% len(i) is the path length TX->refl(1)->refl(2)->RX. xr are the 
% points with xr(i,1,:) and xr(i,2,:) for the first and the second bounce.
% xrWall is the same as xr, just in terms of wall coordinates with
% xrWall(i,1,1) in direction walls(wallInd(i,1)).v1 with [0,1] and 
% xrWall(i,1,2) in direction walls(wallInd(i,1)).v2 with [0,1]. Similar
% for xrWall(i,2,:) for walls(wallInd(i,2)).
 wallInd = zeros(0,2);
 xr = zeros(0,2,3);
 xrWall = zeros(0,2,2);
 len = zeros(0);
 for i = 1:length(walls) % first wall reflection TX->walls(i)
  for l = 1:length(walls) % second wall reflection TX->walls(i)->walls(l)->RX
   if (l~=i)  % do not intend to check a double reflection on the same wall!!
    % we consider always walls(i) as the first reflection (so TX is
    % mirrored) while for walls(l) is the second reflection, so RX is
    % mirrored.
    xm_t = calcMirrorPoint(walls(i),tx(:));  % calculate mirror of TX to wall walls(i)
    xm_r = calcMirrorPoint(walls(l),rx(:));  % calculate mirror of RX to wall walls(l)
     % now get the reflection point on the wall walls(i) so on the first
     % wall seen from the TX
    r1 = [walls(i).v1(:),walls(i).v2(:),-(xm_r-xm_t)]\(xm_t-walls(i).s(:)); % r = [m,n,lambda] with g: x = xr + lambda * (xt-xr)
     % now get the reflection point on the wall pairs(i,2) so on the second
     % wall seen from the TX
    r2 = [walls(l).v1(:),walls(l).v2(:),-(xm_r-xm_t)]\(xm_t-walls(l).s(:)); % r = [m,n,lambda] with g: x = xr + lambda * (xt-xr) 
%      % for visual inspection
%     plotRoom(walls);  % plot room
%     plot3(rx(1),rx(2),rx(3),'xr'); plot3(tx(1),tx(2),tx(3),'xg'); ntx = unitVec(tx-rx);  % norm vector to transmitter
%     quiver3(rx(1),rx(2),rx(3),ntx(1),ntx(2),ntx(3),'r');  % plot LoS for direct
%     xr(:,:) = [(walls(i).s(:)+r1(1)*walls(i).v1(:)+r1(2)*walls(i).v2(:)).';...
%                 (walls(l).s(:)+r2(1)*walls(l).v1(:)+r2(2)*walls(l).v2(:)).'];  % save reflection locations
%     plot3([xr(:,1);rx(1)],[xr(:,2);rx(2)],[xr(:,3);rx(3)],'-ob');
%     plot3([xm_t(1),xm_r(1)],[xm_t(2),xm_r(2)],[xm_t(3),xm_r(3)],'-og');
     % check if point is still on the walls
    if ((r1(1)>=0) && (r1(1)<=1)) && ((r1(2)>=0) && (r1(2)<=1)) && ((r1(3)>=0) && (r1(3)<=1)) && ...
       ((r2(1)>=0) && (r2(1)<=1)) && ((r2(2)>=0) && (r2(2)<=1)) && ((r2(3)>=0) && (r2(3)<=1))
       % found a valid double reflection
     wallInd(end+1,:) = [i,l];  % get wall pairing
     xr(end+1,:,:) = [(walls(i).s(:)+r1(1)*walls(i).v1(:)+r1(2)*walls(i).v2(:)).';...
                      (walls(l).s(:)+r2(1)*walls(l).v1(:)+r2(2)*walls(l).v2(:)).'];  % save reflection locations, point per line
     xrWall(end+1,:,:) = [r1(1),r1(2); r2(1),r2(2)];  % points in wall coordinates [0,1], point per line
     len(end+1) = norm(xm_t-xm_r);  % double reflection length                 
    end;
   end; 
  end;
 end;
return;


function [xmir] = calcMirrorPoint(wall,p)
% Function to calculate the mirror point of p by plane wall defined
% as x = walls(i).s + m*walls(i).v1 + n*walls(i).v2.
   % calculate the mirrored point of p by the wall wall
   %  the wall is defined as x = wall.s + m*wall.v1 + n*wall.v2
 v = cross(wall.v1(:),wall.v2(:));  % calculate orthogonal vector to plain
 r = [wall.v1(:),wall.v2(:),-v]\(p(:)-wall.s(:));  % r = [m,n,lambda] with g: x = p + lambda * v
 xmir = p(:)+2*r(3)*v;   % calculate mirrored point on line g -> so rx is mirrored on plane as point xm
return;


function y = simAr(ar,vare,N)
% Function to simulate an AR process of N samples with filter coefficients
% ar and driving noise variance e.
% See also, "Autoregressive Modeling for Fading Channel Simulation",
% K. E. Baddour and N. C. Beaulieu, IEEE Trans. on Wire. Com., 
% vol. 4, No. 4, Jul. 2005, p. 1650-1662
% PLEASE NOTE: THIS DOES NOT REALLY WORK, SO WE USE nlags AND
% DISCARD THEM LATER ON.
  % initialisation of the AR(p) process according to paper, p. 1655
 minCyy = 1e-4;  % minimum covariance value to define nlag 
 [cyy,lags,idx] = calcACF(length(ar)+1,ar,vare,minCyy); % covariance function, own routine instead of poly2ac()
 nlags = idx+length(ar)+100;  % set up burn process according to covariance
 %cyy = poly2ac(ar,vare);  % calculate covariance function from AR(p)
%  v = cyy(1);
%  x_init(1) = cplxNormRnd(1,1,0,v);  % generate first sample
%  for i = 2:(length(ar)-1)
%   [ar_i,vare_i] = ac2poly(cyy(1:i));  % get limited AR(i-1)
%   v = vare_i;  % should be the same as in the next line
% %   v = (1-abs(ar_i(end))^2)*v;      % iteratively calculate variance
%   x_init(i) = cplxNormRnd(1,1,-ar_i(2:end)*x_init(:),v);  % draw new value
%  end;
%  y = filter(1,ar,cplxNormRnd(N+nlags,1,0,vare),x_init);  % calculate filering
 y = filter(1,ar,cplxNormRnd(N+nlags,1,0,vare));  % calculate filering, ignore now the paper at all!!!
 y = y(end-N+1:end);   % take only the last N values
 if (any(abs(y/sqrt(cyy(1)))>100))  % verification, should not happen
  disp('Big amplitude, please verify!');
  keyboard;
 end;
return;


function [cyy,lag,idx] = calcACF(maxLag,ar,vare,cyyTol)
% Function to calculate the autocovariance function for a given
% AR(p) process provided by the coefficients ar and the noise variance
% vare. The process is defined by 
%     X(t)*ar(1)+...+X(t-p)*ar(p+1)=E(t).
% The minimum length of cyy is restricted to p+1, so will never be smaller
% than that even if maxLag<p+1!
% If cyyTol is provided, the function tries to find the lag, where
% the covariance abs(cyy(l)/cyy(0)) is permanent below. Here 10 values are used.
% The value will be returned as idx.
%  -- Function can be checked by poly2ac(ar,vare) of Matlab! 
 minLags = 10; % minimum lags below cyyTol
 if exist('cyyTol','var') && (maxLag<minLags)
  maxLag = minLags;  % set minimum range for evaluation
 end;
 p = length(ar)-1;  % order of the process
 ar = ar(:).';   % produce a line vector
 A = zeros(2*(p+1),2*p+1);   % init matrix
 for i = 1:p+1   % go through each line
  A(i,:) = [real(ar(i:end)),zeros(1,i-1),imag(ar(i+1:end)),zeros(1,i-1)];
  A(p+1+i,:) = [imag(ar(i:end)),zeros(1,i-1),-real(ar(i+1:end)),zeros(1,i-1)];
  if (i>1)     % not in the first line
   A(i,:) = A(i,:)+[0,real(ar(i-1:-1:1)),zeros(1,p+1-i),-imag(ar(i-1:-1:1)),zeros(1,p+1-i)];
   A(p+1+i,:) = A(p+1+i,:)+[0,imag(ar(i-1:-1:1)),zeros(1,p+1-i),real(ar(i-1:-1:1)),zeros(1,p+1-i)];
  end;
 end;
 r = A\[vare;zeros(2*(p+1)-1,1)];  % calculate r = inv(A)*cye, where elements of r are divided into real and imaginary part
 cyy = [r(1);r(2:p+1)+1i*r(p+2:end)]; % construct cyy covariance vector
 lag = 0:max(maxLag,p);
  % up to here, cyy is calculated until lag=p, starting with lag=0
 if maxLag>p  % need to extend ?
  cyy(maxLag+1) = 0;  % add one for initialisation
  for l = (p+1):maxLag   % calculated already lags 0..p
   cyy(l+1) = -(1/ar(1))*sum(ar(2:end).*cyy(l:-1:(l-p+1)).');
  end;
 end;
 if (exist('cyyTol','var'))
   % Now we need to find the lag, where the covariance is very low for
   % minLags values.
  idx = minLags;
  cov = cyy(1:max(minLags,p));
  while (any(abs(cov/cyy(1))>cyyTol) && (idx<1e7))
   idx = idx+1;   % increase last lag
   cov = [cov(2:end); -(1/ar(1))*sum(ar(2:end).*cov(end:-1:(end-p+1)).')]; % calculate next value and save into vector
  end;
  if (idx==1e7)
   disp('Iteration too long!');
   keyboard;
  end;
 end;
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


function x = unitVec(x)
% Calculates the unit vector.
 x = x/norm(x);
return;

function d = norm_col(mat)
% Function calculates the norm over all coloumns in matrix mat.
% d is then [1 x n] if mat is [m x n] with d(i) = norm(mat(:,i))
 d = sqrt(sum(abs(mat).^2,1));  % norm over coloumn
return;


function deg = rad2deg(ang)
% Function to convert from radian to degree.
 deg = ang*180/pi;
return;

function ang = deg2rad(deg)
% Function to convert from degree into radian.
 ang = deg*pi/180;
return;


%% ----- Routines for random variables ----- %%

function x = drawRnd(sizeOut,dist,par)
% Function to draw from a pdf defined as string dist with
% parameters par. sizeOut = [m,n].
 if strcmpi(dist,'lognpdf')  
   % draw from log-normal distribution
  x = exp(randn(sizeOut).*par(2)+par(1));  % par(1)=mean, par(2)=std for log(Normal(mean,std))
 elseif strcmpi(dist,'Weibull')
   % draw from Weibull distribution
  x = par(1).*(-log(rand(sizeOut))).^(1./par(2)); 
 elseif strcmpi(dist,'Normal')
   % draw from Gaussian distribution
  x = randn(sizeOut).*par(2) + par(1); % draw normal rvs
 elseif strcmpi(dist,'truncNormal')
   % draw from truncated Gaussian distribution
  [x] = truncGauss(sizeOut(1)*sizeOut(2),par(1),par(2));  % can do only [0,inf] truncated normal
  x = reshape(x,sizeOut(1),sizeOut(2));
 elseif strcmpi(dist,'geometric')
   % draw from geometric distribution
  x = ceil(abs(log(rand(sizeOut))./log(1-par(1)))-1); 
  x(x < 0) = 0;
 elseif strcmpi(dist,'neg. binom.')
   % draw from a negative binomial distribution
   % from R. Saucier, "Computer Generation of Statistical Distributions",
   % Army Research Laboratory, ARL-TR-2168, March 2000
   % Works only for par(1) as integer !!!
%    x = sum(drawRnd([sizeOut(1)*sizeOut(2),par(1)],'geometric',par(2)),2); % use the geometric distribution 
%    x = reshape(x,sizeOut(1),sizeOut(2));  
  x = nbinrnd(par(1),par(2),sizeOut); 
 elseif strcmpi(dist,'dicrete')
   % draw from a discrete distribution defined by the 
   % probability mass function par with values 0:length(par)-1.
  x = discrRnd(sizeOut(1),sizeOut(2),par,0:length(par)-1);
 elseif strcmpi(dist,'fixed')
   % this is a wrapper for using a fixed value, so no random variable
   % will be created. Simply x = par for all values.
  x = par(1)*ones(sizeOut);    
 elseif strcmpi(dist,'truncExpPDF')
   % draw from truncated exponential PDF in ]0,1]
  r = rand(sizeOut)*(1-exp(-par(1)));
  x = -log(1-r)*1/par(1);
  while any(x==0)   % check if any values are exactly "0", because of ]0,1]
   r = rand(1,sum(double(x==0)))*(1-exp(-par(1)));
   x(x==0) = -log(1-r)*1/par(1);
  end;
 elseif strcmpi(dist,'truncWeibullPDF')
   % draw from truncated Weibull PDF in ]0,d]
  d = 1;   % upper parameter limit
  K = 1/(1-exp(-(d/par(1))^par(2)));  % multiplicative normalisation constant, using par(1)=lam, par(2)=k 
   % http://en.wikipedia.org/wiki/Weibull_distribution
  r = rand(sizeOut);
  x = par(1)*(-log(1-r/K)).^(1/par(2));  % draw random variables
 elseif strcmpi(dist,'invGauss')
   % draw from an inverse Gaussian distribution (also called Wald
   % distribution), see http://en.wikipedia.org/wiki/Inverse_Gaussian_distribution
   y = (randn(sizeOut)).^2;
   x = par(1)+par(1)^2*y/(2*par(2))-(par(1)/(2*par(2)))*sqrt(4*par(1)*par(2)*y+(par(1)*y).^2);
   z = rand(sizeOut);
   x(z>(par(1)./(par(1)+x))) = par(1)^2./x(z>(par(1)./(par(1)+x)));    
 elseif strcmpi(dist,'mixUniMisesPDF')
   % draw from mixture between uniform and von Mises distribution
%   i = discrRnd(sizeOut(1)*sizeOut(2),1,[par(1),1-par(1)]);  % draw discrete variables first
%   x(i==1) = uniRnd(sum(i==1),[-180,180]);
%   x(i==2) = circ_vmrnd(par(2)*pi/180,par(3),sum(i==2))*180/pi;  % von Mises distribution in degree
 else
  fprintf('Distribution asking for is ''%s''\n',dist);
  error('Not aware of this distribution!');
 end;
return;

function x = truncGauss(N,mu,sig)
% Implementation of the paper
 alpSig =  1.13671779106551;  % alpha*sigma in Appendix C, alpSig = fzero(@(alpSig)(exp(1)*sqrt(2/pi)-alpSig*exp(alpSig.^2/2)),1)
 mua = ((1-alpSig^2)/alpSig)*sig;  
 mub = 0;
 muc = sqrt(pi*sig^2/2);
 %c = sqrt(pi*sig^2/2)*(1+erf(mu/sqrt(2*sig^2))); % constant for the normal distribution
 alpha = (-mu+sqrt(mu^2+4*sig^2))/(2*sig^2); % constant for exp-distribution (4), a bit different to the paper!!
 x = zeros(1,N);  % initialise all to be negative
 lx = 0;
 while lx<N      % check if we still need values
  num = N-lx;    
  if (mu<=mua)
    % take distribution 4 (exponential), region 1
   z = -log(1-rand(1,num))/alpha;  % draw z~g(z)
   rho = exp(-(z-mu).^2/(2*sig^2)-alpha*(mu-z+alpha*sig^2/2));
  elseif ((mu>mua) && (mu<=mub))
    % take distribution 3, region 2
   z = abs(sig*randn(1,num))+mu;
   rho = double(z>=0);       % take only positive z as rho(z)
  elseif ((mu>mub) && (mu<=muc))
    % take distribution 2, region 3
   d = (rand(1,num)<mu/(mu+sqrt(pi/2)*sig));  % if d==1, then uniform, otherwise normal distribution
   u = rand(1,num)*mu;              % uniform distribution for 0<u<mu
   r = abs(sig*randn(1,num))+mu;    % Gaussisn for x>mu
   z = double(d).*u+(1-double(d)).*r;   % mixture of uniform and Gaussian
   rho = double(d).*exp(-(z-mu).^2/(2*sig^2))+(1-double(d)).*ones(1,num);
  elseif (mu>muc)
    % take distribution 1, region 4
   z = sig*randn(1,num)+mu;
   rho = double(z>=0);       % take only positive z as rho(z)
  else
   disp('This should not happen, please check!');
   keyboard;
  end;  
  z = z(rand(1,num)<rho);         % perform rejection
  if (lx+length(z)<=N)
   x(lx+1:lx+length(z)) = z;
  else
   x(lx+1:end) = z(1:N-lx);
  end;
  lx = lx+length(z);
 end;
return;



function x = uniRnd(m,int)
% Function to draw from a uniform distribution.
% x is [m x 1] in the interval [int(1)..int(2)].
% For each x(i) a different interval may be used as [int(i,1)..int(i,2)]
 if (~exist('int','var'))
  int = [0,1];  % standard
 end;
 if isscalar(int) || ((size(int,1)~=2) && (size(int,2)~=2)) || ...
     ((size(int,1)*size(int,2)>2) && (size(int,1)*size(int,2)/2~=m))
  error('Please provide a correct interal!');
 end;
 if (size(int,1)*size(int,2)>2)
  if (size(int,2)~=2)
   int = int.';   % need to transpose the matrix
  end;
 end;
 int = sort(int,2);  % sort via row
 x = (int(:,2)-int(:,1)).*rand(m,1)+int(:,1); % draw uniform from standard
return;

function x = discrRnd(m,n,pmf,values)
% function to create a discrete random variable with
% probability mass function pmf. If ~exist('values','var'), then
% calues will be 1..length(pmf), otherwise values will be used instead.
% Automatically it will be assumed that P(values(i))=pmf(i).
 u = rand(m,n); % draw uniform
 x = zeros(m,n);  % pre-allocate memory
 for i = 1:length(pmf) 
  if (i==1)
   x(u<sum(pmf(1:i))) = i;
  else
   x((u>=sum(pmf(1:i-1))) & (u<sum(pmf(1:i)))) = i;
  end;
 end;
 if exist('values','var') && (length(values)==length(pmf))
  x = values(x);
 end;
return;

function x = cplxNormRnd(m,n,mu,cov)
% Function to draw random variables from a complex
% cyclic normal distribution with mean mu and covariance matrix
% cov. x will be of size [m x n]. Real and imaginary part are considered
% to be independent with the same covariance.
 x = (randn(m,n)+1i*randn(m,n)).*(sqrtm(cov)/2) + mu; % draw normal rvs
return;



%% ----- Testing Routines ----- %%

function testRoutine(par)
% Function used for various tests
  % testing of AR(p) process simulation
%  roots = [0.5*exp(j*pi/4),0.96*exp(-j*pi/8),j*0.6,0.3];  % complex filter
%  vare = 1;   % driving noise variance
%  cyyTol = 1e-1;
%  ar = poly(roots);
%  y = simAr(ar,vare,1000);
%  z = filter(ar,1,y);  % inverse filtering
%  [cyy,lag,idx] = calcACF(400,ar,vare,cyyTol);
%  plot(lag,abs(cyy/cyy(1)),'.-b',lag,cyyTol*ones(1,length(lag)),'-r',idx*ones(1,2),[0,1],'g');
%  figure; plot(1:length(y),real(y),'.-b',1:length(y),imag(y),'.-r'); title('AR process');
%  figure; plot(1:length(z),real(z),'.-b',1:length(z),imag(z),'.-r'); title('Inverse AR process');
%  figure; plot(real(roots),imag(roots),'x',cos(0:0.1:2*pi+0.1),sin(0:0.1:2*pi+0.1),'-r'); grid on; axis equal; title('Pole distribution');
%  figure; plot(lag,abs(cyy),'.-'); title('Covariance of the process');
% ------------------------------
%   % test the window pattern for the amplitude of processes
 N = 10000;  % should be one within this number of points, so equal to "1" 
 [y,lag] = par.winPattern(N);   % return pattern
 plot(lag*1e-3,y,'.-'); grid on;
% ------------------------------
  % test for RVs and distributions
%  x = drawRnd([1,100000],'mixUniMisesPDF',[0.4,2,100]);  % generate mixed uniform von Mises distribution RVs
%  plotpdf(x,[],[],@(x)mixUniMisesPDF(x,0.4,2,100));
%  x = drawRnd([1,100000],'truncExpPDF',[10]);  % generate truncated EXP [0,1]
%  plotpdf(x,[],[],@(x)truncExpPDF(x,10,[0,1]));
%  x = drawRnd([1,10000],'truncNormal',[0.2,2]);  % generate truncated Normal RVs [0,\inf]
%  plotpdf(x,[],[],@(x)truncNormPDF(x,0.2,2));
%  y = 0:0.01:10;
%  plot(y,kden(y,x,[0,100]),'c','LineWidth',2);
%  x = drawRnd([1,10000],'truncWeibull',[1,5]);  % generate truncated Weibull RVs ]0,1]
%  plotpdf(x,[],[],@(x)truncWeibullPDF(x,1,5));
%  y = 0:0.001:1;
%  plot(y,kden(y,x,[0,1]),'c','LineWidth',2);
%  x = drawRnd([1,10000],'invGauss',[1,3]);  % generate inverse Gaussian RVs ]0,\inf]
%  plotpdf(x,[],[],@(x)invGaussPDF(x,1,3));
%  y = 0:0.001:5;
%  plot(y,kden(y,x,[0,inf]),'c','LineWidth',2);

 x = drawRnd([1,10000],'neg. binom.',[5,0.1]);  % generate neg. binomial RVs ]0,\inf]
 plotpdf(x,[0:1:max(x)],[],@(x)negBinomialPDF(x,5,0.1));
 y = 0:1:max(x);
 plot(y,kden(y,x,[0,inf]),'c','LineWidth',2);
 
%  % test a pdf for integration
%  x = 0:0.01:10000;
%  y = truncNormPDF(x,2,10);
%  sum(y)*0.01  % perform integration
% ------------------------------
%  % 2d checking of ellipse calculation 
%  p = 2; f1 = [-1,0,0]; f2 = [1,0,0]; t = 0:0.1:90;  % angle in degree
%  x = calcEllipsePoint(f1,f2,p,t);  % 2d point on x/y plane
%  figure; plot(x(:,1),x(:,2),'.-'); grid on; title('calcEllipsePoint() test');
% ------------------------------
%  % 3d checking of ellipse calculation  
%  alpha = 0;
%  p = 10.96; f1 = [29.5,0,5.2]; f2 = [14.3,0.5,13.6]; t = 0:0.1:360;
%  x = calcEllipsePoint(f1,f2,p,t,0);  % 3d points, with alpha=0
%  x2 = calcEllipsePoint(f1,f2,p,t,alpha,true);  % 3d points, with alpha=0
%  figure; plot3(x(:,1),x(:,2),zeros(size(x,1),1),'.-'); grid on; title('calcEllipsePoint() test'); hold on;
%  plot3(x(:,1),x(:,2),x(:,3),'og'); plot3([f1(1),f2(1),f1(1),f2(1)],[f1(2),f2(2),f1(2),f2(2)],[f1(3),f2(3),0,0],'xr');
%  plot3(x2(:,1),x2(:,2),x2(:,3),'.-r');
 
 keyboard;
return;


function simLoSMilan(gridSize)
% First simulation, LOS of Milan's paper compared to the model
%  c = 299792458;  % speed of light in [m/s]
 % define a room!!
  % define the size of the room (origin will be in the center of the room)
 par.room.dx = 10;  % x-dimension [m]
 par.room.dy = 10;  % y-dimension [m]
 par.room.dz = 3;   % z-dimension [m] (height)
 par.room.winWidth = 1;  % window width [m]
 par.room.winUpper = 0.3; % distance to upper edge for the windows [m]
 par.room.winLower = 1.2; % distance to lower edge for the windows [m]
 par.room.winNums = 5; % number of windows
 par.room.winWall = 1; % wall number where the windows are placed
 par.room.roofWall = 5; % wall with the transmission coefficient par.room.Troof
 par.bandwidth = 5e6;   % bandwidth used in [Hz]
 par.del_spread = 50e-9;  % delay spread in [s]
  % define electromagnetic properties and rec. / trans. positions
 par.freq = 2e9;  % carrier frequency in [Hz]
 par.room.Twin = 10^(-5/20);   % transmission factor for windows [linear]
 par.room.Twall = 10^(-11/20); % transmission factor for walls [linear] 
 par.room.Troof = 10^(-23/20); % transmission factor for roof [linear]
 walls = createRoom(par.room);  % create the structure of the room
%  gridRec = c/freq;         % grid for calculation of gamma0
 gridRec = gridSize;         % grid for calculation of gamma0
 elev = 25:1:85;   % elevation in [deg] for the transmitter
 azi = 0;          % azimuth in [deg] for the transmitter
 trans = 20e6*[cosd(azi)*cosd(elev); cosd(elev)*sind(azi); sind(elev)]; % transmitter positions
 rec = [par.room.dx/2-1,0,-par.room.dz/2+1].';  % receiver position within the room [1m before window, 1m height]
  % perform the simulation
 normLoS = false; % set normalisation to LoS to false for all subroutines
 P = zeros(1,length(trans));  % init output
 fprintf('Checking against Milan''s model with gridSize=%g m \n',gridRec);
%  tic;
 for i = 1:size(trans,2) 
  fprintf('  Calc. at %i/%i',i,size(trans,2)); tic;
  gamma0 = calcPowerDelay(trans(:,i),walls,par.freq,par.bandwidth,gridRec);  % calculate gamma0
  E = calcField(trans(:,i),rec,walls,par.freq,par.bandwidth,normLoS);  % calculate field contribution from all walls
  P(i) = calcNarrowBandPower(E,gamma0,par.del_spread,par.bandwidth);  % calculate narrowband power
  % P(i) is the power including free space loss, therefore normalise it
  P(i) = P(i)*norm(trans(:,i)-rec)^2;  % normalise to free space loss as P/Pfsl
  fprintf('   needed %g s\n',toc);
 end;
%  fprintf('Overall simulation took %g s\n',toc);
  % create figure; 
 loss = -10*log10(P);  % loss calculation in terms of narrowband power
 p = polyfit(elev,loss,1);   % fit linear regression
 fprintf('LOS case : got B=%g A=%g\n\n',p(1),p(2));
 figure; plot(elev,loss,'b',elev,p(1)*elev+p(2),'r'); grid on;
 title('LOS case');
 xlabel('Elevation in [deg]'); ylabel('Loss in [dB]');
 plotGeometry(walls,rec,trans);
 keyboard;
return; 

function P = calcNarrowBandPower(E,gamma0,spread,bandwidth) 
% Function to calculate the narrowband power from the entry loss.
% E is the actual entry loss at the receiver position, while 
% gamm0 and spread describe the exp-PDP
 P = 0.5*gamma0*(sqrt(1+4*bandwidth^2*spread^2)-1)+abs(E).^2; % using exp. PDP
 % replacing the value gamma0 at tau=0 by abs(E)^2
return;


% ---- Additional subroutines for helping etc. ---- %

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
 for i = 1:size(rx,2)
  t = 10*tx(:,i)/norm(tx(:,i));
  plot3(rx(1,i)+[0,t(1)],rx(2,i)+[0,t(2)],rx(3,i)+[0,t(3)],'-b');
 end;
return;


function fig = plotRoom2d(walls,fig)
% function to plot the room geometry in 2d (x,y-plane)
 color = 'grcmyk';
 if ~exist('fig','var')
  fig = figure;  % create a new figure
 else
  figure(fig);   % set the figure 
 end;
 hold on; grid on;
 plot(0,0,'ob');  % plot a point in the origin
 for i = 1:length(walls)  % for each wall
  n = cross(walls(i).v1,walls(i).v2);  % wall normal
  if abs(n(3))~=norm(n)
    % plot only those walls where the normals are not directed via z
   plot([walls(i).s(1),walls(i).s(1)+walls(i).v1(1)], [walls(i).s(2),walls(i).s(2)+walls(i).v1(2)],color(i),'LineWidth',2);  % plot along v1 direction
   plot([walls(i).s(1),walls(i).s(1)+walls(i).v2(1)], [walls(i).s(2),walls(i).s(2)+walls(i).v2(2)],color(i),'LineWidth',2);  % plot along v2 direction
   disp(sprintf(' Wall %i : color is %c',i,color(i)));
  end;
 end;
 axis equal;
 xlabel('x-axis [m]'); ylabel('y-axis [m]'); zlabel('z-axis [m]');
return;


function fig = plotScatt2d(pScatt,limits,color,fig)
% Function to plot the scatterers in a 2d plot including the 
% opening angle. The figure fig is used if provided.
% pScatt is one scatterer per row.
 len = 1;   % default len for opening cone
 if ~exist('color','var') || isempty(color)
  color = 'b';   % default value
 end;
 if ~exist('fig','var')
  fig = figure;  % create a new figure
 else
  figure(fig);   % set the figure 
 end; 
 plot(pScatt(:,1),pScatt(:,2),strcat('.',color));  % plot scatterers
 for i = 1:size(limits,1)
  plot([pScatt(i,1),pScatt(i,1)+len*cosd(limits(i,1))],[pScatt(i,2),pScatt(i,2)+len*sind(limits(i,1))],strcat('-',color));  % plot scatterers
  plot([pScatt(i,1),pScatt(i,1)+len*cosd(limits(i,2))],[pScatt(i,2),pScatt(i,2)+len*sind(limits(i,2))],strcat('-',color));  % plot scatterers
 end;
return;


function plotpdf(x,bins,supKS,func)
% Function to plot an estimate of the pdf for the data values x.
% bins are the center for the histogram as far as provided.
% supKS is the support value for the ksdensity estimation 
% function of Matlab. 
 if ~exist('bins','var') || isempty(bins)
  bins = linspace(min(x),max(x),length(x)/10);
 end;
 dbin = bins(2)-bins(1);  % bin size
 x = x(:);
 n = histc(x,bins)/(length(x)*dbin);  % histogram estimate
 bar(bins,n); grid on; hold on;
 if exist('supKS','var') && ~isempty(supKS)
  [f,xi] = ksdensity(x,bins(1):dbin/5:bins(end),'support',supKS);
 else
  [f,xi] = ksdensity(x,bins(1):dbin/5:bins(end));
 end;
 plot(xi,f,'g.-','LineWidth',2);
 if (~exist('func','var')) || isempty(func)
  plot(xi,normpdf(xi,mean(x),std(x)),'b.-','LineWidth',2);
  legend('Histogram','Kernel Estimator','Gauss');
 else
  plot(xi,func(xi),'r-','LineWidth',2);
  legend('Histogram','Kernel Estimator','Data Fit');
 end;
return


function plotcdf(x,supKS,func)
% Function to plot an estimate of the cdf for the data values x.
% supKS is the support value for the ksdensity estimation 
% function of Matlab. func is a function pointer to the cdf
% to be displayed instead of a Gaussian cdf.
 %cdfplot(x); hold on;
 [f,x_val,flo,fup] = ecdf(x);  % calculate empirical cdf with lower and upper bounds
 stairs(x_val,f,'LineWidth',2); hold on; grid on;  % plot empirical cdf
 stairs(x_val,flo,'r:','LineWidth',2);  % lower bound
 stairs(x_val,fup,'r:','LineWidth',2);  % upper bound 
 xval = min(x):(max(x)-min(x))/1000:max(x);  % x-axis to be plotted
 if exist('supKS','var') && ~isempty(supKS)
  [f,xi] = ksdensity(x,xval,'support',supKS,'function','cdf');
 else
  [f,xi] = ksdensity(x,xval,'function','cdf');
 end;
 plot(xi,f,'g-','LineWidth',2);
 if (~exist('func','var')) || isempty(func)
  plot(xi,normcdf(xi,mean(x),std(x)),'b-','LineWidth',2);  % plot the standard Gaussian
  legend('Emp. CDF','ECDF lower bound 95%','ECDF upper bound 95%','Kernel Estimator','Gauss');
 else
  plot(xi,func(xi),'r-','LineWidth',2);
  legend('Emp. CDF','ECDF lower bound 95%','ECDF upper bound 95%','Kernel Estimator','Data Fit');
 end;
return



%% ----- Different PDFs for Testing ----- %%

function f = mixUniMisesPDF(x,p,mu,kap)
% Function to calculate the pdf of a mixture between 
% a uniform distribution [0,360] and a von Mises distribution
% with parameters mu (mean) and kap as
% f = p*U(0,360)+(1-p)*M(x|mu,kap). x is understood as angle
% in degree from 0..360.
 % for testing
%  sum(mixUniMisesPDF(-180:1e-1:180,par(1),par(2),par(3)))*1e-1
 x = x(:); p=p(1); mu=mu(1); kap=kap(1);
%  x = mod(x,360);    % take care that angles are within 0..2*pi
 bess = besseli(0,kap,1); % with third parameter is I_0(k)*e^{-k}
 if (isnan(bess) || isinf(bess))
  disp('error');
  keyboard;
 end;
 C = 1/(360*bess);  % pre-factor of von Mises pdf
 f = p/(360)+(1-p)*C*exp(kap*(cosd(x-mu)-1));
return;

function f = truncNormPDF(x,mu,sig)
% Function to calculate the pdf of a truncated normal pdf
% with K*norm(x|mu,sig)*u(x) with u(x)=1 for x>=0.
 K = 1/(1-normcdf(0,mu,sig));   % pre-factor to make integral to 1
 f = K/(sig*sqrt(2*pi))*exp(-(1/(2*sig^2))*(x-mu).^2);
 f(x<0) = 0;  % truncation
return;

function f = truncExpPDF(x,mu,lim)
% Function to calculate the pdf of a truncated exp-pdf
% with K*exp(-mu*x) within [lim(1)...lim(2)].
 mu = mu(1);
 lim = sort(lim); lim = lim(1:2);
 K = mu/(exp(-mu*lim(1))-exp(-mu*lim(2)));  % pre-factor for integral to be 1
 f = K*exp(-mu*x);
 f((x<lim(1))|(x>lim(2))) = 0;  % truncation
return;

function f = truncWeibullPDF(x,lam,k)
% Function to calculate the pdf of a truncated Weibull 
% distribution, the truncation will be [0,1]
 K = 1/(1-exp(-(1/lam)^k));  % pre-factor to normalise integral to one
 f = K*(k/lam)*(x/lam).^(k-1).*exp(-(x/lam).^k);  % truncated Weibull
 f(x<0) = 0;  % do truncation for limits
 f(x>1) = 0;
return;

function f = truncWeibullCDF(x,lam,k)
% Function to calculate the cdf of a truncated Weibull 
% distribution, the truncation will be [0,1]
 K = 1/(1-exp(-(1/lam)^k));  % pre-factor to normalise integral to one
 f = K*(1-exp(-(x/lam).^k));  % truncated Weibull
 f(x<0) = 0;  % do truncation for limits
 f(x>1) = 1;
return;

function f = invGaussPDF(x,mu,lam)
% Function to calculate the pdf of an inverse Gaussian distribution
%  see http://en.wikipedia.org/wiki/Inverse_Gaussian_distribution
 f = zeros(size(x));
 f(x>0) = sqrt(lam./(2*pi*x(x>0).^3)).*exp(-lam*(x(x>0)-mu).^2./(2*mu^2*x(x>0))); 
return;

function f = negBinomialPDF(x,r,p)
% Function to calculate the pdf of a negative binomial distribution.
% The geometric distribution is a special case with negBinomialPDF(x,1,p).
 f = zeros(size(x));
 k = find(x>=0);     % only for values x>=0
 f(k) = exp(gammaln(r+x(k)) - gammaln(x(k)+1) - gammaln(r) + r*log(p) + x(k)*log(1-p));
return;





