function walls = createRoom(par,show)
% Function to create a single room for simulation.
%
%   walls = createRoom(par)
% 
% Input: par is a struct of certain parameters.
%         par.dx = x-dimension [m]
%         par.dy = y-dimension [m]
%         par.dz = z-dimension [m] (height)
%         par.winWidth = width of windows [m]
%         par.winUpper = distance to upper edge of the wall for the windows [m]
%         par.winLower = distance to lower edge of the wall for the windows [m]
%         par.winNums = number of windows on par.winWall
%         par.winWall = wall index number where the windows are placed
%         par.roofWall = wall with the transmission coefficient par.Troof
%         par.Twin = transmission factor for windows [linear]
%         par.Twall = transmission factor for walls [linear] 
%         par.Troof = transmission factor for roof [linear]
%        show = for testing to show the room as figure if show=true (optional)
% Output: walls is a struct as room definition
%          walls(nr).s = lower left point of the wall seen from inside
%          walls(nr).v1 = vector of first axis of the wall
%          walls(nr).v2 = vector of second axis of the wall
%          walls(nr).n = normal vector of the wall
%          walls(nr).k = affine transformation parameters (k.A, k.b) for
%                          the wall
%          walls(nr).Twall = transmission coefficient (linear) for the wall
%          walls(nr).win = struct for the windows
%           walls(nr).win.s,walls(nr).win.v1,walls(nr).win.v2,walls(nr).win.n,walls(nr).win.k = as for the wall itself
%           walls(nr).win.Twin = transmission coefficient (linear) for the
%                                 window
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

 if ~exist('par','var') | ~isstruct(par)
  error('function createRoom: Please check the par parameter before calling this function!');
 end;

  % create the walls for the room 
 walls = createWalls(par);  % create the wall plains 
 walls = addWindows(par,walls);  % place windows
 
 if exist('show','var') & islogical(show) & show
  plotRoom(walls);
 end;

return;



%% ----- Subroutines to Define the Room ----- %%

function wall = addWindows(par,wall)
% Function to add windows to a wall wall. par.winUpper (default=0.1m) is the 
% distance to the upper edge, while par.winLower (default=0.75m) is the distance
% to the lower edge. par.winNums is the number of windows and
% par.winWidth the window size. Upper and lower is defined by 
% vector wall.v2, assuming wall.v2 going upwards on the wall.
% wallNr is the number of the wall, where the windows will be added.
 if ~exist('par','var') || ~isstruct(par) || ~getfield(par,'winNums') || ~getfield(par,'winWidth') || ~exist('wall','var') || ~isstruct(wall)
  error('In routine addWindows: Parameters wrong!');
 end;
 if ~getfield(par,'winUpper')
  upper = 0.1;
 else
  upper = par.winUpper;
 end;
 if ~getfield(par,'winLower')
  lower = 0.75;
 else
  lower = par.winLower;
 end;
  % lower is now the distance on the wall to the lower edge of the wall
  % upper is now the distance on the wall to the upper edge of the wall
 wallNr = par.winWall;  % wall where the windows are placed on
 wallLen = norm(wall(wallNr).v1);  % wall length
 wallHeight = norm(wall(wallNr).v2);  % height of the wall
 if (wallLen<par.winWidth*par.winNums)
  error('In routine addWindows: Too many windows, wall is too small!');
 end;
 dist = (wallLen-par.winNums*par.winWidth)/par.winNums;  % space between the windows
 for i = 1:par.winNums  % now place the windows
  wall(wallNr).win(i).s = wall(wallNr).s+(dist/2+(i-1)*(par.winWidth+dist))*unitVec(wall(wallNr).v1)+lower*unitVec(wall(wallNr).v2); % left lower corner of the window
  wall(wallNr).win(i).v1 = par.winWidth*unitVec(wall(wallNr).v1);  % vector v1 to define window plane
  wall(wallNr).win(i).v2 = (wallHeight-upper-lower)*unitVec(wall(wallNr).v2);  % vector v1 to define window plane
  wall(wallNr).win(i).n = wall(wallNr).n;  % sharing the same normal
  wall(wallNr).win(i).k = wall(wallNr).k;  % sharint the same transformation
  wall(wallNr).win(i).Twin = par.Twin;  % set the transmission factor
 end;
return;

function walls = createWalls(par)
% Function calculates the walls given the room size par.dx,par.dy and par.dz.
% Directions are from inside the room. Normal vector n is pointing
% outside the room. par.Twall is the transmission coefficient of walls.
% par.roofWall is a special wall with transmission coefficient 
% par.Troof.
  % wall on the front side of the room (y,z plane, x constant)
 walls(1).s = [par.dx/2,par.dy/2,-par.dz/2].';  % reference point to the plain (lower left corner)
 walls(1).v1 = [0,-par.dy,0].'; % first vector of plain
 walls(1).v2 = [0,0,par.dz].';  % second vector of plain
 walls(1).n = -unitVec(cross(walls(1).v1,walls(1).v2));
 walls(1).k = coorTransform([0,0,0; 1,0,0; 0,1,0; 0,0,1].',...
    repmat(walls(1).s,1,4)+[zeros(3,1),walls(1).n,unitVec(walls(1).v1),unitVec(walls(1).v2)]);
 walls(1).Twall = par.Twall;  % set the transmission factor   
  % wall at the right of the room (x,z plane, y constant)
 walls(2).s = [par.dx/2,-par.dy/2,-par.dz/2].';  % reference point to the plain (lower left corner)
 walls(2).v1 = [-par.dx,0,0].'; % first vector of plain
 walls(2).v2 = [0,0,par.dz].';  % second vector of plain
 walls(2).n = -unitVec(cross(walls(2).v1,walls(2).v2));
 walls(2).k = coorTransform([0,0,0; 1,0,0; 0,1,0; 0,0,1].',...
    repmat(walls(2).s,1,4)+[zeros(3,1),walls(2).n,unitVec(walls(2).v1),unitVec(walls(2).v2)]);
 walls(2).Twall = par.Twall;  % set the transmission factor      
  % wall at the back side of the room (y,z plane, x constant)
 walls(3).s = [-par.dx/2,-par.dy/2,-par.dz/2].';  % reference point to the plain (lower left corner)
 walls(3).v1 = [0,par.dy,0].'; % first vector of plain
 walls(3).v2 = [0,0,par.dz].'; % second vector of plain
 walls(3).n = -unitVec(cross(walls(3).v1,walls(3).v2));
 walls(3).k = coorTransform([0,0,0; 1,0,0; 0,1,0; 0,0,1].',...
    repmat(walls(3).s,1,4)+[zeros(3,1),walls(3).n,unitVec(walls(3).v1),unitVec(walls(3).v2)]);
 walls(3).Twall = par.Twall;  % set the transmission factor      
  % wall at the left side of the room (x,z plane, y constant)
 walls(4).s = [-par.dx/2,par.dy/2,-par.dz/2].';  % reference point to the plain (lower left corner)
 walls(4).v1 = [par.dx,0,0].'; % first vector of plain
 walls(4).v2 = [0,0,par.dz].'; % second vector of plain  
 walls(4).n = -unitVec(cross(walls(4).v1,walls(4).v2));
 walls(4).k = coorTransform([0,0,0; 1,0,0; 0,1,0; 0,0,1].',...
    repmat(walls(4).s,1,4)+[zeros(3,1),walls(4).n,unitVec(walls(4).v1),unitVec(walls(4).v2)]); 
 walls(4).Twall = par.Twall;  % set the transmission factor      
  % roof (x,y plane, z constant)
 walls(5).s = [par.dx/2,par.dy/2,par.dz/2].';  % reference point to the plain (lower left corner)
 walls(5).v1 = [0,-par.dy,0].'; % first vector of plain  
 walls(5).v2 = [-par.dx,0,0].'; % second vector of plain
 walls(5).n = -unitVec(cross(walls(5).v1,walls(5).v2));
 walls(5).k = coorTransform([0,0,0; 1,0,0; 0,1,0; 0,0,1].',...
    repmat(walls(5).s,1,4)+[zeros(3,1),walls(5).n,unitVec(walls(5).v1),unitVec(walls(5).v2)]); 
 walls(5).Twall = par.Twall;  % set the transmission factor      
  % floor (x,y plane, z constant)
 walls(6).s = [-par.dx/2,par.dy/2,-par.dz/2].';  % reference point to the plain (lower left corner)
 walls(6).v1 = [0,-par.dy,0].'; % first vector of plain
 walls(6).v2 = [par.dx,0,0].';  % second vector of plain   
 walls(6).n = -unitVec(cross(walls(6).v1,walls(6).v2));
 walls(6).k = coorTransform([0,0,0; 1,0,0; 0,1,0; 0,0,1].',...
    repmat(walls(6).s,1,4)+[zeros(3,1),walls(6).n,unitVec(walls(6).v1),unitVec(walls(6).v2)]);  
 walls(6).Twall = par.Twall;  % set the transmission factor      
 walls(par.roofWall).Twall = par.Troof;  % set special roof transmission coefficient
%   % verification of values                           
%  plotRoom(walls);
%  i = 5;   % wall under test
%  x = doTransform(walls(i).s,walls(i).k,false);  % check for x(1)==0 ??
%  disp(sprintf('Got       s = %s',mat2str(x(:).',3)));
%  x = doTransform(walls(i).s+walls(i).v1,walls(i).k,false);  % check for x(1)==0 ??
%  disp(sprintf('Got    s+v1 = %s -> in positive y-Richtung',mat2str(x(:).',3)));
%  x = doTransform(walls(i).s+walls(i).v2,walls(i).k,false);  % check for x(1)==0 ??
%  disp(sprintf('Got    s+v2 = %s -> in positive z-Richtung',mat2str(x(:).',3)));
%  x = doTransform(walls(i).s+walls(i).v2+walls(i).v1,walls(i).k,false);  % check for x(1)==0 ??
%  disp(sprintf('Got s+v1+v2 = %s -> in positive y & z-Richtung',mat2str(x(:).',3)));
%  x = doTransform(zeros(3,1),walls(i).k);
%  disp(sprintf('Got %s as origin dist to s in old=%g in new=%g',mat2str(x(:).',3),...
%        norm(walls(i).s),norm(x-doTransform(walls(i).s,walls(i).k,false))));
%  keyboard;
return; 


%% ----- Supplementary Functions ----- %%

function fig = plotRoom(walls,fig)
% function to plot the room geometry.
 color = 'grcmyk';
 if ~exist('fig','var')
  fig = figure;  % create a new figure
 else
  figure(fig);   % set the figure 
 end;
 hold on; grid on;
 plot3(0,0,0,'ob');  % plot a point in the origin
 for i = 1:length(walls)  % for each wall
  point = [walls(i).s(:).';
           walls(i).s(:).'+walls(i).v2(:).';
           walls(i).s(:).'+walls(i).v2(:).'+walls(i).v1(:).';
           walls(i).s(:).'+walls(i).v1(:).'];
  patch(point(:,1),point(:,2),point(:,3),color(i),'FaceAlpha',0.2);
  midPoint = walls(i).s+(walls(i).v1+walls(i).v2)/2; % middle point of plane
  quiver3(midPoint(1),midPoint(2),midPoint(3),...
          midPoint(1)+walls(i).n(1),midPoint(2)+walls(i).n(2),midPoint(3)+walls(i).n(3),color(i));
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


function k = coorTransform(xnew,x)
% function provides a coordinate transformation
% for vectors as coloumns in xpi (new coordinate system)
% and xi (old coordinate system. Output will be k.A and
% k.b with the linear relation xpi=k.A*xi+k.b.
% Calculated as
% [xpi;ones(1,size(xpi,2)]=[A,b;0...0,1]*[xi;ones(1,size(xi,2)].
 xnew = [xnew;ones(1,size(xnew,2))];
 x = [x;ones(1,size(x,2))];
 if (abs(det(x))<3*eps)
  warning('In coorTransform: Matrix badly scaled!');
 end;
 A = xnew/x;
 A(abs(A)<10*eps)=0;  % set low numbers to zero
 k.A = A(1:end-1,1:end-1);
 k.b = A(1:end-1,end);
return;


function l = unitVec(l)
% Function normalises the vector l to a length of 1
 l = l/norm(l);
return;