function data = static_receiver(elev,azi)
%STATIC_RECEIVER Simulate channel for a static receiver with varying angles.
%   DATA = STATIC_RECEIVER(ELEV,AZI) computes the channel impulse response
%   for a receiver at a fixed position while the transmitter elevation and
%   azimuth angles change over time.  ELEV and AZI must be vectors of equal
%   length containing the angles in degrees for each time instant.  The
%   function returns a struct array DATA where each element contains the
%   delay, amplitude and path type of the channel for the corresponding
%   time instant.
%
%   Example:
%       elev = linspace(20,40,5);
%       azi  = linspace(80,120,5);
%       data = static_receiver(elev,azi);
%
%   The receiver position is fixed at 1 m height and 1 m in front of the
%   window centre, while the satellite (transmitter) moves according to the
%   provided angles.
%
%   This is a helper function derived from dop_test.m for scenarios where
%   the receiver does not move but the satellite does.

% Copyright 2016 Thomas Jost, German Aerospace Center (DLR)
%
% This file is part of the satellite-to-indoor channel simulator.
%
% The satellite-to-indoor channel simulator is free software: you can
% redistribute it and/or modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation, either version 2 of
% the License, or (at your option) any later version.
%
% The satellite-to-indoor channel simulator is distributed in the hope that
% it will be useful, but WITHOUT ANY WARRANTY; without even the implied
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with the program.  If not, see <http://www.gnu.org/licenses/>.

if nargin~=2
    error('static_receiver: Please provide elevation and azimuth vectors.');
end
if numel(elev) ~= numel(azi)
    error('static_receiver: ELEV and AZI must have the same length.');
end

% --- Room and simulation parameters ---
par.room.dx = 10;
par.room.dy = 10;
par.room.dz = 3;
par.room.winWidth = 1;
par.room.winUpper = 0.3;
par.room.winLower = 1.2;
par.room.winNums = 5;
par.room.winWall = 1;
par.room.roofWall = 5;
par.bandwidth = 100e6;
par.del_spread = 30e-9;
par.freq = 1.5e9;
par.room.Twin = 10^(-2/20);
par.room.Twall = 10^(-15/20);
par.room.Troof = 10^(-25/20);
par.Nscatt = 10000;

% Receiver position (1 m before window, 1 m height)
rec = [3;0;(-par.room.dz/2+1)];

% Prepare room once
walls = createRoom(par.room);

numT = numel(azi);

data(numT) = struct('delay',[],'amp',[],'type',[]);
    disp(numT)
for t = 1:numT
    % Transmitter position for current angles
    disp(t)
    
    par.trans_pos = 20e6*[cosd(azi(t))*cosd(elev(t)); ...
                          cosd(elev(t))*sind(azi(t)); ...
                          sind(elev(t))];

    % Initialise channel with current transmitter position
    obj = sim_chan(par,walls);

    % Compute channel for static receiver
    cir = sim_chan(obj,rec);

    data(t).delay = cir(:,obj.info.delay);
    data(t).amp   = cir(:,obj.info.amp);
    data(t).type  = cir(:,obj.info.typeCol);
end

end