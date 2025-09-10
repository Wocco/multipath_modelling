Satellite-to-Indoor Channel Simulator

In indoor environments the accuracy of positioning by global navigation 
satellite systems suffers significantly from signal blockage, reflection 
and diffraction. To develop advanced receiver position algorithms working 
in harsh propagation environments, accurate channel simulators are 
necessary. Therefore, DLR investigated the satellite-to-indoor propagation 
channel in its very detail by a measurement campaign in 2008. As an 
outcome, we proposed a novel and accurate wideband satellite-to-indoor 
channel model. Compared to the state of the art, the proposed channel 
model is able to reproduce the spatial characteristics of the wideband 
propagation channel for a moving receiver. The model is based on a hybrid 
approach combining physical-deterministic and stochastic methods: 
(a) waves diffracted and penetrated by walls, windows and doors are 
considered by using physical-deterministic near field methods, 
(b) in order to model multipath components occurring due to reflections 
on walls, a hybrid approach is used, 
(c) the behavior of scattered waves is stochastically modeled. 
The proposed channel model accurately models satellite-to-indoor 
propagation effects, and, thus, can be used for testing and validating 
range estimation algorithms for positioning.


Table of Contents
--------------------------------------------------------------------------
1. Installation
2. Example Simulation File
3. Input Parameters
   3.1. Initialization Phase
   3.2. Runtime Phase
4. Output Parameters
   4.1. Initialization Phase
   4.2. Runtime Phase
5. Further Documentation
6. Contact
7. Literature


1. Installation
--------------------------------------------------------------------------
The satellite-to-indoor channel simulator can be copied to any directory. Please make sure that either all files are within the same directory or the Matlab path variable is set accordingly. In case that the program shall be run on any other systems than Windows 64 Bit or Linux 64 Bit Matlabs, then the mex-files window_integ.cpp and kden.cpp need to be compiled for the specific system. 
Please run (within Matlab command line):

mex window_integ.cpp
mex kden.cpp

to create the compiled versions. Please make sure that a C++ compiler is available within the Matlab environment. The mex-files have been compiled using Microsoft Visual Studio and GCC.


2. Example Simulation File
--------------------------------------------------------------------------
The Matlab file simu.m shows an example of how the satellite-to-indoor channel simulator can be used. Please have a look.


3. Input Parameters
--------------------------------------------------------------------------
The satellite-to-indoor channel simulator is splitted into two parts. First an initialization has to be performed. In a second step that can be arbitrary often repeated, the channel impulse response for a specific receiver position is calculated. Please see also the help in sim_chan.m and createRoom.m for further information.

3.1. Initialization Phase
The input parameters for the "Initialization Phase" to the program sim_chan.m are:

 par.bandwidth - Basic bandwidth in [Hz] of the simulation that is used to define the sampled spatial power delay profile. Please be aware that there is no sampling of the multipath in terms of delay is performed, i.e. multipath components have arbitrary real valued delays.
 
 par.del_spread - Delay spread [s] to define the slope of the exponential decaying spatial power delay profile.
 
 par.freq - Carrier frequency in [Hz] to be used for the simulations. Please be aware that the model has been validated only for L-band frequencies.
 
 par.trans_pos - Three dimensional transmitter position in [m] in the Cartesian coordinate system defined by the scenery. The parameter shall be given as [x,y,z] coordinates.
 
 par.Nscatt - Optional parameter allowing to vary the number of scattered components to be simulated. The default value is 10000.
 
 walls - Struct to define the scenery of the simulated environment. For a complete definition see the help of the output of createRoom.m. Otherwise, create the parameter walls as done in simu.m by the help of createRoom.m using:
 
 room.dx = 10;  % x-dimension of the room [m]
 room.dy = 10;  % y-dimension of the room [m]
 room.dz = 3;   % z-dimension of the room [m] (height)
 room.winWidth = 1;  % window width [m]
 room.winUpper = 0.3; % distance to upper edge for the windows [m]
 room.winLower = 1.2; % distance to lower edge for the windows [m]
 room.winNums = 5;  % number of windows on wall room.winWall
 room.winWall = 1;  % wall number where the windows are placed
 room.roofWall = 5; % wall with the transmission coefficient room.Troof
 room.Twin = 10^(-2/20);   % perpendicular transmission factor for windows [linear]
 room.Twall = 10^(-15/20); % perpendicular transmission factor for walls [linear] 
 room.Troof = 10^(-25/20); % perpendicular transmission factor for roof [linear]
 
 walls = createRoom(par.room);  % use script to create the parameter walls for sim_chan
  
 
3.2. Runtime Phase
During the runtime phase, the satellite-to-indoor channel simulator can be called repeatedly in order to calcualte channel impulse responses for specific receiver positions. Therefore, the input parameters during this phase are:

 obj - Struct as output of sim_chan.m during the initialization phase. No user modification are needed.
 
 rx - Three dimensional receiver position in [m] in the Cartesian coordinate system defined by the scenery. The parameter shall be given as [x,y,z] coordinates.
 

4. Output Parameters
--------------------------------------------------------------------------
The satellite-to-indoor channel simulator is splitted into two parts. First an initialization has to be performed. In a second step that can be arbitrary often repeated, the channel impulse response for a specific receiver position is calculated. Please see also the help in sim_chan.m and createRoom.m for further information.

4.1. Initialization Phase
The output parameter for the "Initialization Phase" of the program sim_chan.m is:

 obj - Struct that defines internal parameters of the satellite-to-indoor channel simulator. No user modification are needed.

 
4.2. Runtime Phase
During the runtime phase, the satellite-to-indoor channel simulator can be called repeatedly in order to calcualte channel impulse responses for specific receiver positions. Therefore, the output parameters during this phase are:

 cir - Channel impulse related parameters as matrix. Each row i defines a multipath component (MPC). 
        cir(i,obj.info.delay) - Delay in [s] of the MPC normalized to the transmitter-receiver direct propagation time.
        cir(i,obj.info.amp) - Complex amplitude of the MPC as electrical field strength.
        cir(i,obj.info.typeCol) - Type of MPC, see obj.info.type.description{cir(i,obj.info.typeCol)} for more information.
        
 p - Coordinate of the last interaction point for the MPC. One path per row i, provided as [x,y,z] in meters in the Cartesian coordinate system of the scenery.


5. Further Documentation
--------------------------------------------------------------------------
Matlab functions are documented according to Matlab standard, please use "help <function>" or "doc <function>". The included mex-files will provide help if called without any parameter.


6. Contact
--------------------------------------------------------------------------
Thomas Jost
Deutsches Zentrum für Luft- und Raumfahrt (DLR)
Institut für Kommunikation und Navigation, Nachrichtensysteme
Oberpfaffenhofen-Wessling
thomas.jost@dlr.de

The software can also be downloaded under GPL v2.0 at:
http://www.dlr.de/kn/channel_models


7. Literature
--------------------------------------------------------------------------
The following publications describe the modeling approach for the satellite-to-indoor propagation channel:

- T. Jost, W. Wang, U.-C. Fiebig, and F. Pérez-Fontán
  A Wideband Satellite-to-Indoor Channel Model for Navigation Applications
  IEEE Trans. Antennas Propag., vol. 62, no. 10, pp. 5307-5320, Oct. 2014
- T. Jost, G. Carrié, F. Pérez-Fontán, W. Wang, and U.-C. Fiebig
  A Deterministic Satellite-to-Indoor Entry Loss Model 
  IEEE Trans. Antennas Propag., vol. 61, no. 4, pp. 2223-2230, Apr. 2013
- ITU Report P.2145-2, 09/2017 
  "Model parameters for the physical-statistical wideband models in Recommendation ITU-R P.681"
- T. Jost, Satellite-to-Indoor Wave Propagation for Positioning Applications, 
  ser. Kommunikationstechnik. Verlag Dr. Hut, 2014, PhD dissertation University of Vigo

The further publications describe the measurement campaign and the data analysis:

- T. Jost, W. Wang, U.-C. Fiebig, and F. Pérez-Fontán
  Detection and Tracking of Mobile Propagation Channel Paths
  IEEE Trans. Antennas Propag., vol. 60, no. 10, pp. 4875-4883, Oct. 2012 
- W. Wang and T. Jost
  A Low-Cost Platform for Time-Variant Wireless Channel Measurements With Application to Positioning
  IEEE Trans. Instrum. Meas., vol. 61, no. 6, pp. 1597-1604, June 2012 
- T. Jost, W. Wang, U.-C. Fiebig, and F. Pérez-Fontán
  Movement of Equivalent Scatterers in Geometry-Based Stochastic Channel Models
  IEEE Antennas Wireless Propag. Lett., vol. 11, pp. 555-558, May 2012
- T. Jost, W. Wang, U.-C. Fiebig, and F. Pérez-Fontán
  Comparison of L- and C-Band Satellite-to-Indoor Broadband Wave Propagation for Navigation Applications
  IEEE Trans. Antennas Propag., vol. 59, no. 10, pp. 3899-3909, Oct. 2011 



