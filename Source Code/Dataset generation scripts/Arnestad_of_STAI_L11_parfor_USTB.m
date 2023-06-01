%% Computation of a STAI dataset with Field II to simulate specular tissue using parameters of an L11-4v 128 element Verasonics Transducer and beamforming with USTB and using MATLAB's parfor function for parallelization
%
% This example shows how to load the data from a Field II simulation of specular tissue into
% USTB objects, and then beamform it with the USTB routines.
% This example uses the L11-4v 128 element Verasonics Transducer
% The Field II simulation program (<field-ii.dk>) should be in MATLAB's path.
%
% This tutorial assumes familiarity with the contents of the
% <../../fresnel/linear_array/html/CPWC_linear_array.html 'CPWC simulation with the USTB built-in Fresnel
% simulator'> tutorial. Please feel free to refer back to that for more
% details.
%
% _by Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>, Ole Marius Hoel
% Rindal <olemarius@olemarius.net> and Arun Asokan Nair <anair8@jhu.edu> 09.05.2017_

%% Clear old workspace and close old plots

clear all;
close all;
% Dette må endres mellom brukerne
if ispc
    addpath(genpath('\\hume.uio.no\student-u55\helenewo\pc\downloads\Field_II_ver_3_30_windows.tar'));
    addpath(genpath('\\hume.uio.no\student-u55\helenewo\pc\Dokumenter\USTB'));
    addpath(genpath('\\hume.uio.no\student-u55\helenewo\MasterThesis'));
elseif isunix
    addpath(genpath('/hom/dsb/field'));
    addpath(genpath('/uio/hume/student-u55/helenewo/pc/Dokumenter/USTB'));
    addpath(genpath('/uio/hume/student-u55/helenewo/MasterThesis/MasterThesis'));
    addpath(genpath('/uio/hume/studentg-u55/helenewo/MasterThesis/datasets'));
end

filename = 'specklesim_v2.uff';
data_path = fullfile('..','datasets');

number_of_scatterers = 1e6;


p = gcp('nocreate'); % If no pool, do not create new one.

if isempty(p)
%     parpool
%     parpool(16)
    myCluster=parcluster('local');
    myCluster.NumWorkers=58;
    parpool(myCluster,myCluster.NumWorkers)
% else
%     poolsize = 16;
end

%% Basic Constants
%
% Our first step is to define some basic constants for our imaging scenario
% - below, we set the speed of sound in the tissue, sampling frequency and
% sampling step size in time.

c0=1540;     % Speed of sound [m/s]
fs=100e6;    % Sampling frequency [Hz]
dt=1/fs;     % Sampling step [s]

%% field II initialisation
%
% Next, we initialize the field II toolbox. Again, this only works if the
% Field II simulation program (<field-ii.dk>) is in MATLAB's path. We also
% pass our set constants to it.

field_init(0);
set_field('c',c0);              % Speed of sound [m/s]
set_field('fs',fs);             % Sampling frequency [Hz]
set_field('use_rectangles',1);  % use rectangular elements

%% Transducer definition L11-4v, 128-element linear array transducer
%
% Our next step is to define the ultrasound transducer array we are using.
% For this experiment, we shall use the L11-4v 128 element Verasonics
% Transducer and set our parameters to match it.

probe = uff.linear_array();
f0                      = 2.56e6;      % Transducer center frequency [Hz]
lambda                  = c0/f0;           % Wavelength [m]
probe.element_height    = 5e-3;            % Height of element [m]
probe.pitch             = 0.300e-3;        % probe.pitch [m]
kerf                    = 0.05e-03;        % gap between elements [m]
probe.element_width     = probe.pitch-kerf;% Width of element [m]
lens_el                 = 60e-3;           % position of the elevation focus
probe.N                 = 64;             % Number of elements %CHANGED
pulse_duration          = 2.5;             % pulse duration [cycles]
z_focus = 40/1000;  %CHANGED

%% pulse definition
pulse = uff.pulse();
pulse.center_frequency = f0;
pulse.fractional_bandwidth = 0.65;        % probe bandwidth [1]
t0 = (-1/pulse.fractional_bandwidth/f0): dt : (1/pulse.fractional_bandwidth/f0);
impulse_response = gauspuls(t0, f0, pulse.fractional_bandwidth);
impulse_response = impulse_response-mean(impulse_response); % To get rid of DC

te = (-pulse_duration/2/f0): dt : (pulse_duration/2/f0);
excitation = square(2*pi*f0*te+pi/2);
one_way_ir = conv(impulse_response,excitation);
two_way_ir = conv(one_way_ir,impulse_response);  
[~, lag] = max(abs(hilbert(two_way_ir)))

% show the pulse to check that the lag estimation is on place (and that the pulse is symmetric)
figure;
plot((1:(length(two_way_ir)))*dt -lag*dt,two_way_ir); hold on; grid on; axis tight
plot((1:(length(two_way_ir)))*dt -lag*dt,abs(hilbert(two_way_ir)),'r')
plot([0 0],[min(two_way_ir) max(two_way_ir)],'g');
legend('2-ways pulse','Envelope','Estimated lag');
title('2-ways impulse response Field II');

%% Aperture Objects
% Next, we define the the mesh geometry with the help of Field II's
% *xdc_linear_array* function.

noSubAz=round(probe.element_width/(lambda/8));        % number of subelements in the azimuth direction
noSubEl=round(probe.element_height/(lambda/8));       % number of subelements in the elevation direction
Th = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]);
Rh = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]);

% We also set the excitation, impulse response and baffle as below:
xdc_excitation (Th, excitation);
xdc_impulse (Th, impulse_response);
xdc_baffle(Th, 0);
xdc_center_focus(Th,[0 0 0]);
xdc_impulse (Rh, impulse_response);
xdc_baffle(Rh, 0);
xdc_center_focus(Rh,[0 0 0]);

%% Speckle Phantom
%
% In our next step, we define our phantom. Here, our goal is to simulate
% speckle so we have a 100 scatterrers with axial and lateral coordinates
% randomly drawn from a uniform distribution and scatterer amplitudes
% randomly drawn from a normal distribution.


xxp_speckle=random('unif',-35e-3,35e-3,number_of_scatterers,1);
zzp_speckle=random('unif',1e-3,60e-3,number_of_scatterers,1);
sca = [xxp_speckle zeros(length(xxp_speckle),1) zzp_speckle];  % list with the scatterers coordinates [m]
amp=randn(length(sca),1);                   % list with the scatterers amplitudes
cropat=round(1.1*2*sqrt((max(sca(:,1))-min(probe.x))^2+max(sca(:,3))^2)/c0/dt);   % maximum time sample, samples after this will be dumped

%%
% figure()
% % scatter3(sca(:,1), sca(:,2), sca(:,3), amp)
% plot(amp, '.')

%% Set up transmit sequence %CHANGED

no_transmits=128; 
transmit_angles = linspace(-30,30,no_transmits);
R_focus = z_focus;


%% Output data
%
% We define the variables to store our output data

t_out=0:dt:((cropat-1)*dt);                 % output time vector
STA=zeros(cropat,probe.N, no_transmits);   %CHANGED (probe.N -> no_transmits) % impulse response channel data 
%% Compute STA signals
% Now, we finally reach the stage where we generate a STA (Synthetic
    % Transmit Aperture) dataset with the help of Field II.
    
disp('Field II: Computing STA dataset');
disp('No waitbar possible for parfor, so just be patient :)');


parfor n=1:no_transmits %CHANGED probe.N -> no_transmits
    disp(['Wave no ', num2str(n)])
%     tic
%     toc
    %Since we are using parfor, we have to initate Field II and the arrays
    %for every worker as well.
%     tstart = tic;
    field_init(0);
%     tstop = toc(tstart)
%     disp(['test etter init field, init time ', num2str(tstop)])
    Th = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]);
    Rh = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]);
    
    xdc_excitation (Th, excitation);
    xdc_impulse (Th, impulse_response);
    xdc_baffle(Th, 0);
	xdc_center_focus(Th,[0 0 0]);
    xdc_impulse (Rh, impulse_response);
    xdc_baffle(Rh, 0);
    xdc_center_focus(Rh,[0 0 0]);
    
    x_focus = sind(transmit_angles(n)).*R_focus; %CHANGED
    z_focus = cosd(transmit_angles(n)).*R_focus;
    xdc_center_focus (Th, [0 0 0]);
    xdc_focus (Th, 0, [x_focus 0 z_focus]);

    xdc_focus_times(Rh, 0, zeros(1,probe.N));
    
    % do calculation
    [v,t]=calc_scat_multi(Th, Rh, sca, amp);
    
    % save data -> with parloop we need to pad the data
    if size(v,1)<cropat
        STA(:,:,n)=padarray(v,[cropat-size(v,1) 0],0,'post');    
    else
        STA(:,:,n)=v(1:cropat,:);
    end
    
    % Sequence generation
    seq(n)=uff.wave();
    seq(n).probe=probe;
    seq(n).source.xyz=[ x_focus 0  z_focus];
    seq(n).sound_speed=c0;
    seq(n).delay = -lag*dt+t; % t0 and center of pulse compensation
%     seq(n).source.azimuth = deg2rad(transmit_angles(n));
end
    
%% Channel Data
%
% In this part of the code, we creat a uff data structure to specifically
% store the captured ultrasound channel data.

channel_data = uff.channel_data();
channel_data.sampling_frequency = fs;
channel_data.sound_speed = c0;
channel_data.initial_time = 0;
channel_data.pulse = pulse;
channel_data.probe = probe;
channel_data.sequence = seq;
channel_data.data = STA./max(STA(:));
    
%% Save UFF dataset
%
% Finally, we save the data into a UFF file.

channel_data.write([data_path filesep filename],'channel_data');

%%
p = gcp('nocreate'); % If no pool, do not create new one.
delete(p)