%% Computation of a FI dataset with Field II and beamforming with USTB
%
% This example shows how to load the data from a Field II simulation into 
% USTB objects, and then beamformt it with the USTB routines. 
% This example uses the P4-2v 64 element Verasonics Transducer
% The Field II simulation program (field-ii.dk) should be in MATLAB's path.
%
% This example is imaging with focused transmit waves (Focused Imaging-FI).
% The example also demonstates calculation of the coherence factor and some
% functionality to plot the images using built in USTB routines, MATLAB
% commands and some details on scan conversion.
% 
%
% authors:  Ole Marius Hoel Rindal <olemarius@olemarius.net>
%           Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
%
% Last updated: 15.01.2020

% Added "doSparse" flag. If sat, a sparse Tx and Rx-design 

function channel_data = fieldii_generate_dataset(path, filename, pointscatters, h)
%     h.dt = 0;
    %% basic constants
    
    if ~isfield(h, 'c0'); h.c0=1540; end    % Speed of sound [m/s]
    if ~isfield(h, 'fs'); h.fs=100e6; end   % Sampling frequency [Hz]
    if ~isfield(h, 'f0'); h.f0=2.56e6; end                              % Transducer center frequency [Hz]
    if ~isfield(h, 'bw'); h.bw=0.67; end                               % probe bandwidth [1]
    if ~isfield(h, 'elm_height'); h.elm_height = 5e-3; end
    if ~isfield(h, 'probepitch'); h.probepitch = 0.300e-3; end
    if ~isfield(h, 'kerf'); h.kerf=0.050e-3; end                         % gap between elements [m]
    if ~isfield(h, 'lens'); h.lens_el=60e-3; end                          % position of the elevation focus
    if ~isfield(h, 'N_elements'); h.N_elements = 64; end
    if ~isfield(h, 'pulse_duration'); h.pulse_duration=2.5; end                     % pulse duration [cycles]
    if ~isfield(h, 'z_focus'); h.z_focus =40/1000; end                       % Transmit focus
    if ~isfield(h, 'frac_bw');  h.frac_bw=0.65; end                       % Transmit focus
    if ~isfield(h, 'N_transmits'); h.N_transmits =128; end                       % Transmit focus
    if ~isfield(h, 'angle'); h.angle = 30; end                       % Transmit focus
    if ~isfield(h, 'R_focus'); h.R_focus=40/1000; end                       % Transmit focus
    if ~isfield(h, 'degs_separated'); h.degs_separated=-3; end                       % Transmit focus
    if ~isfield(h, 'probe_type'); h.probe_type="linear_array"; end                       % Transmit focus
    if ~isfield(h, 'D'); h.D=8; end                       % Transmit focus
%     if ~isfield(h, 'save'); h.save=1; end                       % Transmit focus


    h.dt=1/h.fs;     % Sampling step [s] 
%     doSparse = 0;
    %% field II initialisation
    field_init(0);
    set_field('c',h.c0);              % Speed of sound [m/s]
    set_field('fs',h.fs);             % Sampling frequency [Hz]
    set_field('use_rectangles',1);  % use rectangular elements
    
    %% transducer definition P4-2v Verasonics 64-element phased
    % 
    % Our next step is to define the ultrasound transducer array we are using.
    % For this experiment, we shall use the L11-4v 128 element Verasonics
    % Transducer and set our parameters to match it.
    if h.probe_type == "linear_array"
        probe = uff.linear_array();
    end
    lambda=h.c0/h.f0;                           % Wavelength [m]
    probe.element_height=h.elm_height;              % Height of element [m]
    
    probe.pitch =h.probepitch;                  % pitch [m]
    probe.element_width=probe.pitch-h.kerf;   % Width of element [m]
    probe.N=h.N_elements;                             % Number of elements
    
    %% pulse definition
    pulse = uff.pulse();

    pulse.center_frequency = h.f0;
    pulse.fractional_bandwidth = h.frac_bw;        % probe bandwidth [1]
    t0 = (-1/pulse.fractional_bandwidth/h.f0): h.dt : (1/pulse.fractional_bandwidth/h.f0);
    impulse_response = gauspuls(t0, h.f0, pulse.fractional_bandwidth);
    impulse_response = impulse_response-mean(impulse_response); % To get rid of DC
    
    te = (-h.pulse_duration/2/h.f0): h.dt : (h.pulse_duration/2/h.f0);
    excitation = square(2*pi*h.f0*te+pi/2);
    one_way_ir = conv(impulse_response,excitation);
    two_way_ir = conv(one_way_ir,impulse_response);  
    [~, lag] = max(abs(hilbert(two_way_ir)))
    
    %% aperture objects
    % definition of the mesh geometry
    noSubAz=round(probe.element_width/(lambda/h.D));        % number of subelements in the azimuth direction
    noSubEl=round(probe.element_height/(lambda/h.D));       % number of subelements in the elevation direction
    Th = xdc_linear_array (probe.N, probe.element_width, probe.element_height, h.kerf, noSubAz, noSubEl, [0 0 Inf]); 
    Rh = xdc_linear_array (probe.N, probe.element_width, probe.element_height, h.kerf, noSubAz, noSubEl, [0 0 Inf]); 
    
    % setting excitation, impulse response and baffle
    xdc_excitation (Th, excitation);
    xdc_impulse (Th, impulse_response);
    xdc_baffle(Th, 0);
    xdc_center_focus(Th,[0 0 0]);
    xdc_impulse (Rh, impulse_response);
    xdc_baffle(Rh, 0);
    xdc_center_focus(Rh,[0 0 0]);
    
    %% Set up the transmit sequence
%     h.N_transmits = 4;
    no_transmits=h.N_transmits; 
    transmit_angles = linspace(-h.angle,h.angle,no_transmits);
%     h.R_focus = 60/1000;
    
    
    %% Create a phantom to image
    degs_separated = h.degs_separated;
    n_points = size(pointscatters,1);
    if h.amp == "rand"
        phantom_amplitudes(1:n_points) = rand(n_points,1);
    elseif h.amp == "ones"
        phantom_amplitudes(1:n_points) = ones(n_points,1);
    end
    %% output data
    cropat=round(1.1*2*sqrt((max(pointscatters(:,1))-min(probe.x))^2+max(pointscatters(:,3))^2)/h.c0/h.dt);   % maximum time sample, samples after this will be dumped
    data=zeros(cropat,probe.N,no_transmits);    % impulse response channel data
    
    %% Compute STA signals
    fprintf('Field II: Computing FI dataset \n \n');
    disp('~')
    for n=1:no_transmits
        s = sprintf('\nSimulating transmit %d / %d',n,no_transmits);
        b = repmat('\b', [1, length(s)]);
        fprintf(1, [b, s]);
    
        x_focus = sind(transmit_angles(n)).*h.R_focus;
        z_focus = cosd(transmit_angles(n)).*h.R_focus;
    
	    % Set the focus for this direction with the proper reference point
        xdc_center_focus (Th, [0 0 0]);
        xdc_focus (Th, 0, [x_focus 0 z_focus]);
%         if doSparse
%             xdc_apodization (Th, 0, sparseArray);
%         else
            xdc_apodization (Th, 0, ones(1,probe.N));
%         end
        
        % receive aperture    
%         if doSparse
%             xdc_apodization (Rh, 0, sparseArray);
%         else
            xdc_apodization (Rh, 0, ones(1,probe.N));
%         end
    
        xdc_focus_times(Rh, 0, zeros(1,probe.N));
    
        % do calculation
        [v,t]=calc_scat_multi(Th, Rh, pointscatters, phantom_amplitudes');
    
        % save data -> with parloop we need to pad the data
        if size(v,1)<cropat
            data(:,:,n)=padarray(v,[cropat-size(v,1) 0],0,'post');    
        else
            data(:,:,n)=v(1:cropat,:);
        end
        
        %% SEQUENCE GENERATION
        seq(n)=uff.wave();
        seq(n).probe=probe;
        seq(n).source.xyz=[x_focus 0 z_focus];
        seq(n).sound_speed=h.c0;
        seq(n).delay = -lag*h.dt+t;
    end
    
    
    %% CHANNEL DATA
    channel_data = uff.channel_data();
    channel_data.sampling_frequency = h.fs;
    channel_data.sound_speed = h.c0;
    channel_data.initial_time = 0;
    channel_data.pulse = pulse;
    channel_data.probe = probe;
    channel_data.sequence = seq;
    channel_data.data = data./max(data(:)) + 1000*eps*randn(size(data));
    clear data
    %%
    if ~isfield(h, 'save')
        answer = questdlg('Would you like to save the dataset?', ...
	        'Save Dataset', ...
	        'Yes','No','Yes');
        % Handle response
        switch answer
            case 'Yes'
                channel_data.write(strcat(path,'\', filename, '.uff'),'channel_data')
%                 disp([answer ' coming right up.'])
%                 dessert = 1;
            case 'No'
                disp('Did not save dataset')
        end
    else
        if h.save
            channel_data.write(strcat(path,'\', filename, '.uff'),'channel_data')
        end
    end
end