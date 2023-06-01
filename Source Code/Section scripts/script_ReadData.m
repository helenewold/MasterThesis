% Sjekker om filen inneholder ferdig BeamFormet data 
display=true;
content = uff.index(filename,'/',display);

% Leser channel data hvis ikke (?)
channel_data=uff.read_object(filename,'/channel_data');
channel_data.data = channel_data.data(:,:,:,1);


% Definerer scan
azimuth_axis=zeros(channel_data.N_waves,1);

depth_axis = linspace(0e-3, 58e-3, 512).'; %    z_axis=linspace(1e-3,55e-3,512).';


for n=1:channel_data.N_waves
    azimuth_axis(n)=channel_data.sequence(n).source.azimuth;
end

channel_data.N_frames = 1;
