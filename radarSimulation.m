clear all
clc;

%% Radar Specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
%% User Defined Range and Velocity of target
% define the target's initial position and velocity. 
% bNote : Velocity remains contant

% intial Range and Velocity of target
% Range cannot exceed 200m
% Velocity can be in range +/- 70m/s

tgt_range = 80;
tgt_velocity = -20;



%% FMCW Waveform Generation
% Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

%Operating carrier frequency of Radar
fc = 77e9;                          % operating frequency (Hz)
c = 3e8;                            % speed of light

d_res = 1;                          % range resolution metres
R_max = 200;                        % max range of radar 200 metres


B = c/(2*d_res);                        % Bandwidth
Tchirp = 5.5 * (2 * R_max / c);         % sweep/chrip time (time for signal to travle max radar range)
slope = B/Tchirp;                      % slope of the chirp signal


%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation.
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp.
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time.

for i=1:length(t)

    %For each time stamp update the Range of the Target for constant velocity.
    r_t(i) = tgt_range + (tgt_velocity*t(i));
    % Calculate delay time
    td(i) = 2*r_t(i)/c;

    %For each time sample we need update the transmitted and
    %received signal.
    Tx(i) = cos(2*pi*(fc*t(i)+0.5*slope*t(i)^2));
    Rx(i) = cos(2*pi*(fc*(t(i)-td(i))+0.5*slope*(t(i)-td(i))^2));

    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i);

end



%% RANGE MEASUREMENT

%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix = reshape(Mix, [Nr, Nd]);

%run the FFT on the beat signal along the range bins dimension (Nr) and
signal_fft = fft(Mix);

% Take the absolute value of FFT output & normalise
P2 = abs(signal_fft/Nr);

% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
% reject mirror image
P1 = P2(1: Nr / 2 + 1);


%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

% plot FFT output
f = Nr * (0 : (Nr / 2))/ Nr;
plot(f, P1);
xline(tgt_range-10, '-r', 'Lower Margin')
xline(tgt_range+10, '-r', 'Upper Margin')
title('1D FFT - Range of target');
xlabel('range (metres)');
ylabel('normalised FFT');
axis ([0 200 0 1]);
saveas(gcf, '1Dfft.png');

%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure('Name','2D FFT Range Doppler Response'),surf(doppler_axis,range_axis,RDM);
xlabel('velocity (m/s)');
ylabel('range (m)');
zlabel('signal strength (dB)')
title('2D FFT - Range Doppler Map');
saveas(gcf, '2Dfft.png');





%% CFAR implementation

%Slide Window through the complete Range Doppler Map

%Select the number of Training Cells in both the dimensions.
Tr = 14; % range axis
Td = 14; % doppler axis

%Select the number of Guard Cells in both dimensions around the Cell under
%test (CUT) for accurate estimation
Gr = 4; % range axis
Gd = 4; % doppler axis

% offset the threshold by SNR value in dB
% Define Offset (Adding room above noise threshold for the desired SNR)
offset = 20;

%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);

%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.


% Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
% CFAR

grid_size = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1);
guard_region = (2*Gr+1)*(2*Gd+1);
training_cells = grid_size - guard_region;

% Initialize Vector to hold final signal after thresholding
signal_cfar = zeros(size(RDM));

% Slide across the complete 2D matrix
% loop over the range cells, protecting training and guard
for i = Tr+Gr+1:Nr/2-(Gr+Tr)
    % loop over the doppler cells
    for j = Td+Gd+1:Nd-(Gd+Td)
        
        %select the matrix which includes the training, guard & CUT
        % x = left of CUT to right of CUT
        % y = top of CUT to bottom of CUT
        grid = RDM(i-Gr-Tr : i+Gr+Tr, j-Gd-Td : j+Gd+Td);
        % specify the guard area and set levels to 0
        grid(i-Gr : i+Gr, j-Gd : j+Gd) = 0;
        
        % convert training cells from log to linear with db2pow
        power = db2pow(grid);
        
        % measure average noise level across all training cells
        noise_level = sum(power) / training_cells;
        
        % convert back to log using pow2db
        db = pow2db(noise_level);
        
        %X = [' Noise: ', num2str(db)];
        %disp(X)
        
        % add the noise threshold offset
        threshold = db + offset;
        
        % compare the CUT against threshold
        if RDM(i,j) > threshold
            signal_cfar(i, j) = 1;
            Y = ['RDM value :',num2str(RDM(i,j))];
            %disp(Y)
             
        end
    end
end

% The process above will generate a thresholded block, which is smaller
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0.

%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure('Name','CA-CFAR Filtered Range Doppler Map'),surf(doppler_axis,range_axis,signal_cfar);
xlabel('velocity (m/s)');
ylabel('range (m)');
zlabel('signal strength (dB) Normalised')
title('CA-CFAR Filtered Range Doppler Map');
colorbar;
saveas(gcf, 'CFAR.png');
