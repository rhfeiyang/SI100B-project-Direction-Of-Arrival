clear all; close all;
%% Load data

load("Observation_nb.mat");                      % load data
X = y ;                                          % 4-channel received signals
fs = fs;                                         % sample rate (Hz)
L=length(X(:,1));
%% Plot waveform

subplot(2,1,1);
ts=(0:L-1)/fs;
plot(ts,real(X))
title("time domain")
xlabel('time(s)')
legend('signal1','signal2','signal3','signal4')

subplot(2,1,2);
f=(-fs/2:fs/L:(-fs/2+fs/L*(L-1)));
plot(f,abs(fftshift(fft(X))))
title("frequency domain")
xlabel('frequency(Hz)')
legend('signal1','signal2','signal3','signal4')
%% Array setup

[Frame,nSensors] = size(X);
J = nSensors;                                      % number of sensors
dx = 0.034;                                        % inter-sensor distance in x direction (m)
dy = 0;                                            % sensor distance in y direction (m)
c = 340;                                           % sound velocity  (m)
n_source = 2;                                      % number of sources
Index = linspace(0,J-1,J);
p = (-(J-1)/2 + Index.') * [dx dy];                % sensor position
%% Plot sensor positions

linspec = {'rx','MarkerSize',12,'LineWidth',2};
figure
plot(p(:,1),p(:,2),linspec{:});
title('Sensor positions');
xlabel('x position in meters');
ylabel('y position in meters');
disp('The four microphones are ready !');
%% DoA estimation (MUSIC)

stride =0.5;                                                 % determine the angular resolution(deg)
theta = -90:stride:90;                                       % grid
f_c =2300 ;                                                  % By inspection, we get the center frequency(Hz) from from the ferquency-domain plot
R_x = X'*X/L;                                                % autocorrelation estimate
v =  [sin(theta*pi/180);-cos(theta*pi/180)];                 % direction vector
a_theta = exp(-2j*pi*(Index')*(dx*sin(theta*pi/180)) / (c/f_c));                  % steer vector

% implement eigen-decomposition and obtain the pseudo spectrum
[EV,D] = eig(R_x);
EVA = diag(D);
[EVA,ind] = sort(EVA);
EV = EV(:,ind);
Un = EV(:,1:J-n_source);                                        % noise subspace (columns are eigenvectors), size: J*(J-n_source)
P_sm = 1./diag(a_theta'*(Un*Un')*a_theta);                    % pseudo music power
%% Plot the MUSIC pseudo power spectrum

figure;
linspec = {'b-','LineWidth',2};
plot(theta, 10*log10(abs(P_sm)), linspec{:});
title('MUSIC pseudo power spectrum')
xlabel('Angle in [degrees]');
ylabel('Power spectrum in [dB]');
xlim([-90,90]);
%% Find the local maximum and visualization

P_middle = abs(P_sm(2:end-1));
P_front = abs(P_sm(1:end-2));
P_back = abs(P_sm(3:end));
logic_front = (P_middle - P_front)>0;
logic_back = (P_middle - P_back)>0;
logic = logic_front & logic_back;
P_middle(~logic) = min(P_middle);
P_local = [abs(P_sm(1));P_middle;abs(P_sm(end))];
[~,doa_Idx] = maxk(P_local,n_source);
doa = theta(doa_Idx);
[~,minIdx] = min(abs(doa));
doa_source = doa(minIdx);
[~,maxIdx] = max(abs(doa));
interfer = doa(maxIdx);
disp(['The desired source DOA with MUSIC is: ',num2str(doa_source),' deg']);
disp(['The interfering DOA with MUSIC is: ',num2str(interfer),' deg']);