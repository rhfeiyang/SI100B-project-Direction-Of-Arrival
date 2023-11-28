clear all; close all;
%% Load data

load("Observation_wb.mat");                                         % load data

X = X ;                                          % 4-channel received signals
fs = fs;                                          % sample rate (Hz)
L=length(X(:,1));
dx=0.025;
J = 4;                                             % number of sensors                                        % inter-sensor distance in x direction (m)
dy = 0;                                            % sensor distance in y direction (m)
c = 340;                                           % sound velocity  (m)
n_source = 2;                                      % number of sources
Index = linspace(0,J-1,J);
p = (-(J-1)/2 + Index.') * [dx dy];                % sensor position

soundsc(real(X(:,1)),fs);
pause(3)
%% Plot waveform

subplot(2,1,1);
ts=(0:L-1)/fs;
plot(ts,real(X(:,1)))
title("time domain")
xlabel('time(s)')
legend('signal1')

subplot(2,1,2);
f=(-fs/2:fs/L:(-fs/2+fs/L*(L-1)));
plot(f,abs(fftshift(fft(X(:,1)))))
title("frequency domain")
xlabel('frequency(Hz)')
legend('signal1')
%% 
% STFT

window = 512;
noverlap = window/2;
nfft = 512;
stride=0.5;

[ss,fc]=doastft(X,window,noverlap,nfft,fs);
Pmusic=zeros(361,1);
for x=(1:window/2+1)
    Xr=squeeze(ss(x,:,:));
    f_c=fc(x);
    Pmusic=Pmusic+music(dx,Index,f_c,Xr,n_source,stride);
end
figure;
theta = -90:stride:90;
P_sm=nfft./Pmusic;
linspec = {'b-','LineWidth',2};
plot(theta, 10*log10(abs(P_sm)), linspec{:});
title('MUSIC pseudo power spectrum')
xlabel('Angle in [degrees]');
ylabel('Power spectrum in [dB]');
xlim([-90,90]);
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

for x=1:length(doa)
    disp(['The desired source',num2str(x),' DOA with MUSIC is: ',num2str(doa(x)),' deg']);
end
%% 
M_1 = audioread('Array_output_01.wav');
M_2 = audioread('Array_output_02.wav');
M_3 = audioread('Array_output_03.wav');
M_4 = audioread('Array_output_04.wav');
M = zeros(length(M_1),4);
M(:,1) = M_1;
M(:,2) = M_2;
M(:,3) = M_3;
M(:,4) = M_4;
LM = length(M(:,1));
dx=0.025;
J = 4;                                             % number of sensors                                        % inter-sensor distance in x direction (m)
dy = 0;                                            % sensor distance in y direction (m)
c = 340;                                           % sound velocity  (m)
n_source = 2;                                      % number of sources
Index = linspace(0,J-1,J);
p = (-(J-1)/2 + Index.') * [dx dy];                % sensor position
soundsc(real(M(:,1)),fs);
pause(10)
figure()
subplot(2,1,1);
ts=(0:LM-1)/fs;
plot(ts,real(M(:,1)))
title("time domain")
xlabel('time(s)')
legend('signal1')

subplot(2,1,2);
f=(-fs/2:fs/LM:(-fs/2+fs/LM*(LM-1)));
plot(f,abs(fftshift(fft(M(:,1)))))
title("frequency domain")
xlabel('frequency(Hz)')
legend('signal1')

window = 512;
noverlap = window/2;
nfft = 512;
stride=0.5;

[ss_1,fc_1]=doastft(M,window,noverlap,nfft,fs);
Pmusic=zeros(361,1);
for x=(1:window/2+1)
    Xr=squeeze(ss_1(x,:,:));
    f_c=fc_1(x);
    Pmusic=Pmusic+music(dx,Index,f_c,Xr,n_source,stride);
end
figure;
theta = -90:stride:90;
P_sm=nfft./Pmusic;
linspec = {'b-','LineWidth',2};
plot(theta, 10*log10(abs(P_sm)), linspec{:});
title('MUSIC pseudo power spectrum')
xlabel('Angle in [degrees]');
ylabel('Power spectrum in [dB]');
xlim([-90,90]);
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

for x=1:length(doa)
    disp(['The desired source',num2str(x),' DOA with MUSIC is: ',num2str(doa(x)),' deg']);
end




function [ss,fc]=doastft(X,window,overlap,nfft,fs)
    xlen=length(X);
    rown = ceil((1+nfft)/2);
    coln = 1+fix((xlen-window)/(window-overlap));
    ss = zeros(rown, coln,4);
    co=0;
    for i=1:(window-overlap):xlen-window+1
        xtw=X(i:i+window-1,:).*hamming(window);
        s=fft(xtw,nfft);
        co=co+1;
        ss(:,co,:)=s(1:rown,:);
    end
    fc=(0:rown-1)*fs/nfft;
end
function Pmusic=music(dx,Index,f_c,Xr,n_source,stride)
    J=4;
    L=size(Xr,1);
    c=340;
    theta = -90:stride:90;                                       % grid
    R_x = Xr'*Xr/L;                                                % autocorrelation estimate
    a_theta = exp(-2j*pi*(Index')*(dx*sin(theta*pi/180)) / (c/f_c));                  % steer vector
    % implement eigen-decomposition and obtain the pseudo spectrum
    [EV,D] = eig(R_x);
    EVA = diag(D);
    [~,ind] = sort(EVA);
    EV = EV(:,ind);
    Un = EV(:,1:J-n_source);                                        % noise subspace (columns are eigenvectors), size: J*(J-n_source)
    Pmusic = diag(a_theta'*(Un*Un')*a_theta); 
end