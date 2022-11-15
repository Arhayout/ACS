clear all
close all
clearvars
clc
% Time of arrival of each path [s]
tau =[0 0.2 0.5 1.6 2.3 5 ] * 1e-6;
% Average gain of each path [SI]
sigma2a = 10.^([-3 0 -2 -6 -8 -10]./ 10) ;
% Total number of paths
K=length ( tau ) ;
% Bandwidth of the transmitted signal [Hz]
B=19e6 ;
% Sampling period [s]
Ts=1/B ;
% Transmitted block length (should be greater than N1+N2+1)
N=1000;


% Maximum Doppler frequency [Hzfd =100e4 ;
% Doppler filter sampling frequency and period
fs=B;
% Doppler filter sampling period
ts =1/ fs ;
% Doppler filter length
M=1024;
fd =100e4 ;% Number of taps (N1+N2+1)
N1=10;
N2=ceil(tau(end)/Ts)+ 10;

%White noise generation
 Np= M+N;
 var2=1/2;
 x=sqrt(var2/2)*randn(K,Np)+j*randn(K,Np);
 var(x,0,"all");

%% 
%Doppler filter h(n)
m=(0:M-1);
h = gamma(3/4) .*( fd ./ abs( pi .*( m - M/2 ).*ts )).^ (1/4) .* besselj(1/4,fd .* abs( 2 * pi .*( m - M/2).*ts ));
h(M/2+1) = gamma(3/4)/(gamma(5/4)*fd.^0.5);
hnorm=h/norm(h);
var(h,1,"all")
at=filter(h, 1, x, [],2);
a=at(:,M+1:end);
size(a)
z=(1./var(a.')).' .*a;
size(z);
var(a',1);
var(z',1);

e=a(1,:);
N=length(e);
xdft = fft(e);
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/length(e):fs/2;
plot(freq,pow2db(psdx))
grid on
title("Periodogram Using FFT")
xlabel("Frequency (Hz)")
ylabel("Power/Frequency (dB/Hz)")

figure
pwelch(e,[],[],[],1,'centered')
figure
plot(e)
%%
%%1.3.5 : Using the channel with a QPSK signal
msg=2*round(rand(1,M))-1;
u = msg (1:2:end);
v = msg (2:2:end);
d = u + j*v;
scatterplot(d)
U= toeplitz(d,zeros(N1+N2+1,1)).';
y= sum(gt.*U,1);
scatterplot(y(11:end))

