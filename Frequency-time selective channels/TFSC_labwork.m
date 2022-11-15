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
% Number of taps (N1+N2+1)
N1=10;
N2=ceil (tau (end) /Ts) + 10;
%%
x = sqrt(sigma2a/2).'.*randn(K,N);
y = sqrt(sigma2a/2).'.*randn(K,N);
a = x + j * y;
for k=1:K
  mean(a(k,:))
end
%%
alpha = sinc(repmat(tau,N2+N1+1,1)/Ts-repmat([-N1:1:N2].', 1 ,K));
at_c=[1 0 0 0 0 0];
at= repmat(at_c.',1,N);
gt=alpha * at;
size(alpha)
size(gt)
figure
plot(abs(at(:,1:10)))
hold on
t= -N1*Ts:Ts:N2*Ts;
plot(t,abs(gt(:,1:20)))
hold off
%%
freq_res=fft(gt(:,k).',8*size(gt(:,k),1));
psd=abs(freq_res).^2/(size(freq_res,2));
freq= linspace(-B/2,B/2,size(freq_res,2));
figure
plot(freq,pow2db(psd))
%%
msg=2*round(rand(1,2*N))-1;
u = msg (1:2:end);
v = msg (2:2:end);
st = u + j*v;
scatterplot(st)
U= toeplitz(st,zeros(N1+N2+1,1)).';
y= sum(gt.*U,1);
size(y)
size(filter(gt(:,1),1,st))
norm(y- filter(gt(:,1),1,st))
scatterplot(y(11:end))



