%% preprocess
clear, clc, clf, close all;

% load regular_plastic.mat;
% load composite.mat;
% load carrier_d1.mat;
load carrier_d2.mat;

data = Fri_coeff;

L = length(Time);
z2=wextend(1,'sym',data,round(length(data)/2));
wlen=10;
hop=1;
z2=wkeep1(z2,L+1*wlen);

h=hamming(wlen);
f=1:0.5:50;
fs = .1;
[tfr2,f,t2]=spectrogram(z2,h,wlen-hop,f,fs);
tfr2=tfr2*2/wlen*2;
figure(4)
imagesc(t2+0-wlen/fs/2,f,abs(tfr2))
