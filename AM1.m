
array=[-2 0 0; -1 0 0; 0 0 0; 1 0 0; 2 0 0];

directions=[30, 0 ; 35, 0 ; 90, 0];% azimuth and elevation


%----1. Form the patten of the above array---------------------------------

Z1=pattern(array);

%plot2d3d(Z1,[0:180],0,'gain in dB','initial pattern');

%----2. Theoretical Covariance Matrix Formation-----------------------
S=spv(array, directions);
Rmm=eye(3);  
sigma2=0.0001;
Rxx_theoretical=S*Rmm*S'+sigma2*eye(5,5);

%----3. Practical Covariance Matrix-------------------------------------
load Xaudio.mat;
load Ximage.mat;

%sound(abs(X_au(2,:)), 11025);
displayimage(X_im(2,:),image_size, 201,'The received signal at the 2nd antenna');
Rxx_au = X_au * X_au'/length(X_au(1,:));
Rxx_im = X_im * X_im'/length(X_im(1,:));

%----4.Forget that you know that there are 3 sources present----------------
directions=[];
Rmm=[];
S=[];
sigma2=[];

%----5. Detection Problem--------------------------------------------------
input=eig(Rxx_im);

% %----6. Estimation Problem - based on the 'array pattern'-----------------
% Sd=spv(array,[90,0]);
% a=1;
% wopt= a*inv(Rxx_theoretical) *Sd;
% Z=pattern(array, wopt);
% figure,plot2d3d(Z, [0:180], 0, 'gain in dB','W-H array pattern');
% 
% %----7.Repeat instructions 2, 4, 5.1 and 6 but with noise level 10dB below the level of the sources
directions=[30, 0 ; 35, 0 ; 90, 0];
S=spv(array, directions);
Rmm=eye(3);  
sigma2=0.1;
Rxx_theoretical=S*Rmm*S'+sigma2*eye(5,5);
% 
input2=eig(Rxx_theoretical);
Sd=spv(array,[90,0]);
a=1;
wopt= a*inv(Rxx_theoretical) *Sd;
Z=pattern(array, wopt);
figure,plot2d3d(Z, [0:180], 0, 'gain in dB','W-H array pattern');
% 
% %----8. What conclusions can be drawn from 6 and 7 -----------------------
% %we can easily draw a conclusion that when the SNR is quite small, the
% %Wiener Hopf will invalid
% 
% %----9. Estimation Problem------------------------------------------------
% Z = music(array, Rxx_theoretical);
% plot2d3d(Z,[0:180],0,'dB', 'MuSIC spectrum');
% 
% %----10. Repeat instruction 9 by using and Rxx_au Rxx_im
Z1 = music(array, Rxx_au);
figure,plot2d3d(Z1,[0:180],0,'dB', 'MuSIC spectrum');
title('array pattern for Rxxau');

% Z2 = music(array, Rxx_im);
% figure,plot2d3d(Z2,[0:180],0,'dB', 'MuSIC spectrum');
% title('array pattern for Rxxim');
% 
% %----11. Multipaths - Coherent Sources:------------------------------------
% 
% directions=[30, 0 ; 35, 0 ; 90, 0];
% S=spv(array, directions);
% Rmm=[1,1,0;
%      1,1,0;
%      0,0,1;
%     ];
% sigma2=0.0001;
% Rxx_theoretical=S*Rmm*S'+sigma2*eye(5,5);
% Z2 = music(array, Rxx_theoretical);
% plot2d3d(Z2,[0:180],0,'dB', 'MuSIC spectrum');
% 
% 
Rxx=issp(Rxx_theoretical,4);
array1=[-2 0 0; -1 0 0; 0 0 0; 1 0 0];
Z2 = music(array1, Rxx);
plot2d3d(Z2,[0:180],0,'dB', 'MuSIC spectrum');
% 
% 
% %----12. Reception Problem-------------------------------------------------
Sd=spv(array,[90,0]);
wopt=inv(Rxx_au) *Sd;
yt=wopt'*X_au;

%soundsc(real(yt), 11025);

wopt=inv(Rxx_im) *Sd;
yt=wopt'* X_im;
yt=Normalize(yt);
displayimage(yt, image_size, 202,'The received signal at o/p of W-H beamformer');

Z2 = music(array, Rxx_im);
plot2d3d(Z2,[0:180],0,'dB', 'MuSIC spectrum');
% 
% %%%%%%12.3%%%%$----------------------------------------------------------------
% 
% Sdesired=spv(array, [90,0]);
% S_j=spv(array,[30,0;35,0]);
% P_a=eye(5)- (S_j*inv(S_j'*S_j)*S_j');
% w=P_a*Sdesired;
% y=w'*X_im;
% y=Normalize(y);
% displayimage(y, image_size, 202,'The received signal at o/p of W-H beamformer');
% 
% Sdesired=spv(array, [30,0]);
% S_j=spv(array,[35,0;90,0]);
% P_a=eye(5)- (S_j*inv(S_j'*S_j)*S_j');
% w=P_a*Sdesired;
% y=w'*X_im;
% y=Normalize(y);
% displayimage(y, image_size, 202,'The received signal at o/p of W-H beamformer');
% 
% Sdesired=spv(array, [35,0]);
% S_j=spv(array,[30,0;90,0]);
% P_a=eye(5)- (S_j*inv(S_j'*S_j)*S_j');
% w=P_a*Sdesired;
% y=w'*X_im;
% y=Normalize(y);
% displayimage(y, image_size, 202,'The received signal at o/p of W-H beamformer');
% 
% %----13. Practical Detection Criteria:-------------------------------------
% %%%creat new resources with 250 snapshots%%%
% n = 250;
% N_array= 5;   %N
% N_sources= 3; %M
% snr=10;
% Input=(rand(N_sources,n)+j*rand(N_sources,n))/sqrt(2);
% noise=(rand(N_sources,n)+j*rand(N_sources,n))/sqrt(2)/snr;
% array=[-2 0 0; -1 0 0; 0 0 0; 1 0 0; 2 0 0];
% directions=[30, 0 ; 35, 0 ; 90, 0];
% 
% 
% S_i=spv(array,directions);
% X=S_i*Input;
% X1=awgn(X,snr,'measured');
% Rxx=X1*X1'/n;
% 
% L=250;
% N=5;
% d= eig(X1*X1'/L);             %calculate the eig_value
% d=sort(d,'ascend');         %reorder the data from small to lagre
% d_N=1;
% d_n=[];
% 
% for i=1:N
%     d_N=d_N*d(i);
%     d_n=[d_n,d_N];
% end
% d_n=fliplr(d_n).'; 
% 
% d_N=0;
% d_n1=[];
% for i=1:N
%     d_N=d_N+d(i);
%     d_n1=[d_n1,d_N];
% end
% d_n1=fliplr(d_n1).';    %d_n1 is just the right order of matrix(multi part)
% 
% AIC=[];
% k1=(N:-1:1).';
% k2=(0:(N-1)).';
% k3=(2*N:-1:(N+1)).';
% AIC=-2*L.*(log(d_n)+k1.*(log(k1)-log(d_n1)))+2.*k2.*k3;
% MDL=-L.*(log(d_n)+k1.*(log(k1)-log(d_n1)))+0.5.*k2.*k3*log(L);
% 
