Vp% Probability distribution from histogram
% date:1 May 2020
%===================================================

% load data
datax=load('KLSI_return.txt'); nx=length(datax) % data length
figure(1)
subplot(3,1,1); 
plot(datax); xlabel('time, t'); ylabel('return,r(t)');

% histogram
nbins=round(nx*0.05)    % size of bin is 5 % of data length
[Nr,x]=hist(datax,nbins);
subplot(3,1,2); 
plot(x,Nr); xlabel('return, r'); ylabel('Histogram,N(r)');

% Probability
TotalN=sum(Nr);
prob=Nr./TotalN;
subplot(3,1,3); 
plot(x,prob); xlabel('return, r'); ylabel('Probability, P(r)');
