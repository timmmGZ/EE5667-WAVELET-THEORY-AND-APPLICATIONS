function phi(wname)
% phi_daub.m

% Generates successive approximations of the db2 wavelet,
% uses 4-tap low-pass Daubechies synthesis filter 
% user declares required number of iterations and decides 
% (specifying number of ite  rations) which four approximations 
% are to be grphically illustrated.
numits =10;
plot1 = 1;
plot2 =3;
plot3 = 7;
plot4 = 10;
[LoD,HiD,LoR,HiR] = wfilters(wname); 
LoR = qmf(HiR,1) %or  LoR = (-1).^(1:length(HiR)) .*wrev(HiR);
% reverse High-pass,and then odd index of coefficients * -1

L(1) = length(LoR);
g = zeros(numits,((2^numits-1)*(L(1)-1)+1));
for i = 1:length(LoR) 
    g(1,i) = LoR(i);
end
phi = zeros(numits,((2^numits-1)*(L(1)-1)+1));
phi(1,:) = g(1,:);

for i = 2:numits
   
   L(i) = (2^i -1)*(L(1) - 1) + 1;
   
   for n = 1:L(i)
      
      k = 0;
      while (n - 2^(i-1)*k) > 0
         tempindex = n - 2^(i-1)*k;
         g(i,n) = g(i,n) + g(1,k+1)*g(i-1,tempindex);
         k = k + 1;
      end     
      phi(i,n) = (2^((i-1)/2))*g(i,n);         
   end
        
end

dt1 = 1/2^plot1;
time1 = -dt1:dt1:L(plot1)*dt1;
xmin1 = -dt1;
xmax1 = L(plot1)*dt1+2*dt1;
plotphi1 = [0 phi(plot1,1:L(plot1)) 0];
ymin1 = min(plotphi1) - 0.1;
ymax1 = max(plotphi1) + 0.1;

dt2 = 1/2^plot2;
time2 = -dt2:dt2:L(plot2)*dt2;
xmin2 = -dt2;
xmax2 = L(plot2)*dt2+2*dt2;
plotphi2 = [0 phi(plot2,1:L(plot2)) 0];
ymin2 = min(plotphi2) - 0.1;
ymax2 = max(plotphi2) + 0.1;

dt3 = 1/2^plot3;
time3 = -dt3:dt3:L(plot3)*dt3;
xmin3 = -dt3;
xmax3 = L(plot3)*dt3+2*dt3;
plotphi3 = [0 phi(plot3,1:L(plot3)) 0];
ymin3 = min(plotphi3) - 0.1;
ymax3 = max(plotphi3) + 0.1;

dt4 = 1/2^plot4;
time4 = -dt4:dt4:L(plot4)*dt4;
xmin4 = -dt4;
xmax4 = L(plot4)*dt4+2*dt4;
plotphi4 = [0 phi(plot4,1:L(plot4)) 0];
ymin4 = min(plotphi4) - 0.1;
ymax4 = max(plotphi4) + 0.1;

figure(1)
clf

subplot(4,2,1)
plot([xmin1 xmax1], [0 0],'k-')
hold on
stairs(time1,plotphi1)
axis([xmin1 xmax1 ymin1 ymax1]); 
xlabel('Time')
ylabel(strcat('Amplitude Phi(',int2str(plot1),')'))

subplot(4,2,3)
plot([xmin2 xmax2], [0 0],'k-')
hold on
stairs(time2,plotphi2)
axis([xmin2 xmax2 ymin2 ymax2]); 
xlabel('Time')
ylabel(strcat('Amplitude Phi(',int2str(plot2),')'))

subplot(4,2,5)
plot([xmin3 xmax3], [0 0],'k-')
hold on
stairs(time3,plotphi3)
axis([xmin3 xmax3 ymin3 ymax3]); 
xlabel('Time')
ylabel(strcat('Amplitude Phi(',int2str(plot3),')'))

subplot(4,2,7)
plot([xmin4 xmax4], [0 0],'k-')
hold on
stairs(time4,plotphi4)
axis([xmin4 xmax4 ymin4 ymax4]); 
xlabel('Time')
ylabel(strcat('Amplitude Phi(',int2str(plot4),')'))
end
