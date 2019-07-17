function psi( wname )
numits =10;
plot1 = 1;
plot2 =3;
plot3 = 7;
plot4 = 10;

[LoD,HiD,LoR,HiR] = wfilters(wname); 
HiR =-( (-1).^(1:length(LoR)) ).*wrev(LoR);
L(1)=length(HiR) 

g = zeros(numits,((2^numits-1)*(L(1)-1)+1));
F = zeros(1,L(1));

for i = 1:length(HiR) 
    g(1,i) = HiR(i);
    F(i) = LoR(i);
end
psi = zeros(numits,((2^numits-1)*(L(1)-1)+1));
psi(1,:) = g(1,:);

for i = 2:numits
    L(i) = (2^i -1)*(L(1) - 1) + 1;
   for n = 1:L(i) 
      k = 0;
      while (k <L(1))%filter size is L(1)
         tempindex = n-k;
         if(tempindex<=0)%without padding, in beginning we cant convolve
             break;
         end
         if(mod(tempindex,2)== 1)
             a=g(i-1,floor(tempindex/2)+1); 
         elseif(mod(tempindex,2)== 0)
             a=0;  %since upsample fill zeros, now the corresponding smaple of filter meet the 0        
         end
        %disp(tempindex+" convolutioning"+ F(k+1)+'*'+a);
         g(i,n) = g(i,n) + F(k+1)*a;
         k = k + 1;
      end   
      %disp('==========finish a single convolution=========')
      psi(i,n) = (2^((i-1)/2))*g(i,n);
   end        
end

dt1 = 1/2^plot1;
time1 = -dt1:dt1:L(plot1)*dt1;
xmin1 = -dt1;
xmax1 = L(plot1)*dt1+2*dt1;
plotpsi1 = [0 psi(plot1,1:L(plot1)) 0];
ymin1 = min(plotpsi1) - 0.1;
ymax1 = max(plotpsi1) + 0.1;

dt2 = 1/2^plot2;
time2 = -dt2:dt2:L(plot2)*dt2;
xmin2 = -dt2;
xmax2 = L(plot2)*dt2+2*dt2;
plotpsi2 = [0 psi(plot2,1:L(plot2)) 0];
ymin2 = min(plotpsi2) - 0.1;
ymax2 = max(plotpsi2) + 0.1;

dt3 = 1/2^plot3;
time3 = -dt3:dt3:L(plot3)*dt3;
xmin3 = -dt3;
xmax3 = L(plot3)*dt3+2*dt3;
plotpsi3 = [0 psi(plot3,1:L(plot3)) 0];
ymin3 = min(plotpsi3) - 0.1;
ymax3 = max(plotpsi3) + 0.1;

dt4 = 1/2^plot4;
time4 = -dt4:dt4:L(plot4)*dt4;
xmin4 = -dt4;
xmax4 = L(plot4)*dt4+2*dt4;
plotpsi4 = [0 psi(plot4,1:L(plot4)) 0];
ymin4 = min(plotpsi4) - 0.1;
ymax4 = max(plotpsi4) + 0.1;

figure(1)

subplot(4,2,2)
plot([xmin1 xmax1], [0 0],'k-')
hold on
stairs(time1,plotpsi1)
axis([xmin1 xmax1 ymin1 ymax1]); 
xlabel('Time')
ylabel(strcat('Amplitude psi(',int2str(plot1),')'))

subplot(4,2,4)
plot([xmin2 xmax2], [0 0],'k-')
hold on
stairs(time2,plotpsi2)
axis([xmin2 xmax2 ymin2 ymax2]); 
xlabel('Time')
ylabel(strcat('Amplitude psi(',int2str(plot2),')'))

subplot(4,2,6)
plot([xmin3 xmax3], [0 0],'k-')
hold on
stairs(time3,plotpsi3)
axis([xmin3 xmax3 ymin3 ymax3]); 
xlabel('Time')
ylabel(strcat('Amplitude psi(',int2str(plot3),')'))

subplot(4,2,8)
plot([xmin4 xmax4], [0 0],'k-')
hold on
stairs(time4,plotpsi4)
axis([xmin4 xmax4 ymin4 ymax4]); 
xlabel('Time')
ylabel(strcat('Amplitude psi(',int2str(plot4),')'))
end

