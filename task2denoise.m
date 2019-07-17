function denoise(decomp_times,SNR,type,wname)

[LoD,HiD,LoR,HiR] = wfilters(wname); 

%Coefficients for the db4 wavelet
filter =LoD;
[X,map] = imread('Lena.bmp');
Lena = double(ind2gray(X,map));
[M,N] = size(Lena);

%Low-pass filter
h_lp = filter;
g_lp = wrev(h_lp);

%High-pass filter
h_hp = -( (-1).^(1:length(filter)) ).*wrev(h_lp);
g_hp = wrev(h_hp);


lf = length(filter);

%Calculate noise from SNR input and add it to the image
var_s = (std2(Lena))^2;
var_n = 10^(log10(var_s) - (SNR/10));
eta = sqrt(var_n)*randn(M,N);
y = Lena + eta;

%Calculate standard deviation of the noise which will be used to...
%...calculate the threshold value used to denoise the image.
sigma = sqrt(var_n);

%Calculate hard threshold from S. Mallet, "A Wavelet Tour of Signal Processing", pp. 462,
%Academic Press, 1999, 2nd edition.
if strcmp(type,'H')%If using hard thresholding, calculate threshold as follows:
   thr = 3*sigma;
end

%Calculate soft threshold from S. Mallet, "A Wavelet Tour of Signal Processing", pp. 462,
%Academic Press, 1999, 2nd edition.
if strcmp(type,'S')
   thr = 3*sigma/2;
end
disp("thr ="+thr) %check the thr value

%Loop for decomposing Lena
for k = 1:decomp_times
   
	%Zero-pad the image
	%Extend left and right of image with zeros
	[M1,N1] = size(y);
	right_ext = zeros(M1,lf-1);
	left_ext = zeros(M1,lf-1);
	y = [left_ext y right_ext];
	%Extend upper and lower image with zeros
	[M2,N2] = size(y);
	upper_ext = zeros(lf-1,N2);
	lower_ext = zeros(lf-1,N2);
	y = [upper_ext; y; lower_ext];

	%Low-pass filter signal and down-sample
   %Convolve rows
   [M3,N3] = size(y);
	for j = 1:M3
   	approx{k}(j,:) = conv(y(j,:),h_lp);
	end
	[A,B] = size(approx{k});
	for j = 1:B
   	b1{k}(:,j) = conv(approx{k}(:,j)',h_lp)';
   end
	for j = 1:B
   	horiz_details{k}(:,j) = conv(approx{k}(:,j)',h_hp)';
   end
   s = (B-(M1+lf-1))/2;
   first = floor(s)+ 1;
   last = B-ceil(s);
   

	%Down-sample...
	a{k} = b1{k}(first:2:last,first:2:last);  %Save every other column and row for first approximation
   h{k} = horiz_details{k}(first:2:last,first:2:last);  %Save every other column and row for 1st horizontal details
   
      
	%High-pass filter signal and down-sample
	%Convolve rows
	for j = 1:M3
   	details{k}(j,:) = conv(y(j,:),h_hp);
	end
	for j = 1:B
   	c1{k}(:,j) = conv(details{k}(:,j)',h_hp)';
	end
	for j = 1:B
   	vert_details{k}(:,j) = conv(details{k}(:,j)',h_lp)';
	end
	%Down-sample...
	d{k} = c1{k}(first:2:last,first:2:last);  %Save every other column and row for first diagonal details
   v{k} = vert_details{k}(first:2:last,first:2:last);  %Save every other column and row for first vertical details

   y = a{k};
   [Q, R] = size(y);
   q{k} = [Q R];
   
   %Soft thresholding
   if strcmp(type,'S')
      x = abs(h{k});
      h{k} = sign(h{k}).*(x >= thr).*(x - thr);
      x = abs(d{k});
      d{k} = sign(d{k}).*(x >= thr).*(x - thr);
      x = abs(v{k});
     	v{k} = sign(v{k}).*(x >= thr).*(x - thr); 
	   end
   
   %Hard thresholding
   if strcmp(type,'H')
      h{k} = h{k}.*(abs(h{k})>thr);
      d{k} = d{k}.*(abs(d{k})>thr);
		v{k} = v{k}.*(abs(v{k})>thr);

   end
end


%Loop for reconstructing Lenna
for k = decomp_times:-1:1
   %Upsample and filter approximations
   M = q{k}(1);
   N = q{k}(2);
	upsamp_a{k} = zeros(2*M,2*N);
   for m = 1:M
      for n = 1:N
         upsamp_a{k}(2*m,2*n) = y(m,n);
      end
   end
   for j = 1:2*M
   	col_a{k}(:,j) = conv(upsamp_a{k}(:,j)',g_lp)';
	end
	[A,B] = size(col_a{k});
	for j = 1:A
   	rec_a{k}(j,:) = conv(col_a{k}(j,:),g_lp);
   end
   rec_a{k} = rec_a{k}(lf+1:end,lf+1:end);
   %Upsample and filter horizontal details
   upsamp_h{k} = zeros(2*M,2*N);
   for m = 1:M
      for n = 1:N
         upsamp_h{k}(2*m,2*n) = h{k}(m,n);
      end
   end
   for j = 1:2*M
   	col_h{k}(:,j) = conv(upsamp_h{k}(:,j)',g_hp)';
	end
	for j = 1:A
   	rec_h{k}(j,:) = conv(col_h{k}(j,:),g_lp);
	end
   rec_h{k} = rec_h{k}(lf+1:end,lf+1:end);
   %Upsample and filter diagonal details
   upsamp_d{k} = zeros(2*M,2*N);
   for m = 1:M
      for n = 1:N
         upsamp_d{k}(2*m,2*n) = d{k}(m,n);
      end
   end
   for j = 1:2*M
   	col_d{k}(:,j) = conv(upsamp_d{k}(:,j)',g_hp)';
	end
	for j = 1:A
   	rec_d{k}(j,:) = conv(col_d{k}(j,:),g_hp);
	end
   rec_d{k} = rec_d{k}(lf+1:end,lf+1:end);
   %Upsample and filter vertical details
   upsamp_v{k} = zeros(2*M,2*N);
   for m = 1:M
      for n = 1:N
         upsamp_v{k}(2*m,2*n) = v{k}(m,n);
      end
   end
   for j = 1:2*M
   	col_v{k}(:,j) = conv(upsamp_v{k}(:,j)',g_lp)';
	end
	for j = 1:A
   	rec_v{k}(j,:) = conv(col_v{k}(j,:),g_hp);
	end
	rec_v{k} = rec_v{k}(lf+1:end,lf+1:end);

	y = rec_a{k} + rec_h{k} + rec_d{k} + rec_v{k};
end

recon_image = y(1:end-lf+1,1:end-lf+1);
var_noise = (std2(Lena - recon_image))^2;
SNR = 10*log10(var_s/var_noise);
SNR_orig = 10*log10(var_s/var_n);

f=figure;
S = get(0,'MonitorPosition');
set(f,'position',[S(3)*0.3/2 S(4)*0.5/2 S(3)*0.7 S(4)*0.5]);

subplot(131)
imshow(recon_image,[]);
title(['Reconstructed Image, SNR = ',num2str(SNR),'dB']);
subplot(132)
imshow(Lena,[]);
title('Original Image');
subplot(133)
imshow(Lena+eta,[]);
title(['Original Noise Image, SNR = ',num2str(SNR_orig),'dB']);

%noisyLena=Lena+eta
%save ('Lena_plus_noise.mat','noisyLena')
