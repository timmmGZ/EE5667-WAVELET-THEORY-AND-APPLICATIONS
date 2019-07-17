 clear
S = get(0,'MonitorPosition');
w=S(3)*0.6
h=S(4)*0.45
set(figure,'position',[(S(3)-w)/2 (S(4)-h)/2 w h]);
load Lena1
[XDEN,cfsDEN,dimCFS] = func_denoise_dw2d(X);
colormap(gray);
subplot(121);
imagesc(X);
title('Image denoised in the GUI');
subplot(122);
imagesc(XDEN); 
title('Image denoised with generated code');
norm(XDEN-X,2)
