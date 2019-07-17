[X,map] = imread('Lena.bmp');
Lena = double(ind2gray(X,map));
save ('Lena.mat','Lena')