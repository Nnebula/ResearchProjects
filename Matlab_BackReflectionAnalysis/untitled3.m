im1=imread('111119-50000mV.tif');
im2=imread('111519-51000mV.tif');
im3=imread('56344-13000mV.tif');
im4=imread('65277-26000mV.tif');
Nim1=double(im1)/255;
Nim2=double(im2)/255;
Nim3=double(im3)/255;
Nim4=double(im4)/255;
X1=sum(Nim1,1);
Y1=sum(Nim1,1);
X2=sum(Nim2,1);
Y2=sum(Nim2,2);
X3=sum(Nim3,1);
Y3=sum(Nim3,2);
X4=sum(Nim4,1);
Y4=sum(Nim4,2);

cla
plot(X1, 'r');
hold on
plot(X2, 'b');

cla
plot(smooth(X1,60), 'r');
hold on
plot(smooth(X2,60), 'b');

