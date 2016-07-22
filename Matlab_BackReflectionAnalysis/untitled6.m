
320-min(find(smooth(X1,32)>20))
320-min(find(smooth(X2,32)>20))
320-min(find(smooth(X3,32)>20))
320-min(find(smooth(X4,32)>20))

X=X4;
win=320-min(find(smooth(X,32)>10));
cla
imshow(Nim4);
rectangle('Position', [320-win 240-win win*2 win*2], 'EdgeColor', 'b')


X=X1;
win=320-min(find(smooth(X,32)>10));
cla
imshow(Nim1);
rectangle('Position', [320-win 240-win win*2 win*2], 'EdgeColor', 'b')
