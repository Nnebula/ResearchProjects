function [ r ] = FindCircle ( Nim )
%dx2 = (-319:320).^2;
%dy2 = (-239:240).^2;
%dist = sqrt(ones(480,1)*dx2 + dy2'*ones(1,640));

h = [];
for x=1:640
    for y=1:480
        d = round(sqrt((x-320)^2+(y-240)^2));
        if d==0
            continue
        else if length(h)<d
                h(d) = Nim(y,x);
            else
                h(d) = h(d) + Nim(y,x);
            end
        end
    end
end

plot(h);
h2 = smooth(h,25);
r=max(find(h2>20));

% X=sum(Nim,1);
% Y=sum(Nim,2);
% Xwin=320-min(find(smooth(X,32)>20));
% Ywin=240-min(find(smooth(Y,32)>20));
% cla
% imshow(Nim);
% rectangle('Position', [320-Xwin 240-Ywin Xwin*2 Ywin*2], 'EdgeColor', 'b')
% w = Xwin;
% h = Ywin;
end

