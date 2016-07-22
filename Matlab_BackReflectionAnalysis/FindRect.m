function [ w h found] = FindRect( Nim )
X=smooth(sum(Nim,1),32);
Y=smooth(sum(Nim,2),32);
found = true;

Xind1 = find(X(1:end/2)>20, 1);
if isempty(Xind1)
    Xind1 = 0;
    found = false;
else
    Xind1 = Xind1(1);
end

Xind2 = find(X(end/2+1:end)>20, 1, 'last');
if isempty(Xind2)
    Xind2 = length(X);
    found = false;
else
    Xind2 = Xind2(1)+length(X)/2;
end

Yind1 = find(Y(1:end/2)>20, 1);
if isempty(Yind1)
    Yind1 = 0;
    found = false;
else
    Yind1 = Yind1(1);
end

Yind2 = find(Y(end/2+1:end)>20, 1, 'last');
if isempty(Yind2)
    Yind2 = length(Y);
    found = false;
else
    Yind2 = Yind2(1)+length(Y)/2;
end

%[Xind1 Yind1 Xind2 Yind2]

w = Xind2 - Xind1;
h = Yind2 - Yind1;
imshow(Nim);
if(found)
    rectangle('Position', [Xind1 Yind1 w h], 'EdgeColor', 'b', 'Curvature',[1,1]);
end

end

