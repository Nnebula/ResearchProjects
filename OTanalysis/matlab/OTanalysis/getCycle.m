function cycleIndex = getCycle(data);

figure
plot(data);

[x1 y1] = ginput(1);

if x1 <0 || x1> length(data)
    cycle = [];
    return
end

[x2 y2] = ginput(1);

cycleIndex = [ceil(x1)  floor(x2)];

close;
end

