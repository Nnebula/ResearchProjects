function [ force_z0 ] = findZero( force )
%findZero Find force Zero points
    f_smooth = smooth(force,5000);
    %plot(f_smooth);
    f_derv = diff(f_smooth)*100;
    f_derv(size(f_derv,1)+1)=f_derv(size(f_derv,1));
    smooth(f_derv, 2000);
    %hold all;
    %plot(f_derv);
    f2=f_smooth;
    f2(abs(force) > 10 & f_derv>-0.5) = -1000;
    %hold all;
    %plot(f2);
    ind=find(f2==-1000, 1, 'last');
    f3=f_smooth;
    f3(abs(force) > 5 | abs(f_derv)>0.15) = -1000;
    %hold all;
    %plot(f3);
    f3=f3(ind:end);
    f3(f3==-1000) = [];
    force_z0=median(f3);
end

