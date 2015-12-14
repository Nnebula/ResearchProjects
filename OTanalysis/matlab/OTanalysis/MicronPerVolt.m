function [NyPSD_um] = MicronPerVolt(yTrap, NyPSD);

p = polyfit(NyPSD,yTrap,4);
NyPSD_um = polyval(p,NyPSD)-p(5);
plot(NyPSD,yTrap,'*',NyPSD,NyPSD_um,'-');

end
