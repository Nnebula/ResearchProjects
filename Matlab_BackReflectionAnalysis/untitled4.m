
factors=(abs([-320:319])/320.0).^2;
sum(smooth(X1,60)'.*factors)
sum(smooth(X2,60)'.*factors)
sum(smooth(X3,60)'.*factors)
sum(smooth(X4,60)'.*factors)

sum(X1.*factors)
sum(X2.*factors)
sum(X3.*factors)
sum(X4.*factors)
