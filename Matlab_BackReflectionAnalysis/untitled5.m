
modelFun = @(p, x) p(3)*p(2)./(pi*(p(2).^2+(x-p(1)).^2));
startingVals = [1 1 1];
options = statset('Display', 'iter', 'TolFun', 1e-16)
coefEsts = nlinfit(1:320, X1(321:640), modelFun, startingVals, options);

%modelFun = @(p, x) p(3)*poisspdf(x-round(p(2)), p(1));
modelFun =  @(p,x) p(3) .* (x ./ p(1)).^(p(2)-1) .* exp(-(x ./ p(1)).^p(2));
startingVals = [10 1 1];
%options = statset('Display', 'iter', 'TolFun', 1e-16, 'TolX', 1e-16)
options = statset('TolFun', 1e-16, 'TolX', 1e-16)
coefEsts = nlinfit(1:320, X1(321:640), modelFun, startingVals, options);
coefEsts(2)
options = statset('TolFun', 1e-16, 'TolX', 1e-16)
coefEsts = nlinfit(1:320, X2(321:640), modelFun, startingVals, options);
coefEsts(2)
options = statset('TolFun', 1e-16, 'TolX', 1e-16)
coefEsts = nlinfit(1:320, X3(321:640), modelFun, startingVals, options);
coefEsts(2)
options = statset('TolFun', 1e-16, 'TolX', 1e-16)
coefEsts = nlinfit(1:320, X4(321:640), modelFun, startingVals, options);
coefEsts(2)

cla
hold on;
plot(X4, 'r');
plot(321:640, modelFun(coefEsts, 1:320), 'b');
hold off;

coefEsts(3)
