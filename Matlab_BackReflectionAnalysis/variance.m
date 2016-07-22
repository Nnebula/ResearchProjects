function variance(p)

files = dir([p ,'/*.tif']);
x1 = (1:640);
X = repmat(x1,480,1);
X2 = X.*X;
varx=[];
objpos=[];

for i = 1:length(files)
    f=files(i).name;
    filename = strrep(f,'.tif', '');
    im=imread([p, '/', f]);
    Nim=double(im)/255;
    s=sum(Nim(400:480,:))/81;
    Nim2=Nim-ones(480,1)*s;
    
    objpos(end+1)=str2num(strrep(filename,'mV',''));
    
    wX2=X2.*Nim2; 
    varx(end+1)=sum(sum(wX2))/sum(sum(Nim));
end

f6=figure(6);
plot(objpos/1000,sqrt(varx));

end