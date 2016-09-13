function dn = nON(t,n);  

kOFF0 = 0.1; % s-1    -bond detachment rate at focal adhesions (no force)
kON = 1; % s-1        -bond attachment rate
F0 = 2; % pN          -bond strength 
F = 60; % pN          -*TODO: change to your model
nTOT = 75; %          -Total number of focal adhesions 

dn = -kOFF0*exp(F/(F0*n))*n + kON*(nTOT-n); % n=number of attached focal adhesions or bonds (n_on)

end