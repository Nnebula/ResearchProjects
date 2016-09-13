function f=nONbackGenerator(nTOT)
    function dn = nONback(t,n)

        kOFF0 = 0.1; % s-1    -bond detachment rate at focal adhesions (no force)
        kON = 1; % s-1        -bond attachment rate
        F0 = 1; % pN          -bond strength 
        F = 30; % pN           -*TODO: change to your model
        %nTOT = 100; %         -Total number of focal adhesions 

        dn = kOFF0*exp(F/(F0*n))*n - kON*(nTOT-n); % - nON
    end

    f=@nONback;
end
