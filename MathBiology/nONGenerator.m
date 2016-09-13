function f=nONGenerator(nTOT)
    function dn = nON(t,n)

        kOFF0 = 0.1; % s-1    -bond detachment rate at focal adhesions (no force)
        kON = 1; % s-1        -bond attachment rate
        F0 = 1; % pN          -bond strength 
        F = 30; % pN           -*TODO: change to your model
        %nTOT = 100; %         -Total number of focal adhesions 

        dn = -kOFF0*exp(F/(F0*(n+1e-5)))*n + kON*(nTOT-n); % n=number of attached focal adhesions or bonds (n_on)

    end

    f=@nON;
end
