% invasion of cell into extra cellular matrix of collagen and fibronectin
% mathematical cell biology (math 563)
% supervisor: Leah Keshet
% project idea: Jun Allrad, Naghmeh Rezaei

cla
clc

% --constant variables
kT = 4.11; % pNnm

% --model parameters
Koff0 = 0.1; % s-1    -bond detachment rate at focal adhesions (no force)
Kon = 1; % s-1        -bond attachment rate
F0 = 2; % pN          -bond strength 
Fmat = 2; % pN        -*TODO: change to your model
nTotal = 75; %        -number of available surface receptors for focal adhesions
nFocal = 50;
dt = 0.01; % s        -simulation time step
x = 0; %              -initial position 


% --actin parameters
c = 10; % uM                -actin monomer concentration *TODO:check the number
d = 2.72; % nm              -actin monomer contribution to filament length
KonActin = 11.6; % 1/(s.uM) -actin polymerization rate 
KoffActin = 1.4; % s-1      -actin depolymerization rate

% --actin calculations
c0 = KoffActin/KonActin; % actin critical concentration
Fmax = (kT/d)*log(c/c0); % maximum force generated by one actin filament
n0 = Fmat/Fmax; %          minimum number of actin filaments required for invasion
nActin = n0+10; %         number of pushing actin filaments *TODO: calculate the number

% --dynamic variables
adhesions = ones(1,nTotal); % state of each surface bindinpNg site, 0:un-attached, 1:attached

% --setup figures
subplot(3,1,1); hold on; box on;
subplot(3,1,2); hold on; box on;
subplot(3,1,3); hold on; box on;

% --time loop starts
for t=0:dt:1000 
    
    % --Attachment
    pAttach = 1-exp(-Kon*dt); %                                   probability of one un-attached focal adhesion to attach
    attachThese = (1-adhesions).*( rand(1,nTotal) < pAttach); %   returns an array of false and true with probability of pAttach
      
    %--Detachment 
    Koff = Koff0*exp(Fmat/(F0*nFocal)); %        force-dependent detachment rate (Bell's law)
    pDetach = 1-exp(-Koff.*dt); %                probability of one focal adhesion to detach
    detachThese = adhesions.*(rand(1,nTotal) < pDetach)
    
																				 adhesions (detachThese) = 0;
																				 adhesions (attachThese) = 1; 

    nFocal = sum(adhesions); %                   number of focal adhesions (attached)
    
    if (1)
        
        v = d*(KonActin*c*exp(Fmat*d/(kT*nActin))-KoffActin); %nm/s    velocity of invadopodia
        x = v*dt+x; % nm   length of invadopodia in ECM
        
        %--plot results every few time steps
        if (mod(t,10)==0)
            
            subplot(3,1,1);
            plot(t,v,'bx');
            title('velocity');
            xlabel('time');
            ylabel('velocity');
            
            subplot(3,1,2);
            plot(t,x,'r+');
            title('position');
            xlabel('time');
            ylabel('position');
            
            subplot(3,1,3);
            plot(t,nFocal,'g*');
            title('focal adhesions');
            
            drawnow;
            
        end % end of plot
        
    end % end of if(nFocal*F0>Fmat)

end % end of time loop 


