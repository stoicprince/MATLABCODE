%% Running the modelling

% Define the range of time 
range = (0:20);

% Define the initial conditions of ODEs
ICs = [1.0e10, 2.5e8, 5.268e5, 2.25e9, 0, 0, 1073];

% call the user-defined functions
% pass the range of independent variables
% pass the intitial conditions
[tsol,varsol] = ode45(@ode_sys, range, ICs);
 
%% Results

tsol=zeros(40,100)
varsol=zeros(40,100)

tsol(i,:)=tsol:
varsol(i,:)=varsol(:,1)
End

plot(tsol,varsol(:,[1,2,3,4]))


%% Define the functions

function diffeqs = ode_sys(t,var)

% define the unknown variables

    TS = var(1);    % tumour cell population
    NS = var(2);    % concentration of NK cells
    LS = var(3);    % concentration of CD8+T-cells
    CS = var(4);    % concentration of lymphocytes 
    PS = var(5);    % concentration of nanoparticles
    MS = var(6);    % concentration of chemotherapeutic drugs
    IS = var(7);    % concentration of IL-2
    
% define the model parameters

    ap  = 4.31e-1;
    bp  = 1.02e-9;
    cp  = 2.9077e-13;
    KTp = 9.0e-1;
    sigmaTp =1.8328;
    fp  = 1.25e-2;
    ep  = 1.11e-1*1.25e-2;
    pp  = 2.794e-13;
    pNp = 6.68e-2;
    gNp =2.5036e5;
    KNp = 6.75e-2;
    sigmaNp = 1.8328;
    theta = 2.5036e-3;
    mp  = 9.0e-3;
    jp  = 1.245e-2;
    kp  = 2.019e7;
    qp  = 3.422e-10;
    r1p = 2.9077e-11;
    r2p = 5.8467e-13;
    up  = 4.417e-14;
    kapap = 2.5036e3;
    KLp = 4.86e-2;
    sigmaLp = 1.8328;
    pIp = 2.971;
    gIp = 2.5036e3;
    alphap = 2.25e-1*6.3e-3;
    betap = 6.3e-3;
    KCp = 3.4e-2;
    sigmaCp = 1.8328;
    gamma = 5.199e-1;
    muIp = 11.7427;
    phip = 2.38405e-7;
    omegap = 7.874e-2;
    xip  = 2.5036e3;
    dp  = 2.34;
    lp  = 2.09;
    sp  = 3.8e-3;
    krelp = 1.0e-4;
    keN = 2.23e-4;
 
    if (t<=3600) 
        RLp = 1.77e10;
    elseif (t>3600)
        RLp = 0;            % this defines infusing CD8+T-cells in one hour.
    end
    
  if (t<=3600)
        RPp = 2.3869;
    elseif (t>3600)
        RPp = 0;            % this defines infusing nanoparticles in one hour.

    end

    if (t<=3600)
        RIp = 2.7859e6;
    elseif (t>3600)
        RIp=0;              % this defines infusing IL-2 in one hour.
    end

if (t<=3600)                     %defines chemotherapy infusion rate
        IF=02;
  elseif (t>3600)
      IF=0; 

  end 

    
    DS = dp*power(LS/TS,lp)/(sp+power(LS/TS,lp));
      
% define the ODEs    
    
    diffeqs(1,1) = ap*TS*(1-bp*TS)-cp*NS*TS-DS*TS-KTp*(1-exp(-sigmaTp*MS))*TS;
    diffeqs(2,1) = fp*(ep/fp*CS-NS)-pp*NS*TS+pNp*NS*IS/(gNp+IS)-KNp*(1-exp(-sigmaNp*MS))*NS;
    diffeqs(3,1) = theta*mp*LS/(theta+IS)+jp*TS*LS/(kp+TS)-qp*LS*TS+(r1p*NS+r2p*CS)*TS-up*LS*LS*CS*IS/(kapap+IS)-KLp*(1-exp(-sigmaLp*MS))*LS+pIp*LS*IS/(gIp+IS)+RLp;
    diffeqs(4,1) = betap*(alphap/betap-CS)-KCp*(1-exp(-sigmaCp*MS))*CS;
    diffeqs(5,1) = -krelp*PS - keN*PS + RPp;
    diffeqs(6,1) = -gamma*MS + krelp*PS+IF;
    diffeqs(7,1) = -muIp*IS+phip*CS+omegap*LS*IS/(xip+IS)+RIp;

end