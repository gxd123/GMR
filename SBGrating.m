% STOP BAND OPTIMIZATION
% SBGrating.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% UNITS
degrees     = pi/180;
micrometers = 1;
dB          = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOURCE PARAMETERS
wl = 501;                                   % WAVELENGTH SWEEP RESOLUTION
lam0 = linspace(1.2,1.9,wl) * micrometers;
lamc = 1.55 * micrometers;                  % CENTER WAVELENGTH 
thresh = -20 * dB;                          % THRESHOLD VALUE
SRC.theta = 0;
SRC.phi = 0;
SRC.pte = 1;
SRC.ptm = 0;

% DEVICE PARAMETERS
DEV.er1 = sqrt(10);                        % REFLECTION PERMITTIVITY
DEV.ur1 = 1;                        % REFLECTION PERMEABILITY
DEV.er2 = sqrt(5);                % TRANSMISSION PERMITTIVITY
DEV.ur2 = 1;                        % TRANSMISSION PERMEABILITY
nlay = 3;                           % # OF LAYERS
DEV.NP = 1;

% ASSIGN PERMEABILITIES TO LAYERS
DEV.UR = ones(1,nlay);

% OPTIMIZATION OPTIONS
iter = 10000;                         % # OF ITERATIONS TO LOOP THROUGH
nsol = 5;                          % # OF SOLUTIONS TO STORE
nmin = 1;                           % MIN REFRACTIVE INDEX
nmax = 2;                           % MAX REFRACTIVE INDEX
dmin = 0.001 * micrometers;                  % MIN LAYER SIZE
dmax = max(lam0) * micrometers;              % MAX LAYER SIZE

% SAVING OPTIONS
save_flag = 1;                      % SAVE INFORMATION TO TEXT FILE
filep = '.\Saved Data\';
filen = 'HW11';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM RANDOM OPTIMIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE GLOBAL PARAMETERS
AREA        = zeros(nsol,1);
r_len       = AREA;
REFi        = AREA;
ind_REFi    = AREA;
REFm        = AREA;
g_ind       = zeros(nsol,2);
n           = zeros(nsol,nlay);
d           = n;

% INITIALIZE DEVICE PARAMETERS
ER  = rand(iter,nlay)*(sqrt(nmax)-sqrt(nmin)) + sqrt(nmin);
L   = rand(iter,nlay)*(dmax-dmin) + dmin;

% LOOP THROUGH OPTIMIZATION 
for iteration = 1:iter
    
    % DISPLAY CURRENT ITERATION
    clc;
    disp(['Current Iteration is ' num2str(iteration)]);
    
    % INTIALIZE ARRAYS
    REF     = zeros(wl,1);
    TRN     = zeros(wl,1);
    DEV.ER  = ER(iteration,:);
    DEV.L   = L(iteration,:);

    % LOOP THROUGH THE VISIBLE SPECTRUM
    for w = 1:wl
        SRC.lam0 = lam0(w);
        DAT = tmm1d(DEV,SRC);
        REF(w) = DAT.REF;
        TRN(w) = DAT.TRN;
    end

    % CALCULATE CONSERVATION AND CONVERT TO dB
    CON = REF+TRN;
    CON = 10*log10(CON);
    REFdB = 10*log10(REF);

    % DETERMINE MERIT FUNCTION INDICES
    % LOOP THROUGH POINTS ON THE LEFT HAND SIDE 
    lamci = find(lam0 == lamc);         % CENTER WAVELENGTH INDEX
    ind = [0 0];                        % INDEX FOR -20dB
    Rm = [-1000 -1000];                 % MAX REFLECTANCE IN INTERVAL 
    dif = diff(REFdB);                  % DIFFERENCES BETWEEN POINTS
    for w = 1:lamci-2
        cdif = sign(dif(lamci-w));      % CURRENT DIFFERENCE
        Rc = REFdB(lamci-w);            % CURRENT REFLECTANCE
        ndif = sign(dif(lamci-w-1));    % NEXT DIFFERENCE
        if cdif==-1 && ndif==-1 && Rc>thresh
            ind(1) = lamci-w+1;
            break;
        elseif cdif==-1 && ndif==1
            Rm(1) = max(Rm(1),Rc);
        else 
            ind(1) = lamci-w-1;
        end
    end

    % LOOP THROUGH POINTS ON THE RIGHT HAND SIDE 
    for w = 1:lamci-2
        cdif = sign(dif(lamci+w-1));    % CURRENT DIFFERENCE
        ndif = sign(dif(lamci+w));      % NEXT DIFFERENCE
        Rn = REFdB(lamci+w);            % NEXT REFLECTANCE   
        if cdif==1 && ndif==1 && Rn>thresh
            ind(2) = lamci+w-1;
            break;
        elseif cdif==1 && ndif==-1
            Rm(2) = max(Rm(2),Rn);
        else 
            ind(2) = lamci+w+1;
        end
    end

    % RECTANGLE ALGORITHM
    R = REFdB(ind(1):ind(2));   % REDUCE SAMPLE SIZE
    area = 0;
    for w = 1:length(R)
        i = R<R(w);
        l = 0;
        % LOOP TO LEFT OF CURRENT POSITION
        for j = 1:w-1
            if l == 0 && i(w-j) == 1
                l = l + 2;
            elseif i(w-j) == 1
                l = l + 1;
            else 
                break;
            end
        end
        % LOOP TO RIGHT OF CURRENT POSITION
        for j = 1:length(R)-w
            if j == 1 && i(w+j) == 1
                l = l + 2;
            elseif i(w+j) == 1
                l = l + 1;
            else
                break;
            end
        end

        % FIND AREA OF RECTANGLE AND STORE MAXIMUM AREA 
        A = l * abs(R(w));
        if A > area
            area = A;
            l_L = l;
            R_L = R(w);
            ind_L = ind(1) + w - 1;
        end

    end

    % STORE TO GLOBAL PARAMETERS
    if area > AREA(end)
        % STORE VALUES
        AREA(end) = area;
        r_len(end) = l_L;
        REFi(end) = R_L;
        ind_REFi(end) = ind_L;
        REFm(end) = max(Rm);
        g_ind(end,:) = ind;
        n(end,:) = DEV.ER.^2; 
        d(end,:) = DEV.L;

        % SORT VALUES FROM BIGGEST TO SMALLEST AREA
        [AREA,I] = sort(AREA,'descend');
        r_len = r_len(I);
        REFi = REFi(I);
        ind_REFi = ind_REFi(I);
        REFm = REFm(I);
        g_ind = g_ind(I,:);
        n = n(I,:); 
        d = d(I,:);
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% POST-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MAXIMUM REFLECTANCE IN THE BAND
fprintf('Maximum Reflectance is %4.2f dB \n',-max(abs(REFm)));

% SAVE DATA
if save_flag
    % OPEN TEXT FILE 
    fid = fopen([filep filen '.txt'],'a+');
    
    % ORGANIZE DATA TO WRITE
    for w = 1:nsol
        line = '\n  %d \t %6.4f \t %d \t %4.4f \t\t %d \t\t\t %4.4f \t\t %d,  %d \t\t\t\t';
        line = sprintf(line,iter,AREA(w),r_len(w),REFi(w),ind_REFi(w),REFm(w),g_ind(w,1),g_ind(w,2));
        for iteration = 1:nlay
            line = sprintf([line '%4.4f,  '],n(w,iteration));
        end
        line = sprintf([line ' \t\t\t']);
        for iteration = 1:nlay
            line = sprintf([line '%4.4f,  '],d(w,iteration));
        end
        fprintf(fid,line);
    end
    fprintf(fid,'\n');
    
    % CLOSE TEXT FILE 
    fclose(fid);
    
end
