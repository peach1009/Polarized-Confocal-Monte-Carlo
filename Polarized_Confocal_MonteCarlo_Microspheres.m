% Polarized Monte Carlo simulation of a scanning confocal reflectance 
% Mueller-matrix microscope imaging a slab of Mie-scattering spheres. 
% 
% What the code does: 
% For each input polarization state (H,V, P, R) and each focal depth listed
% in DIRECTORY, the simulation: 
%     1. Samples a launch point in the objective pupil and send the photon to 
%     the tissue surface from which it will reach the focus
%     2. Propagates the photon through the slab
%     3. Detects if it leaves the top surface (a) within the objective NA,
%     and (b) back-projects in the focal plane to within the pinhole radius 
%     (physical radius / magnification) of the current scan position.
% This is repeated at every point of a 20 x 20 raster scan. The four input
% states give the 16 signals needed for a Mueller Matrix decomposition. 
%
% Units: 
%   Mie inputs (rad, lambda) - micrometers (um)
%   Geometry (positions, f_obj, depths) - centimeters (cm)
%
% Coordinates: 
%   z = 0 at the tissue surface and +z into the tissue. Focus at z = Fdepth
%   inside the tissue. 
%
% Outputs: 
%   outXI, outXQ, outXU, outXV - NN x NN Stokes images of detected light 
%   outXIDeep, outXQDeep, outXUDeep, outXVDeep - NN x NN x NNz Detected
%   Stokes binned by photon's last scattering position. 


clear all
close all

tic
% Start parallel pool
delete(gcp('nocreate')); % close any pool; parfor opens a new one 
rng(123); % random seed

% -------- User Setup ---------------% 
% Folder must hold this script and the helper functions cd('path/to/helperfunctions')
directory = {'D1'};
format long
ndirs = numel(directory);

nre_med = 1.33;  % medium RI (must set before foc_list)

% Base output directory ('' = current folder), where D1,2,..DN depth
% folders are generated 
baseDir = pwd; % folder for output files (pwd = current folder path)
depth_um = [0 100 200 300 400 500 600]; % [µm] focus depths
Fdepth_cm = depth_um * 1e-4;            % [cm] 1 µm = 1e-4 cm
obj_gap_ref = 3.4; % [cm] objective gap with the focus at the surface
f_obj = obj_gap_ref; % focus at interface [cm]
% Objective gap per depth: moving the objective by d/n moves the focus by d
Foc_List    =  obj_gap_ref - Fdepth_cm / nre_med;

%% Scanning Parameters
XpCtr = 0;    % [cm] Pupil center x 
YpCtr = 0;    % [cm] Pupil center y 

% Lateral constraints on sample
boxHalfSizeX = 0.5; % [cm]
boxHalfSizeY = 0.5; % [cm]

% [X0, Y0] scan centers
step = 3e-4;                            % [cm] step size: 3 µm in cm
vals = (-9.5:9.5) * step;               % set number of scan points (20 x 20)
[X, Y] = meshgrid(vals, vals);
scan_pts = [X(:), Y(:)];
nScan = size(scan_pts, 1);

%% Build parameter list: one parfor iteration per (depth folder, input state)
param_list = [];
for ddir = 1:ndirs
    for jjj = 1:4
        param_list = [param_list; ddir, jjj];
    end
end

%% Main parallel loop
parfor idx = 1:size(param_list, 1)


    ddir = param_list(idx, 1);   % depth index
    jjj  = param_list(idx, 2);   % polarization index

    % ----- Output Folder: baseDir / Dk -----
    depthFolderName = directory{ddir};
    outputDir = fullfile(baseDir, depthFolderName);
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    %----------------- Photon and grid parameters -----------------
    Nphotons  = 1e6;      % photons launched per scan position
    THRESHOLD = 0.01;     % Russian roulette weight threshold
    CHANCE    = 0.1;      % Russian roulette survival probability 
    NN        = 100;      % lateral bins (2D images and 3D maps)
    NNz       = 250;      % depth bins (3D maps)
    ALIVE     = 1;
    DEAD      = 0;
    MM        = NN-1;

     %----------------- Mie scattering parameters -----------------
    rad     = 1.0/2; % [um] sphere radius
    lambda  = 0.783; % [um] wavelength
    mua     = 0.0;   % [cm^-1] absorption coefficient
    nre_p   = 1.59;  % sphere refractive index
    nre_med = 1.33;  % medium refractive index (keep equal to the value at the top)
    nangles = 1000;  % angular samples of the scattering matrix 

     %----------------- Detection optics -----------------
    M = 5;          % set magnification of system (5x)
    rphys = 0.03/2; % [cm] physical pinhole RADIUS (300 um diameter)
    pinholeR = rphys/M; % [cm], back-projected pinhole radius

    %----------------- Objective -----------------
    f_obj       = 4;        % [cm] 5x objective, EFL 40 mm
    obj_gap_ref = 3.4;      % [cm] objective gap with the focus at interface
    obj_gap     = Foc_List(ddir); % [cm] objective gap for this depth
    Fdepth      = (nre_med/1) * (obj_gap_ref - obj_gap);  % [cm] focal depth in the tissue
    z0          = Fdepth;
    Z_pupil     = f_obj;    % [cm] distance from pupil plane to the surface
    alpha_fill  = 0.9;      % alpha fill factor (fill factor of the objective pupil
    NA          = 0.14;     % numerical aperture

    % ---------------- Derived quantities ----------------
    R_p = f_obj * NA;       % [cm] pupil radius 
    w_p = alpha_fill * R_p; % [cm] beam waist profile
    Zp  = -Z_pupil;         % [cm] z of the pupil plane

    % Allocate
    xp = 0; yp = xp; zp = xp;
    Hx = xp; Hy = xp; Hz = xp;


    ttt1  = asin(NA/nre_med);  % [rad] largest angle accepted in the tissue (NA gate)
    w = tan(0.5*asin(NA/nre_med))*Fdepth;

    x_param = 2*pi*rad/(lambda/nre_med);
    A = pi*rad*rad;
    conv = 1;

    % Calculate scattering matrix
    [S, C, ang] = calcmie(rad, nre_p, nre_med, lambda, nangles, 'ConvergenceFactor', conv);
    Q = getEfficiencies(C, rad(end), 3);
    M = getMuellerMatrix(S);
    s11 = squeeze(M(1,1,:));
    s12 = squeeze(M(1,2,:));
    s33 = squeeze(M(3,3,:));
    s43 = squeeze(M(4,3,:));

    xxx = 2*pi*rad/(lambda/nre_med);
    aa = mie(nre_p/nre_med, xxx);
    g = aa(8); % anisotropy

    % ---------------- Slab and output grids ----------------
    mus    = 30;     % [cm^-1] scattering coefficient
    albedo = mus/(mus + mua);
    slabsize = 1;    % [cm] slab thickness
    hw = 0.015;      % [cm] half-width of the image grids 
    dx = 2.0*hw/NN;  % [cm] lateral bin
    dy = 2.0*hw/NN;
    dz = slabsize/NNz; % [cm] depth bin of the 3D maps

    n_tissue = nre_med;

    % Input Stokes vector for this iteration
    if jjj == 1
        S0 = [1 1 0 0]; % H 
    elseif jjj == 2
        S0 = [1 -1 0 0];% V
    elseif jjj == 3
        S0 = [1 0 1 0]; % P 
    else
        S0 = [1 0 0 1]; % R 
    end

    % Initialize accumulators ONCE per (depth, pol)
    TH = zeros(4,1); RH = zeros(4,1); % TX: transmitted Stokes totals 
    TV = zeros(4,1); RV = zeros(4,1); % RX: reflected Stokes totals
    TP = zeros(4,1); RP = zeros(4,1);
    TR = zeros(4,1); RR = zeros(4,1);
    IR = zeros(NN,NN); % detected Stokes images 
    QR = zeros(NN,NN);
    UR = zeros(NN,NN);
    VR = zeros(NN,NN);

    IRDeep = zeros(NN,NN,NNz); % detected Stokes by last-scattering position
    QRDeep = zeros(NN,NN,NNz);
    URDeep = zeros(NN,NN,NNz);
    VRDeep = zeros(NN,NN,NNz);

    %-------------- Loop over scan positions --------------
    for s = 1:nScan


        % Scan center for this scan position
        X0 = scan_pts(s,1);   % [cm]
        Y0 = scan_pts(s,2);   % [cm]

        %----------------- PHOTON LOOP -----------------
        for i_photon = 1:Nphotons

            % Launch radius in the pupil 
            u = rand;
            trunc = exp(-2*R_p^2/w_p^2);
            r = (w_p/sqrt(2)) * sqrt(-log(1 - u*(1 - trunc)));
            I_weight = 1;

            psi = 2*pi*rand;
            x_p = XpCtr + r*cos(psi);
            y_p = YpCtr + r*sin(psi);

            % --- Interface hit (geometry to focus z0) ---
            a  = ( r / (f_obj * n_tissue) )^2;
            a  = min(a, 1 - 1e-12);
            Hr = Fdepth * sqrt( a / (1 - a) );
            Hx = X0 + Hr*cos(psi);
            Hy = Y0 + Hr*sin(psi);
            Hz = 0;

            xp = x_p; yp = y_p; zp = Zp;

            % --- Direction cosines ---
            v_air  = [Hx-xp, Hy-yp, Hz-zp];
            v_air  = v_air / norm(v_air);
            U = v_air;  % no Fresnel 

            % Initial Stokes rotation
            S = S0;
            S2 = rotSphi(S, psi);
            S = [1; S2(2)/S2(1); S2(3)/S2(1); S2(4)/S2(1)];

            % Launch at interface
            x = Hx; y = Hy; z = Hz;

            photon_status = ALIVE;
            W = 1; % photon weight 

            %******** ALIVE cycle *****************/
            while (photon_status == ALIVE)

                %*** HOP */
                % Store last scattering position then move
                xlast = x;
                ylast = y;
                zlast = z;
                Ulast = U;

                s_step = -log(rand)/(mus+mua); % [cm] step length
                x = x + U(1)*s_step;
                y = y + U(2)*s_step;
                z = z + U(3)*s_step;

                % Lateral sample boundary
                if abs(x) > boxHalfSizeX || abs(y) > boxHalfSizeY
                    photon_status = DEAD;
                else

                    %*** ABSORB */
                    absorb = W*(1-albedo);
                    W = W - absorb;


                    %---------------- Surface reflection (top) ----------------
                    if (z <= 0)

                        % Back-project the last straight segment to the focal plane z = Fdepth
                        Xfoc_global = (x - (z - Fdepth) * ((x-xlast)/((z-zlast))));
                        Yfoc_global = (y - (z - Fdepth) * ((y-ylast)/((z-zlast))));

                        Xrel = Xfoc_global - X0; % offset from the scan position
                        Yrel = Yfoc_global - Y0;
                        RRad = sqrt(Xrel*Xrel + Yrel*Yrel);

                        ttt = acos(abs(U(3)));  % exit angle in the tissue

                        if (ttt < ttt1)  % Within NA (NA gate)


                            w1 = sqrt((Xrel)^2 + (Yrel)^2);

                            if (w1 <= pinholeR)  % Within pinhole - DETECTED!

                                % Rotate from the meridian frame to the lab x-y frame
                                phi = atan2(U(2),U(1));
                                S2  = rotSphi(S, phi);

                                philast = atan2(Ulast(2),Ulast(1));
                                S2last  = rotSphi(S, philast);
                                                         
                                % Image binning, only inside the hw grid
                                if abs(Xfoc_global) <= hw && abs(Yfoc_global) <= hw
                                    % focal plane image indices
                                    ix = floor( (Xfoc_global + hw) / dx ) + 1;
                                    iy = floor( (Yfoc_global + hw) / dy ) + 1;

                                    % 3D indices: last scattering position (xlast, ylast, zlast)
                                    iyd = 1+floor((ylast + hw)/dy);
                                    ixd = 1+floor((xlast + hw)/dx);
                                    izd = 1+floor((zlast)/dz);
                                    ixd = max(min(ixd, NN), 1); % clamp from 1 to NN
                                    iyd = max(min(iyd, NN), 1);
                                    izd = max(min(izd, NNz), 1); % clamp from 1 to NNz

                                    if ix >= 1 && ix <= NN && iy >= 1 && iy <= NN
                                        % check that indices are within 1:NN

                                        % save Stokes
                                        IR(iy,ix) = IR(iy,ix) + W*S2(1);
                                        QR(iy,ix) = QR(iy,ix) + W*S2(2);
                                        UR(iy,ix) = UR(iy,ix) + W*S2(3);
                                        VR(iy,ix) = VR(iy,ix) + W*S2(4);


                                        % save 3D Stokes
                                        IRDeep(iyd,ixd,izd) = IRDeep(iyd,ixd,izd) + W*S2last(1);
                                        QRDeep(iyd,ixd,izd) = QRDeep(iyd,ixd,izd) + W*S2last(2);
                                        URDeep(iyd,ixd,izd) = URDeep(iyd,ixd,izd) + W*S2last(3);
                                        VRDeep(iyd,ixd,izd) = VRDeep(iyd,ixd,izd) + W*S2last(4);
                                    end
                                end


                            end
                        end

                        photon_status = DEAD;
                    end

                    %---------------- Bottom transmission ----------------
                    if (z >= slabsize)

                        phi = -atan2(U(2),U(1));
                        S2  = rotSphi(S, phi);

                        if jjj == 1
                            TH = TH + S2*W;
                        elseif jjj == 2
                            TV = TV + S2*W;
                        elseif jjj == 3
                            TP = TP + S2*W;
                        else
                            TR = TR + S2*W;
                        end

                        photon_status = DEAD;
                    end

                    %---------------- SPIN - REJECTION METHOD ----------------
                    if photon_status == ALIVE

                        condition = 1;

                        while condition

                            theta = acos(2*rand-1);
                            phi   = rand*2.0*pi;

                            I0      = s11(1)*S(1) + s12(1)*(S(2)*cos(2*phi) + S(3)*sin(2*phi));
                            ithedeg = 1+floor(theta*nangles/pi);
                            I       = s11(ithedeg)*S(1) + s12(ithedeg)*(S(2)*cos(2*phi)+S(3)*sin(2*phi));

                            xra = rand*I0;
                            if (xra >= I)
                                condition = 1;
                            else
                                condition = 0;
                            end
                        end

                        U2       = updateU(U, phi, theta);
                        costheta = cos(theta);

                        % Scattering Plane
                        S2 = rotSphi(S, phi);

                        % Status of polarization after scattering
                        S(1) =  s11(ithedeg)*S2(1) + s12(ithedeg)*S2(2);
                        S(2) =  s12(ithedeg)*S2(1) + s11(ithedeg)*S2(2);
                        S(3) =  s33(ithedeg)*S2(3) + s43(ithedeg)*S2(4);
                        S(4) = -s43(ithedeg)*S2(3) + s33(ithedeg)*S2(4);

                        % Angle gamma referenced to meridian plane
                        temp = (sqrt(1-costheta*costheta)*sqrt(1-U2(3)*U2(3)));
                        if (temp == 0)
                            cosi = 0;
                        else
                            if ((phi>pi) && (phi<2*pi))
                                cosi = (U2(3)*costheta-U(3))/temp;
                            else
                                cosi = -(U2(3)*costheta-U(3))/temp;
                                cosi = max(min(cosi,1),-1);
                            end
                        end

                        % Final rotation of Stokes vector
                        sini  = real(sqrt(1-cosi*cosi));
                        cos22 = 2*cosi*cosi-1;
                        sin22 = 2*sini*cosi;

                        S2(1) = S(1);
                        S2(2) = S(2)*cos22 - S(3)*sin22;
                        S2(3) = S(2)*sin22 + S(3)*cos22;
                        S2(4) = S(4);

                        S = S2/S2(1);
                        U = U2;


                        % ROULETTE
                        if (W < THRESHOLD)
                            if (rand <= CHANCE)
                                W = W/CHANCE;
                            else
                                photon_status = DEAD;
                            end
                        end
                    end
                end
            end % end ALIVE loop

        end % end photon loop for this scan

    end % end scan-position loop

    %----------------- Save outputs (once per depth and input state) -----------------
    if jjj == 1
        outHI = IR; outHQ = QR; outHU = UR; outHV = VR;
        parsaveHI(fullfile(outputDir,'outHI.mat'), outHI);
        parsaveHQ(fullfile(outputDir,'outHQ.mat'), outHQ);
        parsaveHU(fullfile(outputDir,'outHU.mat'), outHU);
        parsaveHV(fullfile(outputDir,'outHV.mat'), outHV);
        parsaveHIDeep(fullfile(outputDir,'outHIDeep.mat'), IRDeep);
        parsaveHQDeep(fullfile(outputDir,'outHQDeep.mat'), QRDeep);
        parsaveHUDeep(fullfile(outputDir,'outHUDeep.mat'), URDeep);
        parsaveHVDeep(fullfile(outputDir,'outHVDeep.mat'), VRDeep);
      
    elseif jjj == 2
        outVI = IR; outVQ = QR; outVU = UR; outVV = VR;
        parsaveVI(fullfile(outputDir,'outVI.mat'), outVI);
        parsaveVQ(fullfile(outputDir,'outVQ.mat'), outVQ);
        parsaveVU(fullfile(outputDir,'outVU.mat'), outVU);
        parsaveVV(fullfile(outputDir,'outVV.mat'), outVV);
        parsaveVIDeep(fullfile(outputDir,'outVIDeep.mat'), IRDeep);
        parsaveVQDeep(fullfile(outputDir,'outVQDeep.mat'), QRDeep);
        parsaveVUDeep(fullfile(outputDir,'outVUDeep.mat'), URDeep);
        parsaveVVDeep(fullfile(outputDir,'outVVDeep.mat'), VRDeep);

    elseif jjj == 3
        outPI = IR; outPQ = QR; outPU = UR; outPV = VR;
        parsavePI(fullfile(outputDir,'outPI.mat'), outPI);
        parsavePQ(fullfile(outputDir,'outPQ.mat'), outPQ);
        parsavePU(fullfile(outputDir,'outPU.mat'), outPU);
        parsavePV(fullfile(outputDir,'outPV.mat'), outPV);
        parsavePIDeep(fullfile(outputDir,'outPIDeep.mat'), IRDeep);
        parsavePQDeep(fullfile(outputDir,'outPQDeep.mat'), QRDeep);
        parsavePUDeep(fullfile(outputDir,'outPUDeep.mat'), URDeep);
        parsavePVDeep(fullfile(outputDir,'outPVDeep.mat'), VRDeep);

    elseif jjj == 4
        outRI = IR; outRQ = QR; outRU = UR; outRV = VR;
        parsaveRI(fullfile(outputDir,'outRI.mat'), outRI);
        parsaveRQ(fullfile(outputDir,'outRQ.mat'), outRQ);
        parsaveRU(fullfile(outputDir,'outRU.mat'), outRU);
        parsaveRV(fullfile(outputDir,'outRV.mat'), outRV);
        parsaveRIDeep(fullfile(outputDir,'outRIDeep.mat'), IRDeep);
        parsaveRQDeep(fullfile(outputDir,'outRQDeep.mat'), QRDeep);
        parsaveRUDeep(fullfile(outputDir,'outRUDeep.mat'), URDeep);
        parsaveRVDeep(fullfile(outputDir,'outRVDeep.mat'), VRDeep);
    end

end % parfor

toc
disp('done')
toc
