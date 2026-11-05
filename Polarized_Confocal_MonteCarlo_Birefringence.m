% Polarized Monte Carlo simulation of a scanning confocal reflectance 
% Mueller-matrix microscope imaging a birefringent slab of Mie-scattering spheres. 
% Optics, transport, and detection are the same as in
% Polarized_Confocal_MonteCarlo_Spheres.m; this script adds a birefringent
% medium
%
% What the code does: 
% For each input polarization state (H,V, P, R) and each focal depth listed
% in depth_um_list, the simulation: 
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



clear all
close all

tic

Nphotons = 1e5; % photons launched per scan position

% ------- Sample -------------%
depth_um_list = [0 10 20 30 40 50]; % [um] focus depths
slabsize = 0.005; % [cm] slab thickness (50 um)
dn = 3.46e-4;  % material birefringence (ne - no)
% Fiber orientation
gamma = 10*pi/180;  % [rad] fiber angle in the x-y plane
CC = updateU([1 0 0], pi/2, gamma);
disp('Fiber axis CC:');
disp(CC);
no        = 1.393; % ordinary axis refractive index
nDepth        = numel(depth_um_list);

%  MIE SCATTERING PARAMETERS
radius = 0.15;  % [um] sphere radius
mus = 60;       % [cm^-1] scattering coefficient
lambda = 0.783; % [um] wavelength
mua     = 0.0;  % [cm^-1] absorption coefficient
nre_med = 1.0;  % medium refractive index
nre_p   = 1.46; % sphere refractive index
nangles = 1000; % angular samples
conv    = 1; 

% CALCULATE MIE ONCE
[S_mie, C, ang] = calcmie(radius, nre_p, nre_med, lambda, nangles, 'ConvergenceFactor', conv);
Q_eff   = getEfficiencies(C, radius, 3);
M_muel  = getMuellerMatrix(S_mie);
s11 = squeeze(M_muel(1,1,:));
s12 = squeeze(M_muel(1,2,:));
s33 = squeeze(M_muel(3,3,:));
s43 = squeeze(M_muel(4,3,:));
xxx    = 2*pi*radius/(lambda/nre_med);
aa     = mie(nre_p/nre_med, xxx);
g      = aa(8); % anisotropy
albedo = mus/(mus + mua);


%  CONFOCAL OPTICS PARAMETERS
obj_gap_ref = 3.4;              % [cm] focus at interface
f_obj      = 4;                 % [cm] 5x objective, EFL 40 mm
Z_pupil    = f_obj;             % [cm]
alpha_fill = 0.9;               % alpha fill factor (fill factor of the objective pupil
NA          = 0.14;             % numerical aperture
R_p      = f_obj * NA;          % [cm] pupil radius 
w_p      = alpha_fill * R_p;    % [cm] beam waist profile
Zp       = -Z_pupil;            % [cm] pupil plane z 
ttt1     = asin(NA / nre_med);  % [rad] largest angle accepted in the tissue (NA gate)
M_mag    = 5;                   % set magnification of system (5x)
rphys    = 0.03/2;              % [cm], 300 µm diameter
pinholeR = rphys / M_mag;       % [cm], back-projected pinhole radius
n_tissue = nre_med;
% Lateral constraints on sample
boxHalfSizeX = 0.5; % [cm] 
boxHalfSizeY = 0.5; % [cm]

THRESHOLD = 0.01; % Russian roulette weight threshold
CHANCE    = 0.1;  % Russian roulette survival probability
ALIVE     = 1;
DEAD      = 0;

% SCANNING PARAMETERS
XpCtr = 0;    % [cm] pupil center x
YpCtr = 0;    % [cm] pupil center y
% [X0, Y0] scan centers
step = 3e-4;      % [cm] step size: 3 µm in cm
vals = (-9.5:9.5) * step; % set number of scan points (20 x 20)
[X, Y] = meshgrid(vals, vals);
scan_pts = [X(:), Y(:)];
nScan = size(scan_pts, 1);

% create parallel pool
if isempty(gcp('nocreate'))
    parpool('local', 4);
end

outputDir = pwd;   % base output folder (pwd = current folder)
runTag = 'run_1'; % appended to every output file name

parfor idepth = 1:nDepth

    % calculate focal plane depth
    depth_um  = depth_um_list(idepth);
    Fdepth_cm = depth_um * 1e-4;                    % [cm] focal depth
    Foc       = obj_gap_ref - Fdepth_cm / nre_med;  % objective gap for this depth
    obj_gap   = Foc;
    Fdepth    = (nre_med/1) * (obj_gap_ref - obj_gap); % [cm] focal depth in the tissue
    fprintf('\nDEPTH %d/%d: %g um\n', idepth, nDepth, depth_um);
    depthDir = fullfile(outputDir, sprintf('depth%dum', depth_um));
    if ~exist(depthDir, 'dir')
        mkdir(depthDir);
    end

    % initialize polarization states
    for jjj = 1:4
        if jjj == 1
            S0 = [1 1 0 0];   % H
        elseif jjj == 2
            S0 = [1 -1 0 0 ]; % V
        elseif jjj == 3
            S0 = [1 0 1 0];  % P
        else
            S0 = [1 0 0 1];  % R
        end

        % extraordinary refractive index
        ne_loc   = no + dn;
       
        % detector  gridparameters
        hw = 0.015;   % [cm] half-width of the image grid
        NN = 100;     % lateral bins
        NNz = 250;    % no 3D bins in this script
        dx = 2*hw/NN; % [cm] lateral bin
        dy = 2*hw/NN; 
        dz = slabsize/NNz; % no 3D bins in this script

        % Local accumulators
        RR = zeros(1,4);   % total detected Stokes
        TR = zeros(1,4);   % transmitted Stokes
        IR = zeros(NN,NN); % detected Stokes images
        QR = zeros(NN,NN);
        UR = zeros(NN,NN);
        VR = zeros(NN,NN);
       

     
        % center coordinates
        X0=0;
        Y0=0;
        %-------------- Loop over scan positions ------------
        for s = 1:nScan

            % Scan center for this scan position
            X0 = scan_pts(s,1);   % [cm]
            Y0 = scan_pts(s,2);   % [cm]
           
            %-------------- Photon loop --------------
            for i_photon = 1:Nphotons

                % --- Launch photon from pupil ---
                u_r   = rand;
                trunc = exp(-2*R_p^2/w_p^2);
                r     = (w_p/sqrt(2)) * sqrt(-log(1 - u_r*(1 - trunc)));
                psi = 2*pi*rand;
                x_p = XpCtr + r*cos(psi);
                y_p = YpCtr + r*sin(psi);
                a  = (r / (f_obj * n_tissue))^2;
                a  = min(a, 1 - 1e-12);
                Hr = Fdepth * sqrt(a / (1 - a));
                Hx = X0 + Hr*cos(psi);
                Hy = Y0 + Hr*sin(psi);
                Hz = 0;
                v_air = [Hx-x_p, Hy-y_p, Hz-Zp];
                v_air = v_air / norm(v_air);
                U = v_air; % no Fresnel

                % meridian plane rotation
                S  = S0;
                S2 = rotSphi(S, psi);
                S  = [1, S2(2)/S2(1), S2(3)/S2(1), S2(4)/S2(1)];

                x = Hx; y = Hy; z = Hz;
                photon_status = ALIVE;
                W = 1;

                %******** ALIVE cycle *****************/
                while photon_status == ALIVE

                    xlast = x; ylast = y; zlast = z;

                    s_step = -log(rand)/(mus + mua);  % [cm] step length
                    x = x + U(1)*s_step;
                    y = y + U(2)*s_step;
                    z = z + U(3)*s_step;

                    % calculate how much of the step above was in the tissue
                    s_in = s_step;
                    if z <= 0
                        s_in = s_step * zlast / (zlast - z); % crossed top surface
                    elseif z >= slabsize
                        s_in = s_step * (slabsize - zlast) / ( z - zlast); % crossed bottom surface
                    end

                    %---------------- Birefringence ----------------
                    phi=atan2(U(2),U(1)); % position of meridian plane 

                    % Fiber axis angle in photon meridian frame. Frame of upward traveling photon is mirrored, so sign is flipped. 
                    % Only for fibers in x-y plane
                    if U(3) >= 0
                        beta = gamma - phi; % going into the tissue (down)
                    else
                        beta = phi - gamma; % coming back out (up)
                    end

                    S2= rotSphi(S, beta); % stokes vector in the fiber axis frame

                    %Now calculate the incident angle of the U to fiber*/
                    if  dot(CC,U)==0
                        csi=pi/2;
                    else
                        csi = atan2(norm(cross(CC,U)),dot(CC,U));
                    end

                    nsmax = no*ne_loc/sqrt(ne_loc*ne_loc*cos(csi)*cos(csi)+no*no*sin(csi)*sin(csi));%/*maximal refractive index seen by photon traveling in thetax direction with respect to fast axis*/
                    delta0 = (nsmax-no)*2*pi*(s_in*1e4)/(lambda); %/* [rad] retardance ; birefringence=nsmax-nf */

                    % linear retarder
                    S(1) =  S2(1);
                    S(2) =  S2(2);
                    S(3) =  S2(3)*cos(delta0) + S2(4)*sin(delta0);
                    S(4) = -S2(3)*sin(delta0) + S2(4)*cos(delta0);

                    % back to the photon's meridian frame
                    S2   = rotSphi(S, -beta);
                    S    = S2 / S2(1);

                    % Box constraint
                    if abs(x) > boxHalfSizeX || abs(y) > boxHalfSizeY
                        photon_status = DEAD;
                        continue;
                    end

                    % ABSORB
                    W = W * albedo;

                    % bottom surface: photon leaves slab and is terminated
                    % here 
                    if z>=slabsize
                        photon_status = DEAD;
                        continue;
                    end

                    %-------- Top surface (reflectance) --------
                    if z <= 0

                       % back-project the last segment to the focal plane z = Fdepth
                        Xfoc = x - (z - Fdepth)*((x-xlast)/(z-zlast));
                        Yfoc = y - (z - Fdepth)*((y-ylast)/(z-zlast));
                        % relative to scan position in global coordinates
                        Xrel = Xfoc - X0;
                        Yrel = Yfoc - Y0;

                        % NA angle check
                        ttt = acos(abs(U(3)));
                        if ttt < ttt1
                            w1 = sqrt(Xrel^2 + Yrel^2);

                            % pinhole check - DETECTED
                            if w1 <= pinholeR
                                phi = atan2(U(2), U(1));
                                S2  = rotSphi(S, phi);
                                RR  = RR + S2*W;
                              
                                % within detector
                                % image binning, only inside the +/- hw grid
                                if abs(Xfoc) <= hw && abs(Yfoc) <= hw
                                    ix = floor((Xfoc + hw)/dx) + 1;
                                    iy = floor((Yfoc + hw)/dy) + 1;

                                    if ix>=1 && ix<=NN && iy>=1 && iy<=NN


                                        if jjj == 1

                                            IR(iy,ix) = IR(iy,ix) + W*S2(1);
                                            QR(iy,ix) = QR(iy,ix) + W*S2(2);
                                            UR(iy,ix) = UR(iy,ix) + W*S2(3);
                                            VR(iy,ix) = VR(iy,ix) + W*S2(4);



                                        elseif jjj == 2

                                            IR(iy,ix) = IR(iy,ix) + W*S2(1);
                                            QR(iy,ix) = QR(iy,ix) + W*S2(2);
                                            UR(iy,ix) = UR(iy,ix) + W*S2(3);
                                            VR(iy,ix) = VR(iy,ix) + W*S2(4);



                                        elseif jjj == 3

                                            IR(iy,ix) = IR(iy,ix) + W*S2(1);
                                            QR(iy,ix) = QR(iy,ix) + W*S2(2);
                                            UR(iy,ix) = UR(iy,ix) + W*S2(3);
                                            VR(iy,ix) = VR(iy,ix) + W*S2(4);

                                        else

                                            IR(iy,ix) = IR(iy,ix) + W*S2(1);
                                            QR(iy,ix) = QR(iy,ix) + W*S2(2);
                                            UR(iy,ix) = UR(iy,ix) + W*S2(3);
                                            VR(iy,ix) = VR(iy,ix) + W*S2(4);



                                        end

                                    end
                                end

                            end
                        end
                        photon_status = DEAD;
                        continue;
                    end


                    %-------- Bottom surface (transmission) --------
                    if z >= slabsize
                        phi = -atan2(U(2), U(1));
                        S2  = rotSphi(S, phi);                       
                        photon_status = DEAD;
                        continue;
                    end

                    %-------- SPIN --------
                    condition = true;
                    while condition
                        theta   = acos(2*rand - 1);
                        phi     = rand * 2*pi;
                        ithedeg = 1 + floor(theta*nangles/pi);
                        if ithedeg > nangles; ithedeg = nangles; end
                        I0 = s11(1)*S(1) + s12(1)*(S(2)*cos(2*phi) + S(3)*sin(2*phi));
                        I  = s11(ithedeg)*S(1) + s12(ithedeg)*(S(2)*cos(2*phi) + S(3)*sin(2*phi));
                        if rand*I0 < I
                            condition = false;
                        end
                    end

                    U2       = updateU(U, phi, theta);
                    costheta = cos(theta);
                    S2       = rotSphi(S, phi);

                    S(1) =  s11(ithedeg)*S2(1) + s12(ithedeg)*S2(2);
                    S(2) =  s12(ithedeg)*S2(1) + s11(ithedeg)*S2(2);
                    S(3) =  s33(ithedeg)*S2(3) + s43(ithedeg)*S2(4);
                    S(4) = -s43(ithedeg)*S2(3) + s33(ithedeg)*S2(4);

                    temp = sqrt(1-costheta^2) * sqrt(1-U2(3)^2);
                    if temp == 0
                        cosi = 0;
                    else
                        if (phi > pi) && (phi < 2*pi)
                            cosi = (U2(3)*costheta - U(3)) / temp;
                        else
                            cosi = -(U2(3)*costheta - U(3)) / temp;
                            cosi = max(min(cosi,1),-1);
                        end
                    end

                    sini  = real(sqrt(1 - cosi^2));
                    cos22 = 2*cosi^2 - 1;
                    sin22 = 2*sini*cosi;

                    S2(1) = S(1);
                    S2(2) = S(2)*cos22 - S(3)*sin22;
                    S2(3) = S(2)*sin22 + S(3)*cos22;
                    S2(4) = S(4);
                    S = S2 / S2(1);
                    U = U2;



                    %-------- ROULETTE --------
                    if W < THRESHOLD
                        if rand <= CHANCE
                            W = W / CHANCE;
                        else
                            photon_status = DEAD;
                        end
                    end

                end % while ALIVE

            end % photon loop
        end % scan loop



             %----------------- Save outputs -----------------
        % All names use runTag (set at the top of the script)
        if jjj == 1
            outHI = IR; outHQ = QR; outHU = UR; outHV = VR;
            parsaveHI(fullfile(depthDir, ['outHI_' runTag '.mat']), outHI);    
            parsaveHQ(fullfile(depthDir, ['outHQ_' runTag '.mat']), outHQ);
            parsaveHU(fullfile(depthDir, ['outHU_' runTag '.mat']), outHU);
            parsaveHV(fullfile(depthDir, ['outHV_' runTag '.mat']), outHV);
        elseif jjj == 2
            outVI = IR; outVQ = QR; outVU = UR; outVV = VR;
            parsaveVI(fullfile(depthDir, ['outVI_' runTag '.mat']), outVI);
            parsaveVQ(fullfile(depthDir, ['outVQ_' runTag '.mat']), outVQ);
            parsaveVU(fullfile(depthDir, ['outVU_' runTag '.mat']), outVU);
            parsaveVV(fullfile(depthDir, ['outVV_' runTag '.mat']), outVV);
        elseif jjj == 3
            outPI = IR; outPQ = QR; outPU = UR; outPV = VR;
            parsavePI(fullfile(depthDir, ['outPI_' runTag '.mat']), outPI);
            parsavePQ(fullfile(depthDir, ['outPQ_' runTag '.mat']), outPQ);
            parsavePU(fullfile(depthDir, ['outPU_' runTag '.mat']), outPU);
            parsavePV(fullfile(depthDir, ['outPV_' runTag '.mat']), outPV);
        elseif jjj == 4
            outRI = IR; outRQ = QR; outRU = UR; outRV = VR;
            parsaveRI(fullfile(depthDir, ['outRI_' runTag '.mat']), outRI);
            parsaveRQ(fullfile(depthDir, ['outRQ_' runTag '.mat']), outRQ);
            parsaveRU(fullfile(depthDir, ['outRU_' runTag '.mat']), outRU);
            parsaveRV(fullfile(depthDir, ['outRV_' runTag '.mat']), outRV);
        end

        %RR
    end % polarization loop
end % depth loop
toc
