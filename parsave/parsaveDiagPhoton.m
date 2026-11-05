function parsaveDiagPhoton(filename, lastPath, lastMaxZ, lastFirstZ, ...
                                     lastExitAngle, lastFocalDist, lastBigScatters)
% Save the per-photon diagnostic records for detected photons.
%
%   lastPath        : [N x 1]  total path length traveled inside the tissue [cm]
%   lastMaxZ        : [N x 1]  maximum depth reached during the photon's history [cm]
%   lastFirstZ      : [N x 1]  depth of the first scattering event [cm]
%                              (NaN if the photon was never scattered before detection)
%   lastExitAngle   : [N x 1]  angle of the exit direction from the surface normal [rad]
%   lastFocalDist   : [N x 1]  3D distance from focal point (X0, Y0, Fdepth) to
%                              the last-scatter location [cm]
%   lastBigScatters : [N x 1]  number of large-angle (theta > pi/2) scatters this
%                              photon experienced
%
% These are aligned with lastPos / lastScat / lastStokes from outHPhotonData.mat etc.
% Index k in any of these refers to the same detected photon across all arrays.

save(filename, 'lastPath', 'lastMaxZ', 'lastFirstZ', ...
               'lastExitAngle', 'lastFocalDist', 'lastBigScatters', '-v7.3');

end