function parsaveTrack(fname, trk_pos, trk_U, trk_S_pre, trk_S_rot, ...
    trk_S_scat, trk_S_post, trk_theta, trk_phi, trk_W, trk_launch, trk_S0)
%PARSAVETRACK Save single-photon Stokes-evolution history from parfor.
%   Mirrors the parsave* pattern used elsewhere in the MC pipeline so the
%   save call works inside a parfor loop. Stores the full per-scatter
%   history of one detected photon for verification of polarization physics.
%
%   Fields saved:
%     trk_pos    [N x 3]  position (x,y,z) at each scattering event (pre-scatter)
%     trk_U      [N x 3]  direction cosines at each scattering event (pre-scatter)
%     trk_S_pre  [N x 4]  Stokes vector in lab frame BEFORE rotSphi
%     trk_S_rot  [N x 4]  Stokes AFTER rotation into scattering plane
%     trk_S_scat [N x 4]  Stokes AFTER Mie scattering matrix is applied
%     trk_S_post [N x 4]  Stokes AFTER final meridian rotation + normalize
%     trk_theta  [N x 1]  sampled polar scattering angle (rad)
%     trk_phi    [N x 1]  sampled azimuthal scattering angle (rad)
%     trk_W      [N x 1]  photon weight at each event
%     trk_launch struct   launch state (pos, U, S0, X0, Y0)
%     trk_S0     [1 x 4]  initial Stokes for this polarization channel
%
%   See match_confocal3JRR_AS_V2_TR_tracked.m for usage.

save(fname, 'trk_pos', 'trk_U', 'trk_S_pre', 'trk_S_rot', ...
    'trk_S_scat', 'trk_S_post', 'trk_theta', 'trk_phi', 'trk_W', ...
    'trk_launch', 'trk_S0');
end