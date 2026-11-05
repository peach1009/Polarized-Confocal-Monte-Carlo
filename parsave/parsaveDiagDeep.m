function parsaveDiagDeep(filename, NDeep, W2Deep, Q2Deep, U2Deep, V2Deep)
% Save the voxel-level diagnostic arrays for noise-floor / variance analysis.
%
%   NDeep   : uint32 [NN x NN x NNz]  raw integer photon counts per voxel
%   W2Deep  : double [NN x NN x NNz]  sum of W^2 per voxel (for N_eff)
%   Q2Deep  : double [NN x NN x NNz]  sum of W*Q^2 per voxel (for Q variance)
%   U2Deep  : double [NN x NN x NNz]  sum of W*U^2 per voxel
%   V2Deep  : double [NN x NN x NNz]  sum of W*V^2 per voxel
%
% These can be combined with the existing *Deep arrays to compute:
%   N_eff(v)    = IDeep(v)^2 / W2Deep(v)              % effective # samples
%   Var(Q,v)    = Q2Deep(v)/IDeep(v) - (QDeep(v)/IDeep(v))^2
%   DOP_noise   = 1./sqrt(N_eff)                       % expected DOP from noise
%   DOP_signif  = DOP_voxel ./ DOP_noise               % real if >> 1

save(filename, 'NDeep', 'W2Deep', 'Q2Deep', 'U2Deep', 'V2Deep', '-v7.3');

end