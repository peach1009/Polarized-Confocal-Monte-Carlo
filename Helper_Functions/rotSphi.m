% Rotates a Stokes vector around the optical beam by an azimuthal angle,
% aligning it with the meridian plane defined by photon's new direction of
% propagation
function S2=rotSphi(S, phi)

	cos2phi = cos(2*phi);
	sin2phi = sin(2*phi); 
	
	S2(1) = S(1); 
	S2(2) = S(2)*cos2phi+S(3)*sin2phi; 
	S2(3) = -S(2)*sin2phi+S(3)*cos2phi;  
	S2(4) = S(4); 

   
% Rotation matrix applied to stokes vector according to the 
% [I';Q';U';V'] = [1, 0, 0, 0