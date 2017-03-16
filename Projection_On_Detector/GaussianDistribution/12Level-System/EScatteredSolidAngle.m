function E = EScatteredSolidAngle( xAtom, yAtom, zAtom, xScreen, yScreen, zScreen, LaserPhase, LaserMag )
%SPHERICALWAVE Calculates the field of a spherical wave on a round aperture

k = 2*pi;

r = sqrt(bsxfun(@minus,xScreen,permute(xAtom,[3 2 1])).^2 ...
       + bsxfun(@minus,yScreen,permute(yAtom,[3 2 1])).^2 ...
       + bsxfun(@minus,zScreen,permute(zAtom,[3 2 1])).^2);

phi = bsxfun(@plus,k.*r,LaserPhase);
E = bsxfun(@times,-3/4*exp(1i.*phi)./(k.*r),LaserMag);%.*sin(theta);

end