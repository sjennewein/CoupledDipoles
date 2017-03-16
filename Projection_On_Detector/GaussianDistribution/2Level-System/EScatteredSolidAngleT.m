function Esc = EScatteredSolidAngleT( xAtom, yAtom, zAtom, xScreen, yScreen, zScreen)
%SPHERICALWAVE Calculates the field of a spherical wave on a round aperture

% |El|^2+|Esc|^2+EscEl*+Esc*El

%t1 = Esc;
%t2 = EscEl*;
%t3 = Esc*El;


k = 2*pi;
r = sqrt(bsxfun(@minus,xScreen,permute(xAtom,[3 2 1])).^2 ...
       + bsxfun(@minus,yScreen,permute(yAtom,[3 2 1])).^2 ...
       + bsxfun(@minus,zScreen,permute(zAtom,[3 2 1])).^2);

% theta = zeros(size(r));
% for iAtom = 1:length(xAtom)
%     cxAtom = xAtom(iAtom);
%     cyAtom = yAtom(iAtom);
%     czAtom = zAtom(iAtom);
%     atomScreenR = cat(3,xScreen-cxAtom, ...
%                       yScreen-cyAtom, ...
%                       zScreen-czAtom);
%     atomScreenDist = sqrt(atomScreenR(:,:,1).^2 + ...
%                           atomScreenR(:,:,2).^2 + ...
%                           atomScreenR(:,:,3).^2);
%     eDipole = repmat(dipole,size(xScreen));
%     absDipole = sum(dipole .* dipole,3);
%     atomScreenTheta = acos(sum(eDipole.*atomScreenR,3)./(absDipole .* atomScreenDist));
%     theta(:,:,iAtom) = atomScreenTheta;
% end

Esc = -3/4.*exp(1i.*k.*r)./(k.*r);

% trans = bsxfun(@times,-3/4.*exp(1i.*phi)./(k.*r),LaserMag);%.*sin(theta);
% EIntegrated = permute(sum(sum(bsxfun(@times,sinphi.*Wp.*Wt,trans))),[3 2 1]);

end