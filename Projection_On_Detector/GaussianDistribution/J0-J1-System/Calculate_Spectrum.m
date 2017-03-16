clear;
close all;
profile off
profile clear
% profile on;

tic;
%%%%%%%%%%%%%%%%%%%%%
%%%               %%%
%%%     Atoms     %%%
%%%               %%%
%%%%%%%%%%%%%%%%%%%%%

lambda=780E-9; %780nm 
Gamma=1; % 6MHz linewidth of rubidium
k = 2*pi; %units of lambda

cloudX=200E-9 /lambda;
cloudY=cloudX;
cloudZ=1200E-9/lambda;

dipoleAlpha = zeros(3,3,3);
dipoleAlpha(1,:,1) = ones(1,3);
dipoleAlpha(2,:,2) = ones(1,3);
dipoleAlpha(3,:,3) = ones(1,3);

%%%%%%%%%%%%%%%
%%%         %%%
%%%  Laser  %%%
%%%         %%%
%%%%%%%%%%%%%%%

s = 0.1;
polarization = [1 0 0]';
waist = 1.2E-6/lambda; %focal spot of our laser 0.5um
Omega = sqrt(s*Gamma^2/2);

zR = pi * waist^2; %rayleigh range

w = @(z) waist .* sqrt(1+(z ./ zR).^2);
G = @(z) atan(z ./ zR);
R = @(z) z .* (1 + (zR ./ z).^2);

ELaser = @(rho,z) exp(-rho.^2./w(z).^2 ...
    + 1i .* k .* z ...
    + 1i .* k .* rho.^2 ./ (2 .* R(z)) ...
    - 1i .* G(z) ...
    )./ w(z);

ELaserFar = @(rho,z) exp(-rho.^2./w(z).^2 ...
    + 1i .* k .* z ...
    + 1i .* k .* rho.^2 ./ (2 .* R(z)) ...
    - 1i .* pi/2 ...
    )./ w(z);

ELaserFarMag = @(rho,z) exp(-rho.^2./w(z).^2)./ w(z);
ELaserFarPhase = @(rho,r) k .* r - pi/2;

%%%%%%%%%%%%%%%%%%%%%
%%%               %%%
%%%    Detector   %%%
%%%               %%%
%%%%%%%%%%%%%%%%%%%%%
NA = 0.5;

screenR = 100;
[phi,wp]=retgauss(0,2*pi,6,6);
[theta,wt]=retgauss(0,asin(NA),6,6);
[phig,thetag]=ndgrid(phi,theta);
[Wp,Wt]=ndgrid(wp,wt);
sintheta = sin(thetag);


screenZ   = screenR.*cos(thetag);
screenX   = screenR.*sin(thetag).*cos(phig);
screenY   = screenR.*sin(thetag).*sin(phig);
screenRho = screenR.*sin(thetag);
screenR1  = repmat(screenR, size(screenRho));


ELaserScreen      = ELaserFar(screenRho,screenR1);
ELaserScreenMag   = ELaserFarMag(screenRho,screenR1);
ELaserScreenPhase = ELaserFarPhase(screenRho,screenR1);
ILaserIntegr      = screenR^2.*sum(Wp(:).*Wt(:).*sintheta(:).*abs(ELaserScreen(:)).^2);

%%%%%%%%%%%%%%%%%%%%%%
%%%                %%%
%%%  Run Specific  %%%
%%%                %%%
%%%%%%%%%%%%%%%%%%%%%%

% nAtoms = [1 2 3 4 5 6 7 8 9 10:10:200];
nAtoms = [10];
speed = zeros(1,5);
detunings = -20:0.1:20;
realizations = 100;


%%%%%%%%%%%%%%%%%%%%%%%%
%%%                  %%%
%%%  Implementation  %%%
%%%                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%
count = 0;

SofOmegaSpectra     = cell(size(nAtoms));


for iAtom = 1:numel(nAtoms)
    tic;    
    cAtom = nAtoms(iAtom);
    disp(cAtom);
    resultSOmega       = zeros(size(detunings));
    resultSOmegaYZ       = zeros(size(detunings));

    fullDipoleAlpha = cat(3, kron(ones(cAtom)-eye(cAtom),dipoleAlpha(:,:,1)), ...
                             kron(ones(cAtom)-eye(cAtom),dipoleAlpha(:,:,2)), ...
                             kron(ones(cAtom)-eye(cAtom),dipoleAlpha(:,:,3)));
    
    dipoleBeta = permute(dipoleAlpha,[2 3 1]);
    
    fullDipoleBeta = cat(3, kron(ones(cAtom)-eye(cAtom),dipoleBeta(:,:,1)), ...
                            kron(ones(cAtom)-eye(cAtom),dipoleBeta(:,:,2)), ...
                            kron(ones(cAtom)-eye(cAtom),dipoleBeta(:,:,3)));
        
    vec3NPolarization = repmat(polarization,cAtom,1);

    for iRealization = 1:realizations
        tic;
        atomX = cloudX*randn(cAtom,1);
        atomY = cloudY*randn(cAtom,1);
        atomZ = cloudZ*randn(cAtom,1);
        
        atomAtomX = repmat(atomX, 1, cAtom);
        atomAtomY = repmat(atomY, 1, cAtom);
        atomAtomZ = repmat(atomZ, 1, cAtom);
        
        atomAtomR = cat(3, atomAtomX-atomAtomX', ...
                           atomAtomY-atomAtomY', ...
                           atomAtomZ-atomAtomZ');
        
        fullAtomR = cat(3, kron(atomAtomR(:,:,1),ones(3)), ...
                           kron(atomAtomR(:,:,2),ones(3)), ...
                           kron(atomAtomR(:,:,3),ones(3)));
        
        atomAtomDist = sqrt(atomAtomR(:,:,1).^2 + ...
                            atomAtomR(:,:,2).^2 + ...
                            atomAtomR(:,:,3).^2);
        
        fullAtomDist = kron(atomAtomDist,ones(3));
        
        kronecker = kron(ones(cAtom),eye(3));
        
        q = kronecker ...
            -sum(fullDipoleAlpha.*fullAtomR,3) ...
            .*sum(fullDipoleBeta.*fullAtomR,3)./fullAtomDist.^2;
        q(isnan(q)) = 0;
        
        p = kronecker ...
            -3.*sum(fullDipoleAlpha.*fullAtomR,3) ...
            .*sum(fullDipoleBeta.*fullAtomR,3)./fullAtomDist.^2;
        p(isnan(p)) = 0;
        
        V = -3.*Gamma./(4.*(k.*fullAtomDist).^3) ...
            .*(p.*(1i.*k.*fullAtomDist-1) ...
            +q.*(k.*fullAtomDist).^2) ...
            .*exp(1i.*k.*fullAtomDist);
        
        V(isnan(V)) = 0;
        
        vec3NAtomX = kron(atomX,ones(3,1));
        vec3NAtomY = kron(atomY,ones(3,1));
        vec3NAtomZ = kron(atomZ,ones(3,1));
        
        atomAxisRho = sqrt(vec3NAtomX.^2 + vec3NAtomY.^2);
        
        driving = ELaser(atomAxisRho,vec3NAtomZ).*vec3NPolarization;

        EscElX = screenR^2.*EScatteredSolidAngle(dipoleAlpha(1,1,:),atomX,atomY,atomZ,screenX,screenY,screenZ,-ELaserScreenPhase,ELaserScreenMag,polarization,sintheta,Wp,Wt);
        EscElY = screenR^2.*EScatteredSolidAngle(dipoleAlpha(2,2,:),atomX,atomY,atomZ,screenX,screenY,screenZ,-ELaserScreenPhase,ELaserScreenMag,polarization,sintheta,Wp,Wt);
        EscElZ = screenR^2.*EScatteredSolidAngle(dipoleAlpha(3,3,:),atomX,atomY,atomZ,screenX,screenY,screenZ,-ELaserScreenPhase,ELaserScreenMag,polarization,sintheta,Wp,Wt);
        
        for iDetuning = 1:numel(detunings)
            cDetuning = detunings(iDetuning);
            hamiltonian = (eye(size(V)) * (cDetuning+(1i*Gamma/2))) - V;           
            dj = linsolve(hamiltonian,driving);
            dx = dj(1:3:end);
            dy = dj(2:3:end);
            dz = dj(3:3:end);
            somegaX = sum(EscElX.*dx);
            somegaY = sum(EscElY.*dy);
            somegaZ = sum(EscElZ.*dz);

            resultSOmega(iDetuning) = resultSOmega(iDetuning) + somegaX;
            resultSOmegaYZ(iDetuning) = resultSOmegaYZ(iDetuning) + somegaY + somegaZ;
        end

    end
        SofOmegaSpectra{iAtom} = resultSOmega;

   toc;
end


denominator = ILaserIntegr*realizations;
save('Spectra.mat', ...
     'SofOmegaSpectra','realizations', ...
     'detunings','nAtoms','denominator','cloudX','cloudY','cloudZ','waist');


figure('name','S(omega)');
hold on;
plot(detunings,abs(resultSOmega/denominator+1).^2,'x-');

