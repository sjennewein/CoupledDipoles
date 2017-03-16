clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                    %%%
%%%     Atoms                          %%%
%%%     Energies normalized to GAMMA   %%%
%%%     Length normalized to LAMBDA    %%%
%%%                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda=780E-9; %780nm
Gamma=1; % 6MHz linewidth of rubidium
k = 2*pi; %units of lambda

cloudX=200E-9 /lambda;
cloudY=cloudX;
cloudZ=1200E-9/lambda;

eDipole3 = zeros(1,1,3);
eDipole3(1,1,1) = 1; %x
eDipole3(1,1,2) = 0; %y
eDipole3(1,1,3) = 0; %z

V = @(r,theta) -3/4.*1./(k.*r).^3 .* exp(1i.*r.*k) ... % definition from paper observation of suppression
    .* ((1-3.*cos(theta).^2) ...
    .*  (1i.*k.*r - 1) ...
    +   sin(theta).^2.*(k.*r).^2);

%%%%%%%%%%%%%%%
%%%         %%%
%%%  Laser  %%%
%%%         %%%
%%%%%%%%%%%%%%%

s = 0.1;

waist = 1.2E-6/lambda; %focal spot
Omega = sqrt(s/2); %reduced omega

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
    - 1i .* pi/2 ...
    )./ w(z);

ELaserFarMag = @(rho,z) exp(-rho.^2./w(z).^2)./ w(z);
ELaserFarPhase = @(rho,r) k .* r - pi/2;

%%%%%%%%%%%%%%%%%%%%%
%%%               %%%
%%%    Detector   %%%
%%%               %%%
%%%%%%%%%%%%%%%%%%%%%

NA = 0.5; %numerical aperture of our lens
screenR = 100;

[phi,wp]=retgauss(0,2*pi,6,6);
[theta,wt]=retgauss(0,asin(NA),6,6);
[phig,thetag]=ndgrid(phi,theta);
[Wp,Wt]=ndgrid(wp,wt);
sinphi = sin(thetag);

screenZ   = screenR.*cos(thetag);
screenX   = screenR.*sin(thetag).*cos(phig);
screenY   = screenR.*sin(thetag).*sin(phig);
screenRho = screenR.*sin(thetag);
screenR1  = repmat(screenR, size(screenRho));

ELaserScreen      = ELaserFar(screenRho,screenR1);
ELaserScreenMag   = ELaserFarMag(screenRho,screenR1);
ELaserScreenPhase = ELaserFarPhase(screenRho,screenR1);
ILaserIntegr      = screenR^2.*sum(Wp(:).*Wt(:).*sinphi(:).*abs(ELaserScreen(:)).^2);

%%%%%%%%%%%%%%%%%%%%%%
%%%                %%%
%%%  Run Specific  %%%
%%%                %%%
%%%%%%%%%%%%%%%%%%%%%%

% nAtoms = [1 2 3 4 5 6 7 8 9 10:20:200];
nAtoms = [2 10]
detunings = -20:0.1:20;
realizations = 100;

TransmissionSpectra = cell(size(nAtoms));
SofOmegaSpectra     = cell(size(nAtoms));


%%%%%%%%%%%%%%%%%%%%%%%%
%%%                  %%%
%%%  Implementation  %%%
%%%                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%

for iAtom=1:numel(nAtoms)
    tic;
    cAtom = nAtoms(iAtom);
    disp(cAtom);
    resultSOmega       = zeros(size(detunings));
    resultTOmega       = zeros(size(detunings));
    for j=1:realizations
        atomX = cloudX*randn(cAtom,1);
        atomY = cloudY*randn(cAtom,1);
        atomZ = cloudZ*randn(cAtom,1);
        
        atomAtomX = repmat(atomX, 1, cAtom);
        atomAtomY = repmat(atomY, 1, cAtom);
        atomAtomZ = repmat(atomZ, 1, cAtom);
        
        atomAtomR = cat(3, atomAtomX-atomAtomX', ...
            atomAtomY-atomAtomY', atomAtomZ-atomAtomZ');
        
        atomAtomDist = sqrt(atomAtomR(:,:,1).^2 + ...
            atomAtomR(:,:,2).^2 + ...
            atomAtomR(:,:,3).^2);
        
        atomLaserRho = sqrt(atomX.^2 + atomY.^2);
        
        eDipole = repmat(eDipole3,cAtom);
        
        absDipole = sum(eDipole3 .* eDipole3,3);
        atomAtomTheta = acos(sum(eDipole.*atomAtomR,3)./(absDipole .* atomAtomDist));
        atomAtomTheta(logical(eye(cAtom))) = 0;
        
        interaction = V(atomAtomDist,atomAtomTheta);
        interaction(logical(eye(cAtom))) = 0;
        
        
        driving = ELaser(atomLaserRho,atomZ);
        
        
        EscElS = screenR^2.*EScatteredSolidAngleS(atomX,atomY,atomZ,screenX,screenY,screenZ,-ELaserScreenPhase,ELaserScreenMag,sinphi,Wp,Wt); %sphericalWave + phase of laser
        
        Esc = EScatteredSolidAngleT(atomX,atomY,atomZ,screenX,screenY,screenZ); %sphericalWave + phase of laser
        
        for i = 1:numel(detunings)
            iDetuning = detunings(i);
            hamiltonian = (eye(cAtom) * (iDetuning+(1i/2))) - interaction;
            dj = linsolve(hamiltonian,driving);
            sOmega = sum(dj.*EscElS);
            t = abs(sum(bsxfun(@times,Esc,permute(dj,[2 3 1])),3)+ELaserScreen).^2;
            tOmega = screenR^2.*sum(t(:).*sinphi(:).*Wt(:).*Wp(:));
            resultTOmega(i) = resultTOmega(i) + tOmega;
            resultSOmega(i) = resultSOmega(i) + sOmega;
        end
    end
    
    SofOmegaSpectra{iAtom} = resultSOmega;
    TransmissionSpectra{iAtom} = resultTOmega;
    toc;
end;
denominator = ILaserIntegr*realizations;

save('2lvlResponse.mat','detunings','nAtoms','denominator','waist','SofOmegaSpectra','TransmissionSpectra','cloudX','cloudY','cloudZ','screenR','NA');

figure('name','s omega');
hold on;
plot(detunings,abs(resultSOmega/denominator+1).^2,'x-')
ylim([0 1]);

figure('name','transmission');
plot(detunings,resultTOmega/denominator,'o-');
ylim([0 1]);
