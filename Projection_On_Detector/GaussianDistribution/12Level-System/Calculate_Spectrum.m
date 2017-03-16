clear;
close all;

tic;
%%%%%%%%%%%%%%%%%%%%%
%%%               %%%
%%%     Atoms     %%%
%%%               %%%
%%%%%%%%%%%%%%%%%%%%%
%%
lambda=780E-9; %780nm
Gamma=1; % 6MHz linewidth of rubidium
k = 2*pi; %units of lambda

cloudX=200E-9 /lambda;
cloudY=cloudX;
cloudZ=1200E-9/lambda;

%%
dipoleAlpha = zeros(3,3,3);
dipoleAlpha(1,:,1) = ones(1,3);
dipoleAlpha(2,:,2) = ones(1,3);
dipoleAlpha(3,:,3) = ones(1,3);
%%
ClebschAlpha = zeros(3,3,3,5);
ClebschAlpha(1,:,1,1) = sqrt(1/3)*ones(1,3);%transition pi, m=-2
ClebschAlpha(1,:,1,2) = sqrt(8/15)*ones(1,3);%transition pi, m=-1
ClebschAlpha(1,:,1,3) = sqrt(3/5)*ones(1,3);%transition pi, m=0
ClebschAlpha(1,:,1,4) = sqrt(8/15)*ones(1,3);%transition pi, m=1
ClebschAlpha(1,:,1,5) = sqrt(1/3)*ones(1,3);%transition pi, m=2

ClebschAlpha(2,:,2,1) = 1*ones(1,3);%transition sigma-, m=-2
ClebschAlpha(2,:,2,2) = sqrt(2/3)*ones(1,3);%transition sigma-, m=-1
ClebschAlpha(2,:,2,3) = sqrt(2/5)*ones(1,3);%transition sigma-, m=0
ClebschAlpha(2,:,2,4) = sqrt(1/5)*ones(1,3);%transition sigma-, m=1
ClebschAlpha(2,:,2,5) = sqrt(1/15)*ones(1,3);%transition sigma-, m=2

ClebschAlpha(3,:,3,1) = sqrt(1/15)*ones(1,3);%transition sigma+, m=-2
ClebschAlpha(3,:,3,2) = sqrt(1/5)*ones(1,3);%transition sigma+, m=-1
ClebschAlpha(3,:,3,3) = sqrt(2/5)*ones(1,3);%transition sigma+, m=0
ClebschAlpha(3,:,3,4) = sqrt(2/3)*ones(1,3);%transition sigma+, m=1
ClebschAlpha(3,:,3,5) = 1*ones(1,3);%transition sigma+, m=2

ClebschBeta = permute(ClebschAlpha, [2 1 3 4]);
%%
%%%%%%%%%%%%%%%
%%%         %%%
%%%  Laser  %%%
%%%         %%%
%%%%%%%%%%%%%%%

s = 0.1;
polarization = [1 0 0]';%PI
%polarization = [0 1 0]';%SIGMAM
%polarization = [0 0 1]';%SIGMAP

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

screenR = 50;
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


nAtoms = [50];
detunings = -5:0.1:5;
realizations = 100;


%%%%%%%%%%%%%%%%%%%%%%%%
%%%                  %%%
%%%  Implementation  %%%
%%%                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%
count = 0;

SofOmegaSpectra     = cell(size(nAtoms));


for iAtom = 1:numel(nAtoms)
    
    cAtom = nAtoms(iAtom);
    disp(cAtom);
    resultSOmega       = zeros(size(detunings));
    %%
    fullDipoleAlpha = cat(3, kron(ones(cAtom)-eye(cAtom),dipoleAlpha(:,:,1)), ...
        kron(ones(cAtom)-eye(cAtom),dipoleAlpha(:,:,2)), ...
        kron(ones(cAtom)-eye(cAtom),dipoleAlpha(:,:,3)));
    
    dipoleBeta = permute(dipoleAlpha,[2 3 1]);
    
    fullDipoleBeta = cat(3, kron(ones(cAtom)-eye(cAtom),dipoleBeta(:,:,1)), ...
        kron(ones(cAtom)-eye(cAtom),dipoleBeta(:,:,2)), ...
        kron(ones(cAtom)-eye(cAtom),dipoleBeta(:,:,3)));
    
    vec3NPolarization = repmat(polarization,cAtom,1);
    
    for iRealization = 1:realizations
        
        
        atomX = cloudX*randn(cAtom,1);
        atomY = cloudY*randn(cAtom,1);
        atomZ = cloudZ*randn(cAtom,1);
        atomM= randi([1 5],cAtom,1);
        atomPhase=0*2*pi*rand(cAtom,1);
        
        atomAtomX = repmat(atomX, 1, cAtom);
        atomAtomY = repmat(atomY, 1, cAtom);
        atomAtomZ = repmat(atomZ, 1, cAtom);
        atomAtomPhase = repmat(atomPhase, 1, cAtom);
        
        %%
        
        atomClebschAlpha=ClebschAlpha(:,:,:,atomM);
        B=permute(atomClebschAlpha,[1 4 2 3]);
        C=reshape(B,cAtom*3,3,3);
        fullClebschAlpha=repmat(C,1,cAtom);
        
        F = sum(C,3);
        atomClebschAlphaVect = F(:,1);
        
        atomClebschBeta=ClebschBeta(:,:,:,atomM);
        D=permute(atomClebschBeta,[1 4 2 3]);
        E=reshape(D,cAtom*3,3,3);
        fullClebschBeta=repmat(E,1,cAtom);
        %%
        
        atomAtomR = cat(3, atomAtomX-atomAtomX', ...
            atomAtomY-atomAtomY', ...
            atomAtomZ-atomAtomZ');
        
        atomAtomPhaseDiff = atomAtomPhase-atomAtomPhase';
        
        fullAtomR = cat(3, kron(atomAtomR(:,:,1),ones(3)), ...
            kron(atomAtomR(:,:,2),ones(3)), ...
            kron(atomAtomR(:,:,3),ones(3)));
        
        atomAtomDist = sqrt(atomAtomR(:,:,1).^2 + ...
            atomAtomR(:,:,2).^2 + ...
            atomAtomR(:,:,3).^2);
        
        fullAtomDist = kron(atomAtomDist,ones(3));
        
        fullAtomPhaseDiff = kron(atomAtomPhaseDiff,ones(3));
        
        kronecker = kron(ones(cAtom),eye(3));
        
        q = kronecker ...
            -sum(fullDipoleAlpha.*fullAtomR,3) ...
            .*sum(fullDipoleBeta.*fullAtomR,3)./fullAtomDist.^2;
        q(isnan(q)) = 0;
        
        p = kronecker ...
            -3.*sum(fullDipoleAlpha.*fullAtomR,3) ...
            .*sum(fullDipoleBeta.*fullAtomR,3)./fullAtomDist.^2;
        p(isnan(p)) = 0;
        %%
        V = -3.*Gamma./(4.*(k.*fullAtomDist).^3) ...
            .*sum(fullClebschAlpha,3) ...
            .*sum(fullClebschBeta,3) ...
            .*(p.*(1i.*k.*fullAtomDist-1) ...
            +q.*(k.*fullAtomDist).^2) ...
            .*exp(1i.*k.*fullAtomDist) ...
            .*exp(1i.*fullAtomPhaseDiff);
        V(isnan(V)) = 0;
        
        vec3NAtomX = kron(atomX,ones(3,1));
        vec3NAtomY = kron(atomY,ones(3,1));
        vec3NAtomZ = kron(atomZ,ones(3,1));
        
        atomAxisRho = sqrt(vec3NAtomX.^2 + vec3NAtomY.^2);
        
        driving = ELaser(atomAxisRho,vec3NAtomZ).*vec3NPolarization.*atomClebschAlphaVect;
        
        sphericalWavesDipolePI = EScatteredSolidAngle(atomX,atomY,atomZ,screenX,screenY,screenZ,-ELaserScreenPhase,ELaserScreenMag);
        
        for iDetuning = 1:numel(detunings)
            cDetuning = detunings(iDetuning);
            hamiltonian = (eye(size(V)) * (cDetuning+(1i*Gamma/2))) - V;
            
            dj = linsolve(hamiltonian,driving);
            dPI = dj(1:3:end).*atomClebschAlphaVect(1:3:end);
            EscElPI = sum(bsxfun(@times, sphericalWavesDipolePI,permute(dPI,[3 2 1])),3);
            
            somega = screenR^2.*sum(sintheta(:).*Wp(:).*Wt(:).*EscElPI(:));
            
            resultSOmega(iDetuning) = resultSOmega(iDetuning) + somega;
        end
        
    end
    SofOmegaSpectra{iAtom} = resultSOmega;
    
    
end
toc;
denominator = ILaserIntegr*realizations;
SofOmegaSquaredSpectra = abs(resultSOmega/denominator+1).^2;

%%
save('12lvlSpectra_1Atom_mm2_mm1_m0_m1_m2.mat', 'SofOmegaSpectra','SofOmegaSquaredSpectra','realizations', 'detunings','nAtoms','denominator','cloudX','cloudY','cloudZ','waist');


figure;
hold on;
plot(detunings,abs(resultSOmega/denominator+1).^2,'k--');
xlim([min(detunings) max(detunings)]);
ylabel('|s(\omega)|^2');
ylim([0 1]);