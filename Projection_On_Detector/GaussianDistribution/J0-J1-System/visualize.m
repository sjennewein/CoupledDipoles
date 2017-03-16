clear;
close all;

load('Spectra.mat');
k = 2*pi;
alpha = @(w) (6i.*pi/k^3)./(1-2i.*w);
S = @(w) 1+0.93*1i.*k.*alpha(w)./(pi.*waist.^2);
det = linspace(min(detunings),max(detunings),1000);
for iAtom = 1:numel(nAtoms)
    cAtom = nAtoms(iAtom);
    figure('name',[num2str(cAtom) ' J0J1Atoms']);
    hold on;
    plot(detunings,abs(SofOmegaSpectra{iAtom}/denominator+1).^2)
    plot(det,abs(S(det)).^2);
    ylim([0 1]);
    xlim([min(detunings) max(detunings)]);
    legend('|S(\omega)|^2 int. over lens', ...
           'Single atom lens theory', ...
           'Location','SE');
    xlabel('Detuning [\Delta/\Gamma]');
    ylabel('|s(\omega)|^2');
end