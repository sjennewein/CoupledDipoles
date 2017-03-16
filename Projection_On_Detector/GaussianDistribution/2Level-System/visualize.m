clear;
close all;

data = load('2lvlResponse_sigma-sigma+.mat');

for i=1:numel(data.TransmissionSpectra)
    figure('name',num2str(data.nAtoms(i)));
    plot(data.detunings,data.TransmissionSpectra{i}/data.denominator);
    ylim([0 1]);
    xlim([-8 8]);
end