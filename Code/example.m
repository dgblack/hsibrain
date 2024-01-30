%%% Example program showing the use of the unmixing functions

%% First load the endmembers
endmemberPath = 'EndmemberSpectra.csv';
[wavelengths, endmembers, endmemberNames] = loadEndmembers(endmemberPath);

% Plot the endmembers
figure
plot(wavelengths,endmembers,'LineWidth',1.25)
grid on
xlabel('Wavelength (nm)','FontSize',12)
ylabel('Intensity (a.u.)','FontSize',12)
set(gca,'FontSize',12)
xlim([wavelengths(1) wavelengths(end)])
legend(endmemberNames,'FontSize',12)
title('Endmember Spectra','FontSize',14)

%% Generate some spectra
nSpecs = 1000;
[cSim, simSpecs] = SimulateSpectra(endmembers, nSpecs);

% Plot the simulated spectra
figure
plot(wavelengths,simSpecs,'LineWidth',1.25)
grid on
xlabel('Wavelength (nm)','FontSize',12)
ylabel('Intensity (a.u.)','FontSize',12)
set(gca,'FontSize',12)
xlim([wavelengths(1) wavelengths(end)])
title('Simulated Spectra','FontSize',14)

%% Unmix with ISTA
lambdaIsta = 1.4; % Adjust this to obtain the desired performance
tic
cIsta = unmixISTA(simSpecs,endmembers,lambdaIsta);
dt = toc;

[reconErr, abundanceErr, L0, falsePos, SAM] = ComputePerformance(cIsta, cSim, simSpecs, endmembers);

disp(['ISTA runtime per spectrum: ' num2str(dt/nSpecs)]);
disp(['ISTA reconstruction MSE: ' num2str(reconErr(1)) ' +/- ' num2str(reconErr(2))]);
disp(['ISTA abundance MSE: ' num2str(abundanceErr(1)) ' +/- ' num2str(abundanceErr(2))]);
disp(['ISTA mean L0 norm: ' num2str(L0(1)) ' +/- ' num2str(L0(2))]);
disp(['ISTA false positive rate: ' num2str(falsePos(1))]);
disp(['ISTA spectral angle similarity: ' num2str(SAM(1)) ' +/- ' num2str(SAM(2))]);
disp('----')

%% Unmix with SNPR
lambdaSnpr = 0.35; % Adjust this to obtain the desired performance
tic
cSnpr = unmixSNPR(simSpecs,endmembers,lambdaSnpr);
dt = toc;

[reconErr, abundanceErr, L0, falsePos, SAM] = ComputePerformance(cSnpr, cSim, simSpecs, endmembers);

disp(['SNPR runtime per spectrum: ' num2str(dt/nSpecs)]);
disp(['SNPR reconstruction MSE: ' num2str(reconErr(1)) ' +/- ' num2str(reconErr(2))]);
disp(['SNPR abundance MSE: ' num2str(abundanceErr(1)) ' +/- ' num2str(abundanceErr(2))]);
disp(['SNPR mean L0 norm: ' num2str(L0(1)) ' +/- ' num2str(L0(2))]);
disp(['SNPR false positive rate: ' num2str(falsePos(1))]);
disp(['SNPR spectral angle similarity: ' num2str(SAM(1)) ' +/- ' num2str(SAM(2))]);


% Helper function to compute the performance metrics
function [reconErr, abundanceErr, L0, falsePos, SAM] = ComputePerformance(cCalc, cExp, mes, endmembers)
    badIdxs = any(isnan(cCalc));
    cCalc(:,badIdxs) = []; cExp(:,badIdxs) = []; mes(:,badIdxs) = [];

    recon = endmembers * cCalc;

    reconErr = vecnorm(recon-mes,2).^2;
    reconErr = [mean(reconErr) std(reconErr)];

    abundanceErr = mean((cCalc-cExp).^2,2);
    abundanceErr = [mean(abundanceErr) std(abundanceErr)];

    L0 = sum(cCalc>0);
    L0 = [mean(L0) std(L0)];

    falsePos = sum(sum(cCalc>0 & cExp == 0)) / (size(endmembers,2)*size(cCalc,2));

    SAM = dot(recon,mes) ./ (vecnorm(recon,2) .* vecnorm(mes,2));
    SAM = [mean(SAM) std(SAM)];
end