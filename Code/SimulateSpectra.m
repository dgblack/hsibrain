function [cSim, simSpecs] = SimulateSpectra(endmembers, nSpecs) 
% Parameters
nBasis = size(endmembers,2);
m = size(endmembers,1);

% Statistics from human data distributions
cMeans = [3.9100 0.6373 0.3418 0.0432 0.0788 0.0374 0.4721 0.0577 0.0729];
cStds = [7.9966 1.3128 0.4919 0.0939 0.2826 0.0723 0.5615 0.1643 0.1390];

% Generate random concentrations according to the distributions seen in
% human data
cSim = zeros(nBasis,nSpecs);
for i = 1:nBasis
    cSim(i,:) = normrnd(cMeans(i),cStds(i),1,size(cSim,2));

    % Remove negative concentrations and enforce sparsity
    cSim(i,cSim(i,:)<0.15) = 0;

    % Make sure there are no strong PpIX spectra with PpIX620 > PpIX634
    % This never happens in reality
    if i == 2
        while any(cSim(i,:) > cSim(1,:))
            idxs = cSim(i,:) > cSim(1,:);
            cSim(i,idxs) = cSim(i,idxs) / 2;
        end
    end
end

% Create spectra
simSpecs = endmembers * cSim;

% Add noise
mns = mean(simSpecs,2);
for i = 1:nSpecs
    noise = zeros(m,1);
    for j = 1:m
        mn = simSpecs(j,i);
        if mn <= 0.01
            mn = 0.001;
        end
        noise(j) = normrnd(0,sqrt(mn/300));
    end
    noise = abs(movmean(sgolayfilt(noise,15,31),16));
    simSpecs(:,i) = simSpecs(:,i) + noise;
end

end