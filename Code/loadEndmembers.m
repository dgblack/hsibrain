function [wavelengths, endmembers, endmemberNames] = loadEndmembers(path)
    t = readtable(path);
    endmemberNames = t.Properties.VariableNames(2:end);
    wavelengths = table2array(t(:,1));
    endmembers = table2array(t(:,2:end));
end