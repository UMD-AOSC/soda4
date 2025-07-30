function [discharge] = read_crdc_river(fname)
% This function reads CRDC's rivers discharge from text file
opts = detectImportOptions(fname);
T = readtable(fname,opts);
dum = T{:,4};
[m,n] = size(dum);
for i=1:m
    if dum(i) < 0.0
        dum(i) = NaN;
    end
end
discharge = dum;
%
end
