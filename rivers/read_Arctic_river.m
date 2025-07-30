function [discharge] = read_Arctic_river(fname)
% This function reads Arctic's rivers discharge from text file
opts = detectImportOptions(fname);
opts = setvartype(opts,'AverageOfDischarge','double');
T = readtable(fname,opts);
dum = T{:,2};
[m,n] = size(dum);
n=1;
for i =1:m-1
    if (i-1)/13-floor((i-1)/13) ~= 0
       j(n)=i;
       n=n+1;
    end
end
discharge = dum(j);
%
end
