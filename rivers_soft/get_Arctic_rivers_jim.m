function [Aa_full,fa_full2] = get_Arctic_rivers(fa_full);
%
% Load seven monthly Arctic river discharge time series from greatarcticrivers.org.  
% I converted the daily data to monthly and filled gaps in excel.  Units: m3/s.  This is from my excel files.
% To check file content: ncdisp('/aosc/horse/carton/rivers/dai-tren/coastal-stns-Vol-monthly.updated-May2019.nc');
%
[yenise] = textread('/aosc/horse/carton/rivers/arctic_rivers/yenise-1980-2019.txt','%f'); %#7
[lena] = textread('/aosc/horse/carton/rivers/arctic_rivers/lena-1980-2019.txt','%f');% #9
[ob] = load('/aosc/horse/carton/rivers/arctic_rivers/ob-1980-2019.txt','%f');% #13
[mackenzie] = load('/aosc/horse/carton/rivers/arctic_rivers/mackenzie-1980-2019.txt','%f'); % #19
[yukon] = load('/aosc/horse/carton/rivers/arctic_rivers/yukon-1980-2019.txt','%f'); % #24
[pechora] = load('/aosc/horse/carton/rivers/arctic_rivers/pechora-1980-2019.txt','%f'); % #31
[kolyma] = load('/aosc/horse/carton/rivers/arctic_rivers/kolyma-1980-2019.txt','%f'); % #35

Aa_full = zeros(480,7);
Aa_full(:,1) = yenise(:);
Aa_full(:,2) = lena(:);
Aa_full(:,3) = ob(:);
Aa_full(:,4) = mackenzie(:);
Aa_full(:,5) = yukon(:);
Aa_full(:,6) = pechora(:);
Aa_full(:,7) = kolyma(:);
%
% check Yenise: plot(fa_full(:,7));hold on;plot(Aa_full(:,1));hold off
% check Lena: plot(fa_full(:,9));hold on;plot(Aa_full(:,2));hold off
%  check total transport: mean(sum(Aa_full,2)+fa_full(:,31)+fa_full(:,33)+fa_full(:,37)+fa_full(:,64)+fa_full(:,97))*1.3 = 0.102Sv 
% The total should be ~0.10Sv according to Aagard and Carmack (1989).
% inflate the transports and substitute into fa_full.  We follow the guidence of AOMIP (https://www.whoi.edu/page.do?pid=30587) 
% to inflate the discharge of the 12 gauged rivers by 30% to account for the missing ungauged discharge.  
fa_full2 = fa_full;
fa_full2(:,7) = Aa_full(:,1)*1.3;% replace yenise
fa_full2(:,9) = Aa_full(:,2)*1.3;% replace lena
fa_full2(:,13)= Aa_full(:,3)*1.3;% replace ob
fa_full2(:,19)= Aa_full(:,4)*1.3;% replace mackenzie
fa_full2(:,24)= Aa_full(:,5)*1.3;% replace yukon
%fa_full(1:205,31)=fa_full(1:205,31)*1.3;
fa_full2(206:480,31)= Aa_full(206:480,6)*4.90/3.48; %replace part of Pechora (see below)
fa_full2(:,31)= fa_full(:,31)*1.3;% inflate Pechora 
fa_full2(:,35)= Aa_full(:,7)*1.3;% replace kolyma
fa_full2(:,33)=fa_full(:,33)*1.3;% inflate Khatanga
fa_full2(:,37)=fa_full(:,37)*1.3;% inflate Sev. Dvina
fa_full2(:,64)=fa_full(:,64)*1.3;% inflate Indigirka
fa_full2(:,97)=fa_full(:,97)*1.3;% inflate Olenek
fa_full2(:,104)=fa_full(:,104)*1.3;% inflate Yana
%
% Note: D-T bases their analysis on the gauge a t Oksino 141 km from Pechura delta.  I extend the time series using the record from Ust'-Tsil'ma (300+ km from delta)
%      since the discharge is lower at Ust'-Tsil'ma I multiply the U-T-based extension by the ratio of the mean transports at the two gauges, 4.90/3.48=1.41. 
%
%
end
