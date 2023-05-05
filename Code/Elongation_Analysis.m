% Elongation Analysis

clear;
% load('G:\Shared drives\Lim_Lab\Bomyi\elongation_manuscript\Code\initial_traces.mat')
load('/Users/bomyilim/Library/CloudStorage/GoogleDrive-bomyilim@seas.upenn.edu/Shared drives/Lim_Lab/Bomyi/elongation_manuscript/Code/initial_traces.mat')
%% Signal Extraction
thresh = 0;

hbmlp12 = signal_extraction(hbmlp12,4.5,2,'AP',thresh);
hbmlp13 = signal_extraction(hbmlp13,4.5,2,'AP',thresh);
hbmlp14 = signal_extraction(hbmlp14,4.5,2,'AP',thresh); 

spmlp12 = signal_extraction(spmlp12,4.5,2,'DV',thresh);
spmlp13 = signal_extraction(spmlp13,4.5,2,'DV',thresh);
spmlp14 = signal_extraction(spmlp14,4.5,2,'DV',thresh);

ssmlp12 = signal_extraction(ssmlp12,4.5,2,'DV',thresh);
ssmlp13 = signal_extraction(ssmlp13,4.5,2,'DV',thresh);
ssmlp14 = signal_extraction(ssmlp14,4.5,2,'DV',thresh);

sshsogmlp12 = signal_extraction(sshsogmlp12,4.5,2,'DV',thresh);
sshsogmlp13 = signal_extraction(sshsogmlp13,4.5,2,'DV',thresh);
sshsogmlp14 = signal_extraction(sshsogmlp14,4.5,2,'DV',thresh);

sshthsmlp12 = signal_extraction(sshthsmlp12,4.5,2,'DV',thresh);
sshthsmlp13 = signal_extraction(sshthsmlp13,4.5,2,'DV',thresh);
sshthsmlp14 = signal_extraction(sshthsmlp14,4.5,2,'DV',thresh);

sdsnamlp12 = signal_extraction(sdsnamlp12,4.5,2,'DV',thresh);
sdsnamlp13 = signal_extraction(sdsnamlp13,4.5,2,'DV',thresh);
sdsnamlp14 = signal_extraction(sdsnamlp14,4.5,2,'DV',thresh);

sdsnamglp12 = signal_extraction(sdsnamglp12,6,2,'DV',thresh);
sdsnamglp13 = signal_extraction(sdsnamglp13,6,2,'DV',thresh);
sdsnamglp14 = signal_extraction(sdsnamglp14,6,2,'DV',thresh);

sdsnamglgp12 = signal_extraction(sdsnamglgp12,8.5,2,'DV',thresh);
sdsnamglgp13 = signal_extraction(sdsnamglgp13,8.5,2,'DV',thresh);
sdsnamglgp14 = signal_extraction(sdsnamglgp14,8.5,2,'DV',thresh);

sisnamlp12 = signal_extraction(sisnamlp12,4.5,2,'DV',thresh);
sisnamlp13 = signal_extraction(sisnamlp13,4.5,2,'DV',thresh);
sisnamlp14 = signal_extraction(sisnamlp14,4.5,2,'DV',thresh);

sisnamglp12 = signal_extraction(sisnamglp12,6,2,'DV',thresh);
sisnamglp13 = signal_extraction(sisnamglp13,6,2,'DV',thresh);
sisnamglp14 = signal_extraction(sisnamglp14,6,2,'DV',thresh);

sdsogmlp12 = signal_extraction(sdsogmlp12,4.5,2,'DV',thresh);
sdsogmlp13 = signal_extraction(sdsogmlp13,4.5,2,'DV',thresh);
sdsogmlp14 = signal_extraction(sdsogmlp14,4.5,2,'DV',thresh);

sdsogmglp12 = signal_extraction(sdsogmglp12,6,2,'DV',thresh);
sdsogmglp13 = signal_extraction(sdsogmglp13,6,2,'DV',thresh);
sdsogmglp14 = signal_extraction(sdsogmglp14,6,2,'DV',thresh);

sdsogmglgp12 = signal_extraction(sdsogmglgp12,8.5,2,'DV',thresh);
sdsogmglgp13 = signal_extraction(sdsogmglgp13,8.5,2,'DV',thresh);
sdsogmglgp14 = signal_extraction(sdsogmglgp14,8.5,2,'DV',thresh);

sisogmlp12 = signal_extraction(sisogmlp12,4.5,2,'DV',thresh);
sisogmlp13 = signal_extraction(sisogmlp13,4.5,2,'DV',thresh);
sisogmlp14 = signal_extraction(sisogmlp14,4.5,2,'DV',thresh);

sisogmglp12 = signal_extraction(sisogmglp12,6,2,'DV',thresh);
sisogmglp13 = signal_extraction(sisogmglp13,6,2,'DV',thresh);
sisogmglp14 = signal_extraction(sisogmglp14,6,2,'DV',thresh);

% sdthsmlp12 = signal_extraction(sdthsmlp12,4.5,2,'DV',thresh);
% sdthsmlp13 = signal_extraction(sdthsmlp13,4.5,2,'DV',thresh);
sdthsmlp14 = signal_extraction(sdthsmlp14,4.5,2,'DV',thresh);

sdthsmglp12 = signal_extraction(sdthsmglp12,6,2,'DV',thresh);
sdthsmglp13 = signal_extraction(sdthsmglp13,6,2,'DV',thresh);
sdthsmglp14 = signal_extraction(sdthsmglp14,6,2,'DV',thresh);

sdthsmglgp12 = signal_extraction(sdthsmglgp12,8.5,2,'DV',thresh);
sdthsmglgp13 = signal_extraction(sdthsmglgp13,8.5,2,'DV',thresh);
sdthsmglgp14 = signal_extraction(sdthsmglgp14,8.5,2,'DV',thresh);

sithsmlp12 = signal_extraction(sithsmlp12,4.5,2,'DV',thresh);
sithsmlp13 = signal_extraction(sithsmlp13,4.5,2,'DV',thresh);
sithsmlp14 = signal_extraction(sithsmlp14,4.5,2,'DV',thresh);

sithsmglp12 = signal_extraction(sithsmglp12,6,2,'DV',thresh);
sithsmglp13 = signal_extraction(sithsmglp13,6,2,'DV',thresh);
sithsmglp14 = signal_extraction(sithsmglp14,6,2,'DV',thresh);

hbmyp12 = signal_extraction(hbmyp12,5.6,2,'AP',thresh);
hbmyp13 = signal_extraction(hbmyp13,5.6,2,'AP',thresh);
hbmyp14 = signal_extraction(hbmyp14,5.6,2,'AP',thresh); 

sdsogmyp12 = signal_extraction(sdsogmyp12,5.6,2,'DV',thresh);
sdsogmyp13 = signal_extraction(sdsogmyp13,5.6,2,'DV',thresh);
sdsogmyp14 = signal_extraction(sdsogmyp14,5.6,2,'DV',thresh);

sdsogmyep12 = signal_extraction(sdsogmyep12,3,2,'DV',thresh);
sdsogmyep13 = signal_extraction(sdsogmyep13,3,2,'DV',thresh);
sdsogmyep14 = signal_extraction(sdsogmyep14,3,2,'DV',thresh);

phbmlp = signal_extraction(phbmlp,4.5,2,'DV',thresh);
phbmyp = signal_extraction(phbmyp,5.6,2,'DV',thresh);

spsmlp13 = signal_extraction(spsmlp13,4.5,2,'DV',thresh);
spsmlp14 = signal_extraction(spsmlp14,4.5,2,'DV',thresh);

shsmlp14 = signal_extraction(shsmlp14,4.5,2,'DV',thresh);

shsmyep14 = signal_extraction(shsmyep14,3,2,'DV',thresh);

sdsogmmwp14 = signal_extraction(sdsogmmwp14,4.2,2,'DV',thresh);

sdsogmmwep14 = signal_extraction(sdsogmmwep14,3.4,2,'DV',thresh);

sdsogmkep14 = signal_extraction(sdsogmkep14,2.8,2,'DV',thresh);
%% Elongation Rate
hbrate13 = [];
hbrate14 = [];
sprate12 = [];
sprate13 = [];
sprate14 = [];
ssrate12 = [];
ssrate13 = [];
ssrate14 = [];
sshsograte12 = [];
sshsograte13 = [];
sshsograte14 = [];
sshthsrate12 = [];
sshthsrate13 = [];
sshthsrate14 = [];
sdsnamlprate12 = [];
sdsnamlprate13 = [];
sdsnamlprate14 = [];
sdsnamglprate12 = [];
sdsnamglprate13 = [];
sdsnamglprate14 = [];
sdsnamglgprate12 = [];
sdsnamglgprate13 = [];
sdsnamglgprate14 = [];
sisnamlprate12 = [];
sisnamlprate13 = [];
sisnamlprate14 = [];
sisnamglprate12 = [];
sisnamglprate13 = [];
sisnamglprate14 = [];
sdsogmlprate12 = [];
sdsogmlprate13 = [];
sdsogmlprate14 = [];
sdsogmglprate12 = [];
sdsogmglprate13 = [];
sdsogmglprate14 = [];
sdsogmglgprate12 = [];
sdsogmglgprate13 = [];
sdsogmglgprate14 = [];
sisogmlprate12 = [];
sisogmlprate13 = [];
sisogmlprate14 = [];
sisogmglprate12 = [];
sisogmglprate13 = [];
sisogmglprate14 = [];
sdthsmlprate12 = [];
sdthsmlprate13 = [];
sdthsmlprate14 = [];
sdthsmglprate12 = [];
sdthsmglprate13 = [];
sdthsmglprate14 = [];
sdthsmglgprate12 = [];
sdthsmglgprate13 = [];
sdthsmglgprate14 = [];
sithsmlprate12 = [];
sithsmlprate13 = [];
sithsmlprate14 = [];
sithsmglprate12 = [];
sithsmglprate13 = [];
sithsmglprate14 = [];

hbmyprate12 = [];
hbmyprate13 = [];
hbmyprate14 = [];
sdsogmyprate12 = [];
sdsogmyprate13 = [];
sdsogmyprate14 = [];
sdsogmyeprate12 = [];
sdsogmyeprate13 = [];
sdsogmyeprate14 = [];

phbmlprate = [];
phbmyprate = [];

spsrate13 = [];
spsrate14 = [];

shsrate14 = [];

shsmyeprate14 = [];

sdsogmmwprate14 = [];

sdsogmmweprate14 = [];

sdsogmkeprate14 = [];

for j = 1:2
    phbmyprate = [phbmyprate 5.6./[phbmyp(j).dT]];
end

for j = 1:length(sdthsmlp14)
    sdthsmlprate14 = [sdthsmlprate14 4.5./[sdthsmlp14(j).dT]];
end

for j = 1:length(sdsogmyp12)
%     hbrate12(j) = insert/nanmean(hbmlp12(j).dT);
    sdsogmyprate12 = [sdsogmyprate12 5.6./[sdsogmyp12(j).dT]];
    phbmlprate = [phbmlprate 4.5./[phbmlp(j).dT]];
end
for j = 1:length(spmlp14)
    sprate12 = [sprate12 4.5./[spmlp12(j).dT]];
    sprate13 = [sprate13 4.5./[spmlp13(j).dT]];
    sprate14 = [sprate14 4.5./[spmlp14(j).dT]];
    
    ssrate12 = [ssrate12 4.5./[ssmlp12(j).dT]];
    ssrate13 = [ssrate13 4.5./[ssmlp13(j).dT]];
    ssrate14 = [ssrate14 4.5./[ssmlp14(j).dT]];
    
    hbrate13 = [hbrate13 4.5./[hbmlp13(j).dT]];
    hbrate14 = [hbrate14 4.5./[hbmlp14(j).dT]];
    
%     sshsograte12 = [sshsograte12 4.5./[sshsogmlp12(j).dT]];
    sshsograte13 = [sshsograte13 4.5./[sshsogmlp13(j).dT]];
    sshsograte14 = [sshsograte14 4.5./[sshsogmlp14(j).dT]];
    
%     sshthsrate12 = [sshthsrate12 4.5./[sshthsrate12(j).dT]];
%     sshthsrate13 = [sshthsrate13 4.5./[sshthsrate13(j).dT]];
    sshthsrate14 = [sshthsrate14 4.5./[sshthsmlp14(j).dT]];

    sdsnamlprate12 = [sdsnamlprate12 4.5./[sdsnamlp12(j).dT]];
    sdsnamlprate13 = [sdsnamlprate13 4.5./[sdsnamlp13(j).dT]];
    sdsnamlprate14 = [sdsnamlprate14 4.5./[sdsnamlp14(j).dT]];
        
%     sdsnamlprate12 = [sdsnamlprate12 insert./[sdsnamlp12(j).dT]];
    sdsnamglprate13 = [sdsnamglprate13 6./[sdsnamglp13(j).dT]];
    sdsnamglprate14 = [sdsnamglprate14 6./[sdsnamglp14(j).dT]];
    
%     sdsnamglgprate12 = [sdsnamlprate12 insert./[sdsnamglgp12(j).dT]];
    sdsnamglgprate13 = [sdsnamglgprate13 8.5./[sdsnamglgp13(j).dT]];
    sdsnamglgprate14 = [sdsnamglgprate14 8.5./[sdsnamglgp14(j).dT]];
    
    sisnamlprate12 = [sisnamlprate12 4.5./[sisnamlp12(j).dT]];
    sisnamlprate13 = [sisnamlprate13 4.5./[sisnamlp13(j).dT]];
    sisnamlprate14 = [sisnamlprate14 4.5./[sisnamlp14(j).dT]];
    
%     sisnamglprate12 = [sisnamglprate12 6./[sisnamglp12(j).dT]];
    sisnamglprate13 = [sisnamglprate13 6./[sisnamglp13(j).dT]];
    sisnamglprate14 = [sisnamglprate14 6./[sisnamglp14(j).dT]];
    
%     sdsogmlprate12 = [sdsogmlprate12 4.5./[sdsogmlp12(j).dT]];
    sdsogmlprate13 = [sdsogmlprate13 4.5./[sdsogmlp13(j).dT]];
    sdsogmlprate14 = [sdsogmlprate14 4.5./[sdsogmlp14(j).dT]];
    
%     sdsogmglprate12 = [sdsogmglprate12 6./[sdsogmglp12(j).dT]];
    sdsogmglprate13 = [sdsogmglprate13 6./[sdsogmglp13(j).dT]];
    sdsogmglprate14 = [sdsogmglprate14 6./[sdsogmglp14(j).dT]];
    
%     sdsogmglgprate12 = [sdsogmglgprate12 8.5./[sdsogmglgp12(j).dT]];
    sdsogmglgprate13 = [sdsogmglgprate13 8.5./[sdsogmglgp13(j).dT]];
    sdsogmglgprate14 = [sdsogmglgprate14 8.5./[sdsogmglgp14(j).dT]];
    
%     sisogmlprate12 = [sisogmlprate12 4.5./[sisogmlp12(j).dT]];
    sisogmlprate13 = [sisogmlprate13 4.5./[sisogmlp13(j).dT]];
    sisogmlprate14 = [sisogmlprate14 4.5./[sisogmlp14(j).dT]];
    
    %     sisogmglprate12 = [sisogmglprate12 insert./[sisogmglp12(j).dT]];
    sisogmglprate13 = [sisogmglprate13 6./[sisogmglp13(j).dT]];
    sisogmglprate14 = [sisogmglprate14 6./[sisogmglp14(j).dT]];
    
%     sdthsmlprate12 = [sdthsmlprate12 insert./[sdthsmlp12(j).dT]];
%     sdthsmlprate13 = [sdthsmlprate13 4.5./[sdthsmlp13(j).dT]];
    
    
%     sdthsmglprate12 = [sdthsmglprate12 6./[sdthsmglp12(j).dT]];
    sdthsmglprate13 = [sdthsmglprate13 6./[sdthsmglp13(j).dT]];
    sdthsmglprate14 = [sdthsmglprate14 6./[sdthsmglp14(j).dT]];
    
%     sdthsmglgprate12 = [sdthsmglgprate12 8.5./[sdthsmglgp12(j).dT]];
    sdthsmglgprate13 = [sdthsmglgprate13 8.5./[sdthsmglgp13(j).dT]];
    sdthsmglgprate14 = [sdthsmglgprate14 8.5./[sdthsmglgp14(j).dT]];

%     sithsmlprate12 = [sithsmlprate12 insert./[sithsmlp12(j).dT]];
    sithsmlprate13 = [sithsmlprate13 4.5./[sithsmlp13(j).dT]];
    sithsmlprate14 = [sithsmlprate14 4.5./[sithsmlp14(j).dT]];
    
%     sithsmglprate12 = [sithsmglprate12 6./[sithsmglp12(j).dT]];
    sithsmglprate13 = [sithsmglprate13 6./[sithsmglp13(j).dT]];
    sithsmglprate14 = [sithsmglprate14 6./[sithsmglp14(j).dT]];
    
%     hbmyprate12 = [hbmyprate12 5.6./[hbmyp12(j).dT]];
    hbmyprate13 = [hbmyprate13 5.6./[hbmyp13(j).dT]];
    hbmyprate14 = [hbmyprate14 5.6./[hbmyp14(j).dT]];
    
    sdsogmyprate13 = [sdsogmyprate13 5.6./[sdsogmyp13(j).dT]];
    sdsogmyprate14 = [sdsogmyprate14 5.6./[sdsogmyp14(j).dT]];
    
    sdsogmyeprate12 = [sdsogmyeprate12 3./[sdsogmyep12(j).dT]];
    sdsogmyeprate13 = [sdsogmyeprate13 3./[sdsogmyep13(j).dT]];
    sdsogmyeprate14 = [sdsogmyeprate14 3./[sdsogmyep14(j).dT]];
    
    spsrate13 = [spsrate13 4.5./[spsmlp13(j).dT]];
    spsrate14 = [spsrate14 4.5./[spsmlp14(j).dT]];
    
    shsrate14 = [shsrate14 4.5./[shsmlp14(j).dT]];
    
    shsmyeprate14 = [shsmyeprate14 3./[shsmyep14(j).dT]];
    
    sdsogmmwprate14 = [sdsogmmwprate14 4.2./[sdsogmmwp14(j).dT]];
    
    sdsogmmweprate14 = [sdsogmmweprate14 3.4./[sdsogmmwep14(j).dT]];
    
    sdsogmkeprate14 = [sdsogmkeprate14 2.8./[sdsogmkep14(j).dT]];
end

%% discard rates above 10
sprate13(sprate13 > 10) = NaN;
sprate14(sprate14 > 10) = NaN;

ssrate13(ssrate13 > 10) = NaN;
ssrate14(ssrate14 > 10) = NaN;

hbrate13(hbrate13 > 10) = NaN;
hbrate14(hbrate14 > 10) = NaN;

sshsograte14(sshsograte14 > 10) = NaN;

sshthsrate14(sshthsrate14 > 10) = NaN;

sdsnamlprate13(sdsnamlprate13 > 10) = NaN;
sdsnamlprate14(sdsnamlprate14 > 10) = NaN;

sdsnamglprate13(sdsnamglprate13 > 10) = NaN;
sdsnamglprate14(sdsnamglprate14 > 10) = NaN;

sdsnamglgprate13(sdsnamglgprate13 > 10) = NaN;
sdsnamglgprate14(sdsnamglgprate14 > 10) = NaN;

sisnamlprate13(sisnamlprate13 > 10) = NaN;
sisnamlprate14(sisnamlprate14 > 10) = NaN;

sisnamglprate13(sisnamglprate13 > 10) = NaN;
sisnamglprate14(sisnamglprate14 > 10) = NaN;

sdsogmlprate13(sdsogmlprate13 > 10) = NaN;
sdsogmlprate14(sdsogmlprate14 > 10) = NaN;

sdsogmglprate13(sdsogmglprate13 > 10) = NaN;
sdsogmglprate14(sdsogmglprate14 > 10) = NaN;

sdsogmglgprate13(sdsogmglgprate13 > 10) = NaN;
sdsogmglgprate14(sdsogmglgprate14 > 10) = NaN;

sisogmlprate13(sisogmlprate13 > 10) = NaN;
sisogmlprate14(sisogmlprate14 > 10) = NaN;

sisogmglprate13(sisogmglprate13 > 10) = NaN;
sisogmglprate14(sisogmglprate14 > 10) = NaN;

% sdthsmlprate13(sdthsmlprate13 > 10) = NaN;
sdthsmlprate14(sdthsmlprate14 > 10) = NaN;

sdthsmglprate13(sdthsmglprate13 > 10) = NaN;
sdthsmglprate14(sdthsmglprate14 > 10) = NaN;

sdthsmglgprate13(sdthsmglgprate13 > 10) = NaN;
sdthsmglgprate14(sdthsmglgprate14 > 10) = NaN;

sithsmlprate13(sithsmlprate13 > 10) = NaN;
sithsmlprate14(sithsmlprate14 > 10) = NaN;

sithsmglprate13(sithsmglprate13 > 10) = NaN;
sithsmglprate14(sithsmglprate14 > 10) = NaN;

hbmyprate12(hbmyprate12 > 10) = NaN;
hbmyprate13(hbmyprate13 > 10) = NaN;
hbmyprate14(hbmyprate14 > 10) = NaN;

sdsogmyprate12(sdsogmyprate12 > 10) = NaN;
sdsogmyprate13(sdsogmyprate13 > 10) = NaN;
sdsogmyprate14(sdsogmyprate14 > 10) = NaN;

sdsogmyeprate12(sdsogmyeprate12 > 10) = NaN;
sdsogmyeprate13(sdsogmyeprate13 > 10) = NaN;
sdsogmyeprate14(sdsogmyeprate14 > 10) = NaN;

phbmlprate(phbmlprate > 10) = NaN;
phbmyprate(phbmyprate > 10) = NaN;

spsrate13(spsrate13 > 10) = NaN;
spsrate14(spsrate14 > 10) = NaN;

shsrate14(shsrate14 > 10) = NaN;

shsmyeprate14(shsmyeprate14 > 10) = NaN;

sdsogmmwprate14(sdsogmmwprate14 > 10) = NaN;

sdsogmmweprate14(sdsogmmweprate14 > 10) = NaN;

sdsogmkeprate14(sdsogmkeprate14 > 10) = NaN;
%%
spErate13 = 4.5/nanmean([spmlp13.dT]);
spErate14 = 4.5/nanmean([spmlp14.dT]);
ssErate13 = 4.5/nanmean([ssmlp13.dT]);
ssErate14 = 4.5/nanmean([ssmlp14.dT]);
hbErate13 = 4.5/nanmean([hbmlp13.dT]);
hbErate14 = 4.5/nanmean([hbmlp14.dT]);
sshsogErate13 = 4.5/nanmean([sshsogmlp13.dT]);
sshsogErate14 = 4.5/nanmean([sshsogmlp14.dT]);
% sshthsErate13 = 4.5/nanmean([sshthsmlp13.dT]);
% sshthsErate14 = 4.5/nanmean([sshthsmlp14.dT]);
sshsogErate13 = 4.5/nanmean([sshsogmlp13.dT]);
sshsogErate14 = 4.5/nanmean([sshsogmlp14.dT]);

spmad13 = mad(sprate13,1);
spmad14 = mad(sprate14,1);
ssmad13 = mad(ssrate13,1);
ssmad14 = mad(ssrate14,1);
hbmad13 = mad(hbrate13,1);
hbmad14 = mad(hbrate14,1);
sshsogmad13 = mad(sshsograte13,1);
sshsogmad14 = mad(sshsograte14,1);
sdsnamlpmad13 = mad(sdsnamlprate13,1);
sdsnamlpmad14 = mad(sdsnamlprate14,1);
sdsnamglpmad13 = mad(sdsnamglprate13,1);
sdsnamglpmad14 = mad(sdsnamglprate14,1);

figure;
errorbar(1,nanmedian(hbrate14),hbmad14,'ko'); hold on
errorbar(2,nanmedian(sprate14),spmad14,'bo');
errorbar(3,nanmedian(ssrate14),ssmad14,'ro');
errorbar(4,nanmedian(sshsograte14),sshsogmad14,'go');
errorbar(5,nanmedian(sdsnamlprate14),sdsnamlpmad14,'co');

%% Elongation Boxplot 14

data = [hbrate14 sprate14 ssrate14 sshsograte14 sshthsrate14 sdsnamlprate14 sdsnamglprate14 sdsnamglgprate14 sisnamlprate14 sisnamglprate14 sdsogmlprate14 sdsogmglprate14 sdsogmglgprate14 sisogmlprate14 ...
    sisogmglprate14 sdthsmlprate14 sdthsmglprate14 sdthsmglgprate14 sithsmlprate14 sithsmglprate14];
grp = [ones(1, length(hbrate14)) 2*ones(1,length(sprate14)) 3*ones(1,length(ssrate14)) 4*ones(1,length(sshsograte14)) 5* ones(1,length(sshthsrate14)) 6*ones(1,length(sdsnamlprate14)) 7*ones(1,length(sdsnamglprate14)) ...
    8*ones(1,length(sdsnamglgprate14)) 9*ones(1,length(sisnamlprate14)) 10*ones(1,length(sisnamglprate14)) 11*ones(1,length(sdsogmlprate14)) 12*ones(1,length(sdsogmglprate14)) 13*ones(1,length(sdsogmglgprate14)) ...
    14*ones(1,length(sisogmlprate14)) 15*ones(1,length(sisogmglprate14)) 16*ones(1,length(sdthsmlprate14)) 17*ones(1,length(sdthsmglprate14)) 18*ones(1,length(sdthsmglgprate14)) 19*ones(1,length(sithsmlprate14)) 20*ones(1,length(sithsmglprate14))];
figure;
boxplot(data,grp,'Labels',{'Hb','Sna prox','Sna shadow','Sna shadow sog','Sna shadow ths','SD sna mlp', 'SD sna mglp', 'SD sna mglgp', 'SI sna mlp','SI sna mglp','SD sog mlp','SD sog mglp', 'SD sog mglgp', 'SI sog mlp' ...
    'SI sog mglp','SD ths mlp','SD ths mglp', 'SD ths mglgp', 'SI ths mlp','SI ths mglp'}); hold on
ylim([0 7])
set(gca,'XTickLabelRotation',45)

% xCenter = 1; 
% spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
% 
% for i = 1:length(hbrate14)
%     plot(rand(length(hbrate14(i)))*spread -(spread/2) + xCenter, hbrate14(i), 'b.','linewidth', 1)
% end

%% Elongation Boxplot 13

data = [hbrate13 hbrate14 sprate13 sprate14 ssrate13 ssrate14 sdsnamlprate13 sdsnamlprate14 sdsnamglprate13 sdsnamglprate14 sdsnamglgprate13 sdsnamglgprate14 sdsogmlprate13 sdsogmlprate14 sdsogmglprate13 sdsogmglprate14 ...
    sdsogmglgprate13 sdsogmglgprate14 sdthsmlprate13 sdthsmlprate14 sdthsmglprate13 sdthsmglprate14 sdthsmglgprate13 sdthsmglgprate14 sisnamlprate13 sisnamlprate14 sisnamglprate13 sisnamglprate14 sisogmlprate13 sisogmlprate14 ...
    sisogmglprate13 sisogmglprate14 sithsmlprate13 sithsmlprate14 sithsmglprate13 sithsmglprate14];
grp = [ones(1, length(hbrate13)) 2*ones(1,length(hbrate14)) 3*ones(1,length(sprate13)) 4*ones(1,length(sprate14)) 5*ones(1,length(ssrate13)) ...
    6*ones(1,length(ssrate14)) 7*ones(1,length(sdsnamlprate13)) 8*ones(1,length(sdsnamlprate14)) 9*ones(1,length(sdsnamglprate13)) 10*ones(1,length(sdsnamglprate14)) 11*ones(1,length(sdsnamglgprate13)) 12*ones(1,length(sdsnamglgprate14))...
    13*ones(1,length(sdsogmlprate13)) 14*ones(1,length(sdsogmlprate14)) 15*ones(1,length(sdsogmglprate13)) 16*ones(1,length(sdsogmglprate14)) 17*ones(1,length(sdsogmglgprate13)) 18*ones(1,length(sdsogmglgprate14)) ...
    19*ones(1, length(sdthsmlprate13)) 20*ones(1,length(sdthsmlprate14)) 21*ones(1,length(sdthsmglprate13)) 22*ones(1,length(sdthsmglprate14)) 23*ones(1,length(sdthsmglgprate13)) 24*ones(1,length(sdthsmglgprate14)) 25*ones(1,length(sisnamlprate13)) ...
    26*ones(1,length(sisnamlprate14)) 27*ones(1,length(sisnamglprate13)) 28*ones(1,length(sisnamglprate14)) 29*ones(1,length(sisogmlprate13)) 30*ones(1,length(sisogmlprate14)) ...
    31*ones(1,length(sisogmglprate13)) 32*ones(1,length(sisogmglprate14)) 33*ones(1,length(sithsmlprate13)) 34*ones(1,length(sithsmlprate14)) 35*ones(1,length(sithsmglprate13)) 36*ones(1,length(sithsmglprate14))];
% figure;
% boxplot(data,grp,'Labels',{'Hb','14','Sna prox','14','Sna shadow','14','SD sna mlp', '14','SD sna mglp','14', 'SD sna mglgp', '14', 'SD sog mlp','14','SD sog mglp','14', ...
%     'SD sog mglgp', '14', 'SD ths mlp','14','SD ths mglp','14','SD ths mglgp', '14','SI sna mlp','14', ...
%     'SI sna mglp','14','SI sog mlp','14', 'SI sog mlp', '14', 'SI sog mglp','14','SI ths mglp','14'}); hold on
% ylim([0 7])
% set(gca,'XTickLabelRotation',45)

%% Histogram
figure; histogram(log(hbrate14),100); hold on
histogram(log(sprate14),100); hold on
histogram(log(ssrate14),100);
histogram(log(sshsograte14),100);
histogram(log(sdsnamlprate14),100);
histogram(log(sdsnamglprate14),100);
histogram(log(sisnamlprate14),100);
histogram(log(sisnamglprate14),100);
histogram(log(sdsogmlprate14),100);
histogram(log(sdsogmglprate14),100);
histogram(log(sisogmlprate14),100);

%% See Traces of High rates
%{
placeholder = sisnamlp14;
for j = 1:length(placeholder)
    indx{j} = find(placeholder(j).rate > 10);
end
q = 0;
figure;
for j = 1:length(indx)
    for i = 1:length(indx{j})
        q = q +1;
        subplot(5,10,q)
        yyaxis left
        plot(placeholder(j).intT,placeholder(j).intM(:,placeholder(j).active(indx{j}(i))),'g'); hold on
        hline(.3*max(placeholder(j).M(:,placeholder(j).active(indx{j}(i)))),{'color',rgb('green')});
        plot(placeholder(j).actM(indx{j}(i)),.3*max(placeholder(j).M(:,placeholder(j).active(indx{j}(i)))),'go');
        yyaxis right
        plot(placeholder(j).intT,placeholder(j).intP(:,placeholder(j).active(indx{j}(i))),'r'); hold on
        hline(.3*max(placeholder(j).P(:,placeholder(j).active(indx{j}(i)))),{'color',rgb('red')});
        xlim([0 30])
        plot(placeholder(j).actP(indx{j}(i)),.3*max(placeholder(j).P(:,placeholder(j).active(indx{j}(i)))),'ro');
        text(1,1,num2str(placeholder(j).rate(indx{j}(i))))
    end
end
%}

%% Plot traces

marker = sdsnamlp13;
figure;
A = 0.3;
j = 1;
q = 0;
for i = 1:1:200
    q = q +1;
    subplot(5,8,q)
    yyaxis left
    plot(marker(j).intT,marker(j).intM(:,marker(j).active(i)),'g'); hold on
    hline(marker(j).threshM,'g');

%     vline(5);
    plot(marker(j).actM(i),A*max(marker(j).M(:,marker(j).active(i))),'go');
    set(gca,'ycolor',['g'])
%     ylim([0 1])
    yyaxis right
    plot(marker(j).intT,marker(j).intP(:,marker(j).active(i)),'r'); hold on
    hline(marker(j).threshP,'r');
%     yline(A*max(spsmlp14(j).P(:,spsmlp14(j).active(i))),'color',rgb('red'));
    xlim([0 35])
    plot(marker(j).actP(i),A*max(marker(j).P(:,marker(j).active(i))),'ro');
    set(gca,'ycolor',['r'])

end

%% NC
data = [ sdsnamlprate14 sdsogmlprate14 sdthsmlprate14];
grp = [ones(1,length(sdsnamlprate14)) 2*ones(1,length(sdsogmlprate14)) 3*ones(1,length(sdthsmlprate14))];
figure;
boxplot(data,grp,'Labels',{'SD sna mlp','SD sog mlp', 'SD ths mlp'}); hold on
set(gca,'XTickLabelRotation',45)
ylim([0 7])


%% Cross Correlation
data = [[hbmlp14.xcorr] [spmlp14.xcorr] [ssmlp14.xcorr] [sshsogmlp14.xcorr] [sdsnamlp14.xcorr] [sdsnamglp14.xcorr] [sdsnamglgp14.xcorr]  [sdsogmlp14.xcorr] [sdsogmglp14.xcorr] [sdsogmglgp14.xcorr] ...
     [sdthsmlp14.xcorr] [sdthsmglp14.xcorr] [sdthsmglgp14.xcorr] [sisnamlp14.xcorr] [sisnamglp14.xcorr] [sisogmlp14.xcorr] [sisogmglp14.xcorr] [sithsmlp14.xcorr] [sithsmglp14.xcorr]];
grp = [ones(1,length([hbmlp14.xcorr])) 2*ones(1,length([spmlp14.xcorr])) 3*ones(1,length([ssmlp14.xcorr])) 4*ones(1,length([sshsogmlp14.xcorr])) 5*ones(1,length([sdsnamlp14.xcorr])) 6*ones(1,length([sdsnamglp14.xcorr])) ...
    7*ones(1,length([sdsnamglgp14.xcorr])) 8*ones(1,length([sdsogmlp14.xcorr])) 9*ones(1,length([sdsogmglp14.xcorr])) 10*ones(1,length([sdsogmglgp14.xcorr])) 11*ones(1,length([sdthsmlp14.xcorr])) 12*ones(1,length([sdthsmglp14.xcorr])) ...
    13*ones(1,length([sdthsmglgp14.xcorr])) 14*ones(1,length([sisnamlp14.xcorr])) 15*ones(1,length([sisnamglp14.xcorr]))  16*ones(1,length([sisogmlp14.xcorr])) 17*ones(1,length([sisogmglp14.xcorr])) ...
    18*ones(1,length([sithsmlp14.xcorr])) 19*ones(1,length([sithsmglp14.xcorr]))];

figure;
boxplot(data,grp,'Labels',{'Hb','Sna prox','Sna shadow','Sna shadow sog','SD sna mlp', 'SD sna mglp', 'SD sna mglgp', 'SD sog mlp','SD sog mglp','SD sog mglgp', 'SD ths mlp','SD ths mglp',...
    'SD ths mglgp','SI sna mlp','SI sna mglp','SI sog mlp','SI sog mglp','SI ths mlp','SI ths mglp'}); hold on
set(gca,'XTickLabelRotation',45)
ylim([0 7])

figure; plot(hbrate14,[hbmlp14.xcorr],'ko'); xlim([0 10]); ylim([0 10]); xlabel('leading edge'); ylabel('xcorr'); title('hb')
figure; plot(sprate14,[spmlp14.xcorr],'ko'); xlim([0 10]);ylim([0 10]);xlabel('leading edge'); ylabel('xcorr'); title('sna prox')
figure; plot(ssrate14,[ssmlp14.xcorr],'ko'); xlim([0 10]);ylim([0 10]);xlabel('leading edge'); ylabel('xcorr'); title('sna sh')
figure; plot(sshsograte14,[sshsogmlp14.xcorr],'ko'); xlim([0 10]);ylim([0 10]);xlabel('leading edge'); ylabel('xcorr'); title('sna sh sog')

figure; plot(sdsnamlprate14,[sdsnamlp14.xcorr],'ko'); xlim([0 10]);ylim([0 10]);xlabel('leading edge'); ylabel('xcorr'); title('SD sna MLP')
figure; plot(sdsnamglprate14,[sdsnamglp14.xcorr],'ko'); xlim([0 10]);ylim([0 10]);xlabel('leading edge'); ylabel('xcorr'); title('SD sna MgLP')
figure; plot(sdsnamglgprate14,[sdsnamglgp14.xcorr],'ko'); xlim([0 10]);ylim([0 10]);xlabel('leading edge'); ylabel('xcorr'); title('SD sna MgLgP')

figure; plot(sisnamlprate14,[sisnamlp14.xcorr],'ko'); xlim([0 10]);ylim([0 10]);xlabel('leading edge'); ylabel('xcorr'); title('SI sna MLP')
figure; plot(sisnamglprate14,[sisnamglp14.xcorr],'ko'); xlim([0 10]);ylim([0 10]);xlabel('leading edge'); ylabel('xcorr'); title('SI sna MgLP')

figure; plot(sdsogmlprate14,[sdsogmlp14.xcorr],'ko'); xlim([0 10]);ylim([0 10]);xlabel('leading edge'); ylabel('xcorr'); title('SD sog MLP')
figure; plot(sdsogmglprate14,[sdsogmglp14.xcorr],'ko'); xlim([0 10]);ylim([0 10]);xlabel('leading edge'); ylabel('xcorr'); title('SD sog MgLP')
figure; plot(sdsogmglgprate14,[sdsogmglgp14.xcorr],'ko'); xlim([0 10]);ylim([0 10]);xlabel('leading edge'); ylabel('xcorr'); title('SD sog MgLgP')

figure; plot(sisogmlprate14,[sisogmlp14.xcorr],'ko'); xlim([0 10]);ylim([0 10]);xlabel('leading edge'); ylabel('xcorr'); title('SI sog MLP')
figure; plot(sisogmglprate14,[sisogmglp14.xcorr],'ko'); xlim([0 10]);ylim([0 10]);xlabel('leading edge'); ylabel('xcorr'); title('SI sog MgLP')

figure; plot(sdthsmlprate14,[sdthsmlp14.xcorr],'ko'); xlim([0 10]);ylim([0 10]);xlabel('leading edge'); ylabel('xcorr'); title('SD ths MLP')
figure; plot(sdthsmglprate14,[sdthsmglp14.xcorr],'ko'); xlim([0 10]);ylim([0 10]);xlabel('leading edge'); ylabel('xcorr'); title('SD ths MgLP')
figure; plot(sdthsmglgprate14,[sdthsmglgp14.xcorr],'ko'); xlim([0 10]);ylim([0 10]);xlabel('leading edge'); ylabel('xcorr'); title('SD ths MgLgP')

figure; plot(sithsmlprate14,[sithsmlp14.xcorr],'ko'); xlim([0 10]);ylim([0 10]);xlabel('leading edge'); ylabel('xcorr'); title('SI ths MLP')
figure; plot(sithsmglprate14,[sithsmglp14.xcorr],'ko'); xlim([0 10]);ylim([0 10]);xlabel('leading edge'); ylabel('xcorr'); title('SI ths MgLP')

%%
data = [hbrate14 sprate14 ssrate14 sshsograte14 sdsnamlprate14 sdsnamglprate14 sdsnamglgprate14 sdsogmlprate14 sdsogmglprate14 sdsogmglgprate14  ...
     sdthsmlprate14 sdthsmglprate14 sdthsmglgprate14 sisnamlprate14 sisnamglprate14 sisogmlprate14 sisogmglprate14 sithsmlprate14 sithsmglprate14];
grp = [ones(1, length(hbrate14)) 2*ones(1,length(sprate14)) 3*ones(1,length(ssrate14)) 4*ones(1,length(sshsograte14)) 5*ones(1,length(sdsnamlprate14)) 6*ones(1,length(sdsnamglprate14)) ...
    7*ones(1,length(sdsnamglgprate14)) 8*ones(1,length(sdsogmlprate14)) 9*ones(1,length(sdsogmglprate14)) 10*ones(1,length(sdsogmglgprate14)) 11*ones(1,length(sdthsmlprate14)) 12*ones(1,length(sdthsmglprate14)) ...
    13*ones(1,length(sdthsmglgprate14)) 14*ones(1,length(sisnamlprate14)) 15*ones(1,length(sisnamglprate14)) 16*ones(1,length(sisogmlprate14)) 17*ones(1,length(sisogmglprate14)) 18*ones(1,length(sithsmlprate14)) ...
    19*ones(1,length(sithsmglprate14))];
figure;
boxplot(data,grp,'Labels',{'Hb','Sna prox','Sna shadow','Sna shadow sog','SD sna mlp', 'SD sna mglp', 'SD sna mglgp', 'SD sog mlp','SD sog mglp','SD sog mglgp','SD ths mlp','SD ths mglp','SD ths mglgp','SI sna mlp' ...
    'SI sna mglp','SI sog mlp','SI sog mglp','SI ths mlp','SI ths mglp'}); hold on
ylim([0 7])
set(gca,'XTickLabelRotation',45)

%% Paper boxplot
% Cell = {sdsnamlprate14;sisnamlprate14;sdsnamglprate14;sisnamglprate14; sdsogmlprate14; sisogmlprate14; sdsogmglprate14;sisogmglprate14; sdthsmlprate14; sithsmlprate14; ...
%     sdthsmglprate14;sithsmglprate14};

Cell = {[sdsogmlp14.OutputP];[sdsogmglp14.OutputP];[sdsogmyep14.OutputP];[sdsogmyp14.OutputP]};
% Cell = {sdsogmglprate14; sdsogmyeprate14; sdsogmmweprate14; sdsogmyprate14; sdsogmmwprate14};
% Cell = {sdsnamlprate14; sdsogmlprate14; sdthsmlprate14};
% Cell = {hbrate14,sprate14,ssrate14,sdsogmlprate14};
% Cell = {sdsogmlprate14;sdsogmglprate14;sdsogmglgprate14};
dataCell = [];
grp = [];


%get 10 - 90 percentile
for i = 1:length(Cell)
    tmp = prctile([Cell{i}],[10,90]);
    Cell{i}(Cell{i} < tmp(1)) = [];
    Cell{i}(Cell{i} > tmp(2)) = [];
end

for i = 1:length(Cell)
    dataCell=[dataCell, Cell{i}];
    grp=[grp, ones(1, length(Cell{i}))*i];
end

figure,boxplot(dataCell,grp,'Labels',{'sd sog mlp','sd sog mglp','sd sog my(e)p','sd sog myp'},'notch','on')
% figure,boxplot(dataCell,grp,'Labels',{'sna prox?','sna shadow?','sna prox again'},'notch','on')%,'Labels',string)%,'LabelOrientation','inline')
% set(gca,'XTickLabelRotation',45)
% ylim([0 10000])

hold on
% 
% plotSpread(Cell)
% %
xCenter = 1:numel(Cell); 
spread = 0.3; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(Cell)
    plot(rand(size(Cell{(i)}))*spread -(spread/2) + xCenter(i), Cell{(i)}, 'b.','linewidth', 1)
end


%% Quick boxplot

Cell = {sdsnamlprate14; sdsogmlprate14; sdthsmlprate14};

dataCell = [];
grp = [];


%get 10 - 90 percentile
for i = 1:length(Cell)
    tmp = prctile([Cell{i}],[10,90]);
    Cell{i}(Cell{i} < tmp(1)) = [];
    Cell{i}(Cell{i} > tmp(2)) = [];
end

for i = 1:length(Cell)
    dataCell=[dataCell, Cell{i}];
    grp=[grp, ones(1, length(Cell{i}))*i];
end

figure,boxplot(dataCell,grp,'notch','on')
% figure,boxplot(dataCell,grp,'Labels',{'sna prox?','sna shadow?','sna prox again'},'notch','on')%,'Labels',string)%,'LabelOrientation','inline')
% set(gca,'XTickLabelRotation',45)
% ylim([0 10])
hold on
% 
% plotSpread(Cell)
% %
xCenter = 1:numel(Cell); 
spread = 0.3; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(Cell)
    plot(rand(size(Cell{(i)}))*spread -(spread/2) + xCenter(i), Cell{(i)}, 'b.','linewidth', 1)
end

%% Mean and std dev
% celldata = {hbrate14;sprate14;ssrate14; sdsnamlprate14; sisnamlprate14; sdsnamglprate14; sisnamglprate14; sdsnamglgprate14; sdsogmlprate14; sisogmlprate14; sdsogmglprate14; sisogmglprate14; ...
%     sdsogmglgprate14; sdthsmlprate14; sithsmlprate14; sdthsmglprate14; sithsmglprate14; sdthsmglgprate14};
celldata = {sdsnamglgprate14;sdsogmglgprate14;sdthsmglgprate14};
for i = 1:length(celldata)
    meandata(i) = nanmean(celldata{i});
    stddata(i) = nanstd(celldata{i});
    vardata(i) = nanvar(celldata{i});
end
figure;
errorbar([1:3],meandata(1:3),stddata(1:3),'ko')

figure;
plot(stddata,'ko')
xlim([.5 3.5])

%% Output vs Elongation rate
clear output outputM
marker = sdsnamlp14;
markerrate = sdsnamlprate14;

for j = 1:length(marker)
    for i = 1:length(marker(j).active)
        output{j}(i) = trapz(marker(j).intT,marker(j).intM(:,marker(j).active(i)));
    end
end
outputM = [output{:}];
figure;
plot(outputM,markerrate,'ko')
ylim([0 10])
ylabel('Rate (kb/min)');
xlabel('MS2 Output (A.U.)');
title('Tmp')

%%
test = corrcoef(sdsnamlprate14,[sdsnamlp14.OutputM],'Rows','complete');

%% Enhancer Strenght Comparrison (must use same power levels)
%Sog MLP SD and SI
%Ths MgLP SD and SI

%All Sog can be compared and all Sna except Sna MLP

%Ths can be compared within ths except ths MLP

sisnamlp14outputM = nanmean([sisnamlp14.OutputM]);
stdsisnamlp14outputM = nanstd([sisnamlp14.OutputM])/sqrt(length([sisnamlp14.OutputM]));

sdsnamglp14outputM = nanmean([sdsnamglp14.OutputM]);
stdsdsnamglp14outputM = nanstd([sdsnamglp14.OutputM])/sqrt(length([sdsnamglp14.OutputM]));
sisnamglp14outputM = nanmean([sisnamglp14.OutputM]);
stdsisnamglp14outputM = nanstd([sisnamglp14.OutputM])/sqrt(length([sisnamglp14.OutputM]));

sdsnamglgp14outputM = nanmean([sdsnamglgp14.OutputM]);
stdsdsnamglgp14outputM = nanstd([sdsnamglgp14.OutputM])/sqrt(length([sdsnamglgp14.OutputM]));

sdsogmlp14outputM = nanmean([sdsogmlp14.OutputM]);
stdsdsogmlp14outputM = nanstd([sdsogmlp14.OutputM])/sqrt(length([sdsogmlp14.OutputM]));
sisogmlp14outputM = nanmean([sisogmlp14.OutputM]);
stdsisogmlp14outputM = nanstd([sisogmlp14.OutputM])/sqrt(length([sisogmlp14.OutputM]));

sdsogmglp14outputM = nanmean([sdsogmglp14.OutputM]);
stdsdsogmglp14outputM = nanstd([sdsogmglp14.OutputM])/sqrt(length([sdsogmglp14.OutputM]));
sisogmglp14outputM = nanmean([sisogmglp14.OutputM]);
stdsisogmglp14outputM = nanstd([sisogmglp14.OutputM])/sqrt(length([sisogmglp14.OutputM]));

sdthsmlp14outputM = nanmean([sdthsmlp14.OutputM]);
stdsdthsmlp14outputM = nanstd([sdthsmlp14.OutputM])/sqrt(length([sdthsmlp14.OutputM]));
sithsmlp14outputM = nanmean([sithsmlp14.OutputM]);
stdsithsmlp14outputM = nanstd([sithsmlp14.OutputM])/sqrt(length([sithsmlp14.OutputM]));

sdthsmglp14outputM = nanmean([sdthsmglp14.OutputM]);
stdsdthsmglp14outputM = nanstd([sdthsmglp14.OutputM])/sqrt(length([sdthsmglp14.OutputM]));
sithsmglp14outputM = nanmean([sithsmglp14.OutputM]);
stdsithsmglp14outputM = nanstd([sithsmglp14.OutputM])/sqrt(length([sithsmglp14.OutputM]));

ssmlp14outputM = nanmean([ssmlp14.OutputM]);
stdssmlp14outputM = nanstd([ssmlp14.OutputM])/sqrt(length([ssmlp14.OutputM]));
sshsogmlp14outputM = nanmean([sshsogmlp14.OutputM]);
stdsshsogmlp14outputM = nanstd([sshsogmlp14.OutputM])/sqrt(length([sshsogmlp14.OutputM]));
sdsnamlp14outputM = nanmean([sdsnamlp14.OutputM]);
stdsdsnamlp14outputM = nanstd([sdsnamlp14.OutputM])/sqrt(length([sdsnamlp14.OutputM]));
hbmlp14outputM = nanmean([hbmlp14.OutputM]);
stdhbmlp14outputM = nanstd([hbmlp14.OutputM])/sqrt(length([hbmlp14.OutputM]));
spmlp14outputM = nanmean([spmlp14.OutputM]);
stdspmlp14outputM = nanstd([spmlp14.OutputM])/sqrt(length([spmlp14.OutputM]));
sshthsmlp14outputM = nanmean([sshthsmlp14.OutputM]);
stdsshthsmlp14outputM = nanstd([sshthsmlp14.OutputM])/sqrt(length([sshthsmlp14.OutputM]));



figure;
errorbar(1,sdsnamglp14outputM,stdsdsnamglp14outputM,'ko'); hold on %recheck laser power 
errorbar(2,sdsnamglgp14outputM,stdsdsnamglgp14outputM,'ko'); %recheck laser power
errorbar(3,sdsogmlp14outputM,stdsdsogmlp14outputM,'ko'); 
errorbar(4,sdsogmglp14outputM,stdsdsogmglp14outputM,'ko');

errorbar(5,sisnamlp14outputM,stdsisnamlp14outputM,'ko');
errorbar(6,sisnamglp14outputM,stdsisnamglp14outputM,'ko');
errorbar(7,sisogmlp14outputM,stdsisogmlp14outputM,'ko');
errorbar(8,sisogmglp14outputM,stdsisogmglp14outputM,'ko');
xlim([.5 8.5])

figure;
errorbar(1,sdthsmlp14outputM,stdsdthsmlp14outputM,'ko');hold on
errorbar(2,sdthsmglp14outputM,stdsdthsmglp14outputM,'ko');
errorbar(3,sithsmlp14outputM,stdsithsmlp14outputM,'ko');
errorbar(4,sithsmglp14outputM,stdsithsmglp14outputM,'ko');
xlim([.5 4.5])

figure;
errorbar(1,hbmlp14outputM,stdhbmlp14outputM,'ko'); hold on;
errorbar(2,spmlp14outputM,stdspmlp14outputM,'ko');
errorbar(3,ssmlp14outputM,stdssmlp14outputM,'ko');
errorbar(4,sshsogmlp14outputM,stdsshsogmlp14outputM,'ko');
errorbar(5,sshthsmlp14outputM,stdsshsogmlp14outputM,'ko');
errorbar(6,sdsnamlp14outputM,stdsdsnamlp14outputM,'ko');
xlim([.5 6.5])

data = [ [sdsnamglp14.OutputM] [sdsnamglgp14.OutputM]  [sdsogmlp14.OutputM] [sdsogmglp14.OutputM] [sdsogmglgp14.OutputM] ...
       [sisnamlp14.OutputM] [sisnamglp14.OutputM] [sisogmlp14.OutputM] [sisogmglp14.OutputM]];
grp = [ 1*ones(1,length([sdsnamglp14.OutputM])) 2*ones(1,length([sdsnamglgp14.OutputM])) 3*ones(1,length([sdsogmlp14.OutputM])) 4*ones(1,length([sdsogmglp14.OutputM])) 5*ones(1,length([sdsogmglgp14.OutputM]))  ...
     6*ones(1,length([sisnamlp14.OutputM])) 7*ones(1,length([sisnamglp14.OutputM]))  8*ones(1,length([sisogmlp14.OutputM])) 9*ones(1,length([sisogmglp14.OutputM]))];

data = [ [sisnamlp14.OutputM] [sisnamglp14.OutputM] [sisogmlp14.OutputM] [sisogmglp14.OutputM]];
grp = [ 1*ones(1,length([sisnamlp14.OutputM])) 2*ones(1,length([sisnamglp14.OutputM])) 3*ones(1,length([sisogmlp14.OutputM])) 4*ones(1,length([sisogmglp14.OutputM]))];

figure;
boxplot(data,grp,'Labels',{'SI sna mlp', 'SI sna mglp', 'SI sog mlp', 'SI sog mglp'}); hold on
set(gca,'XTickLabelRotation',45)
ylabel('Output')

%% Active Nuclei
for i = 1:4
    totalsdsnamlp14(i) = length(sdsnamlp14(i).M(1,:));
    totalsdsnamglp14(i) = length(sdsnamglp14(i).M(1,:));
    totalsdsnamglgp14(i) = length(sdsnamglgp14(i).M(1,:));
    totalsdsogmlp14(i) = length(sdsogmlp14(i).M(1,:));
    totalsdsogmglp14(i) = length(sdsogmglp14(i).M(1,:));
    totalsdsogmglgp14(i) = length(sdsogmglgp14(i).M(1,:));
    totalsdthsmlp14(i) = length(sdthsmlp14(i).M(1,:));
    totalsdthsmglp14(i) = length(sdthsmglp14(i).M(1,:));
    totalsdthsmglgp14(i) = length(sdthsmglgp14(i).M(1,:));
    totalsisnamlp14(i) = length(sisnamlp14(i).M(1,:));
    totalsisnamglp14(i) = length(sisnamglp14(i).M(1,:));
    totalsisogmlp14(i) = length(sisogmlp14(i).M(1,:));
    totalsisogmglp14(i) = length(sisogmglp14(i).M(1,:));
    totalsithsmlp14(i) = length(sithsmlp14(i).M(1,:));
    totalsithsmglp14(i) = length(sithsmglp14(i).M(1,:));
    totalhbmlp14(i) = length(hbmlp14(i).M(1,:));
    totalspmlp14(i) = length(spmlp14(i).M(1,:));
    totalssmlp14(i) = length(ssmlp14(i).M(1,:));
    totalsshsogmlp14(i) = length(sshsogmlp14(i).M(1,:));
    totalsshthsmlp14(i) = length(sshthsmlp14(i).M(1,:));
    
    activenuc(1,i) = length([sdsnamlp14(i).active])/(totalsdsnamlp14(i)); % alternate variable name nucsdsnamlp
    activenuc(2,i) = length([sdsnamglp14(i).active])/(totalsdsnamglp14(i));
    activenuc(3,i) = length([sdsnamglgp14(i).active])/(totalsdsnamglgp14(i));
    activenuc(4,i) = length([sdsogmlp14(i).active])/(totalsdsogmlp14(i));
    activenuc(5,i) = length([sdsogmglp14(i).active])/(totalsdsogmglp14(i));
    activenuc(6,i) = length([sdsogmglgp14(i).active])/(totalsdsogmglgp14(i));
    activenuc(7,i) = length([sdthsmlp14(i).active])/(totalsdthsmlp14(i));
    activenuc(8,i) = length([sdthsmglp14(i).active])/(totalsdthsmglp14(i));
    activenuc(9,i) = length([sdthsmglgp14(i).active])/(totalsdthsmglgp14(i));
    activenuc(10,i) = length([sisnamlp14(i).active])/(totalsisnamlp14(i));
    activenuc(11,i) = length([sisnamglp14(i).active])/(totalsisnamglp14(i));
    activenuc(12,i) = length([sisogmlp14(i).active])/(totalsisogmlp14(i));
    activenuc(13,i) = length([sisogmglp14(i).active])/(totalsisogmglp14(i));
    activenuc(14,i) = length([sithsmlp14(i).active])/(totalsithsmlp14(i));
    activenuc(15,i) = length([sithsmglp14(i).active])/(totalsithsmglp14(i));
    
    activenuc(16,i) = length([hbmlp14(i).active])/(totalhbmlp14(i));
    activenuc(17,i) = length([spmlp14(i).active])/(totalspmlp14(i));
    activenuc(18,i) = length([ssmlp14(i).active])/(totalssmlp14(i));
    activenuc(19,i) = length([sshsogmlp14(i).active])/(totalsshsogmlp14(i));
    activenuc(20,i) = length([sshthsmlp14(i).active])/(totalsshthsmlp14(i));
end


%Combining all nuclie in one condition
% activenuc(1) = length([sdsnamlp14.active])/sum(totalsdsnamlp14); % alternate variable name nucsdsnamlp
% activenuc(2) = length([sdsnamglp14.active])/sum(totalsdsnamglp14);
% activenuc(3) = length([sdsnamglgp14.active])/sum(totalsdsnamglgp14);
% activenuc(4) = length([sdsogmlp14.active])/sum(totalsdsogmlp14);
% activenuc(5) = length([sdsogmglp14.active])/sum(totalsdsogmglp14);
% activenuc(6) = length([sdsogmglgp14.active])/sum(totalsdsogmglgp14);
% activenuc(7) = length([sdthsmlp14.active])/sum(totalsdthsmlp14);
% activenuc(8) = length([sdthsmglp14.active])/sum(totalsdthsmglp14);
% activenuc(9) = length([sdthsmglgp14.active])/sum(totalsdthsmglgp14);
% activenuc(10) = length([sisnamlp14.active])/sum(totalsisnamlp14);
% activenuc(11) = length([sisnamglp14.active])/sum(totalsisnamglp14);
% activenuc(12) = length([sisogmlp14.active])/sum(totalsisogmlp14);
% activenuc(13) = length([sisogmglp14.active])/sum(totalsisogmglp14);
% activenuc(14) = length([sithsmlp14.active])/sum(totalsithsmlp14);
% activenuc(15) = length([sithsmglp14.active])/sum(totalsithsmglp14);

meanactivenuc = mean(activenuc,2);
stdactivenuc = std(activenuc,[],2);


%elongation mean and std
elongationrates(1) = median(sdsnamlprate14,'omitnan');
elongationrates(2) = median(sdsnamglprate14,'omitnan');
elongationrates(3) = median(sdsnamglgprate14,'omitnan');
elongationrates(4) = median(sdsogmlprate14,'omitnan');
elongationrates(5) = median(sdsogmglprate14,'omitnan');
elongationrates(6) = median(sdsogmglgprate14,'omitnan');
elongationrates(7) = median(sdthsmlprate14,'omitnan');
elongationrates(8) = median(sdthsmglprate14,'omitnan');
elongationrates(9) = median(sdthsmglgprate14,'omitnan');
elongationrates(10) = median(sisnamlprate14,'omitnan');
elongationrates(11) = median(sisnamglprate14,'omitnan');
elongationrates(12) = median(sisogmlprate14,'omitnan');
elongationrates(13) = median(sisogmglprate14,'omitnan');
elongationrates(14) = median(sithsmlprate14,'omitnan');
elongationrates(15) = median(sithsmglprate14,'omitnan');

elongationrates(16) = median(hbrate14,'omitnan');
elongationrates(17) = median(sprate14,'omitnan');
elongationrates(18) = median(ssrate14,'omitnan');
elongationrates(19) = median(sshsograte14,'omitnan');
elongationrates(20) = median(sshthsrate14,'omitnan');

errelongationrates(1,:) = quantile(sdsnamlprate14,[.25 .75]);
errelongationrates(2,:) = quantile(sdsnamglprate14,[.25 .75]);
errelongationrates(3,:) = quantile(sdsnamglgprate14,[.25 .75]);
errelongationrates(4,:) = quantile(sdsogmlprate14,[.25 .75]);
errelongationrates(5,:) = quantile(sdsogmglprate14,[.25 .75]);
errelongationrates(6,:) = quantile(sdsogmglgprate14,[.25 .75]);
errelongationrates(7,:) = quantile(sdthsmlprate14,[.25 .75]);
errelongationrates(8,:) = quantile(sdthsmglprate14,[.25 .75]);
errelongationrates(9,:) = quantile(sdthsmglgprate14,[.25 .75]);
errelongationrates(10,:) = quantile(sisnamlprate14,[.25 .75]);
errelongationrates(11,:) = quantile(sisnamglprate14,[.25 .75]);
errelongationrates(12,:) = quantile(sisogmlprate14,[.25 .75]);
errelongationrates(13,:) = quantile(sisogmglprate14,[.25 .75]);
errelongationrates(14,:) = quantile(sithsmlprate14,[.25 .75]);
errelongationrates(15,:) = quantile(sithsmglprate14,[.25 .75]);

errelongationrates(16,:) = quantile(hbrate14,[.25 .75]);
errelongationrates(17,:) = quantile(sprate14,[.25 .75]);
errelongationrates(18,:) = quantile(ssrate14,[.25 .75]);
errelongationrates(19,:) = quantile(sshsograte14,[.25 .75]);
errelongationrates(20,:) = quantile(sshthsrate14,[.25 .75]);


figure;
yyaxis left
errorbar(1:15,meanactivenuc(1:15),stdactivenuc(1:15),'ko')
xticklabels({'SD sna mlp', 'SD sna mglp', 'SD sna mglgp', 'SD sog mlp','SD sog mglp','SD sog mglgp','SD ths mlp','SD ths mglp','SD ths mglgp','SI sna mlp' ...
    'SI sna mglp','SI sog mlp','SI sog mglp','SI ths mlp','SI ths mglp'})
xticks([1:15])
set(gca,'XTickLabelRotation',45)
ylabel('Percent Active Nuclei')
ylim([0 1])
xlim([0 15.5])

yyaxis right
errorbar(1:15,elongationrates(1:15),elongationrates(1:15) - errelongationrates([1:15],1)',elongationrates(1:15) - errelongationrates([1:15],2)','rs')

legend('Active Nuclei','Elongation Rates')


figure;
yyaxis left
errorbar(1:5,meanactivenuc(16:20),stdactivenuc(16:20),'ko')
yticklabels({'hb mlp','sna prox mlp','sna sh sna mlp','sna sh sog mlp','sna sh ths mlp','SD sna mlp','SD sna mglp','SD sna mglgp'})
xticks([1:15])
set(gca,'XTickLabelRotation',45)
ylabel('Percent Active Nuclei')
ylim([0 1])
xlim([0 5.5])

yyaxis right
errorbar(1:5,elongationrates(16:20),elongationrates(16:20) - errelongationrates([16:20],1)',elongationrates(16:20) - errelongationrates([16:20],2)','rs')

legend('Active Nuclei','Elongation Rates')

%% Average Signal Intensity
for j = 1:length(ssmlp14)
    for i = 1:length(ssmlp14(j).intT)
        avesignalssmlp(j,i) = mean(ssmlp14(j).intM(i,:));
        avesignalspmlp(j,i) = mean(spmlp14(j).intM(i,:));
    end
end

figure;
plot(1:100,mean(avesignalspmlp)), hold on
plot(1:100,mean(avesignalssmlp))
legend('sna prox','sna sh')

for i = 1:4
    maxhb(i) = max(max([hbmlp14(i).intM]));
    maxsp(i) = max(max([spmlp14(i).intM]));
    maxss(i) = max(max([ssmlp14(i).intM]));
    maxsdsnamlp(i) = max(max([sdsnamlp14(i).intM]));
end
maxhb = max(maxhb);
maxsp = max(maxsp);
maxss = max(maxss);
maxsdsnamlp = max(maxsdsnamlp);

for i = 1:4
    testoutputhb{i} = [hbmlp14(i).intM]./maxhb;
end

%% Calulate average trace
marker12 = sdsnamlp12;
marker13 = sdsnamlp13;
marker14 = sdsnamlp14;

for j = 1:length(marker12)
    for i = 1:length(marker12(j).intT)
        avesigM12(j,i) = mean(marker12(j).intM(i,:));
        avesigP12(j,i) = mean(marker12(j).intP(i,:));
    end
end

for j = 1:length(marker13)
    for i = 1:length(marker13(j).intT)
        avesigM13(j,i) = mean(marker13(j).intM(i,:));
        avesigP13(j,i) = mean(marker13(j).intP(i,:));
    end
end

for j = 1:length(marker14)
    for i = 1:length(marker14(j).intT)
        avesigM14(j,i) = mean(marker14(j).intM(i,:));
        avesigP14(j,i) = mean(marker14(j).intP(i,:));
    end
end
figure;
plot(1:100,mean(avesigM12)), hold on
plot(1:100,mean(avesigP12))
legend('MS2','PP7')
title('NC12')
ylim([0 6500])

figure;
plot(1:100,mean(avesigM13)), hold on
plot(1:100,mean(avesigP13))
legend('MS2','PP7')
title('NC13')
ylim([0 6500])

figure;
plot(1:100,mean(avesigM14)), hold on
plot(1:100,mean(avesigP14))
legend('MS2','PP7')
title('NC14')
ylim([0 6500])

%% MS2 PP7 Active Ratio
% sdsnaM12 = length([sdsnamlp12.activeM])/(95+122+102+113);
% sdsnaP12 = length([sdsnamlp12.activeP])/(95+122+102+113);
% sdsnaM13 = length([sdsnamlp13.activeM])/(182+248+108+195);
% sdsnaP13 = length([sdsnamlp13.activeP])/(182+248+108+195);
% sdsnaM14 = length([sdsnamlp14.activeM])/(356+555+328+461);
% sdsnaP14 = length([sdsnamlp14.activeP])/(356+555+328+461);
% figure; bar([sdsnaM12 sdsnaP12; sdsnaM13 sdsnaP13; sdsnaM14 sdsnaP14])
marker12 = sdsnamlp12;
marker13 = sdsnamlp13;
marker14 = sdsnamlp14;

for j = 1:length(marker12)
    activeM12(j) = length([marker12(j).activeM]);
    activeP12(j) = length([marker12(j).activeP]);
end
for j = 1:length(marker13)
    activeM13(j) = length([marker13(j).activeM]);
    activeP13(j) = length([marker13(j).activeP]);
end
for j = 1:length(marker14)
    activeM14(j) = length([marker14(j).activeM]);
    activeP14(j) = length([marker14(j).activeP]);
end

active12 = activeP12./activeM12;
active13 = activeP13./activeM13;
active14 = activeP14./activeM14;

figure; bar([mean(active12); mean(active13); mean(active14)]); hold on
errorbar(1,mean(active12),std(active12),'k');
errorbar(2,mean(active13),std(active13),'k');
errorbar(3,mean(active14),std(active14),'k');

%% Activation over time
for i = 1:4
    if ~isempty(sdsnamlp14(i).active)
        actsdsnamlp14(i,:) = interp1(sdsnamlp14(i).normT,sdsnamlp14(i).activetime./sdsnamlp14(i).activetime(end),[.01:.01:1]);
    else
       actsdsnamlp14(i,:) = NaN;
    end
    
    if ~isempty(sdsnamglp14(i).active)
        actsdsnamglp14(i,:) = interp1(sdsnamglp14(i).normT,sdsnamglp14(i).activetime./sdsnamglp14(i).activetime(end),[.01:.01:1]);
    else
       actsdsnamglp14(i,:) = NaN;
    end
    
    if ~isempty(sdsnamglgp14(i).active)
        actsdsnamglgp14(i,:) = interp1(sdsnamglgp14(i).normT,sdsnamglgp14(i).activetime./sdsnamglgp14(i).activetime(end),[.01:.01:1]);
    else
        actsdsnamglgp14(i,:) = NaN;
    end
    
    
    
    if ~isempty(sdsogmlp14(i).active)
        actsdsogmlp14(i,:) = interp1(sdsogmlp14(i).normT,sdsogmlp14(i).activetime./sdsogmlp14(i).activetime(end),[.01:.01:1]);
    else
        actsdsogmlp14(i,:) = NaN;
    end
    
    if ~isempty(sdsogmglp14(i).active)
        actsdsogmglp14(i,:) = interp1(sdsogmglp14(i).normT,sdsogmglp14(i).activetime./sdsogmglp14(i).activetime(end),[.01:.01:1]);
    else
        actsdsogmglp14(i,:) = NaN;
    end
    
    if ~isempty(sdsogmglgp14(i).active)
        actsdsogmglgp14(i,:) = interp1(sdsogmglgp14(i).normT,sdsogmglgp14(i).activetime./sdsogmglgp14(i).activetime(end),[.01:.01:1]);
    else
        actsdsogmglgp14(i,:) = NaN;
    end
    
    
    if ~isempty(sdthsmlp14(i).active)
        actsdthsmlp14(i,:) = interp1(sdthsmlp14(i).normT,sdthsmlp14(i).activetime./sdthsmlp14(i).activetime(end),[.01:.01:1]);
    else
        actsdthsmlp14(i,:) = NaN;
    end
    
    if ~isempty(sdthsmglp14(i).active)
        actsdthsmglp14(i,:) = interp1(sdthsmglp14(i).normT,sdthsmglp14(i).activetime./sdthsmglp14(i).activetime(end),[.01:.01:1]);
    else
        actsdthsmglp14(i,:) = NaN;
    end
    
    if ~isempty(sdthsmglgp14(i).active)
        actsdthsmglgp14(i,:) = interp1(sdthsmglgp14(i).normT,sdthsmglgp14(i).activetime./sdthsmglgp14(i).activetime(end),[.01:.01:1]);
    else
        actsdthsmglgp14(i,:) = NaN;
    end
    
    
    %SI
    if ~isempty(sisnamlp14(i).active)
        actsisnamlp14(i,:) = interp1(sisnamlp14(i).normT,sisnamlp14(i).activetime./sisnamlp14(i).activetime(end),[.01:.01:1]);
    else
        actsisnamlp14(i,:) = NaN;
    end
    
    if ~isempty(sisnamglp14(i).active)
        actsisnamglp14(i,:) = interp1(sisnamglp14(i).normT,sisnamglp14(i).activetime./sisnamglp14(i).activetime(end),[.01:.01:1]);
    else 
        actsisnamglp14(i,:) = NaN;
    end
        
    
    if ~isempty(sisogmlp14(i).active)
        actsisogmlp14(i,:) = interp1(sisogmlp14(i).normT,sisogmlp14(i).activetime./sisogmlp14(i).activetime(end),[.01:.01:1]);
    else
        actsisogmlp14(i,:) = NaN;
    end
    
    if ~isempty(sisogmglp14(i).active)
        actsisogmglp14(i,:) = interp1(sisogmglp14(i).normT,sisogmglp14(i).activetime./sisogmglp14(i).activetime(end),[.01:.01:1]);
    else
        actsisogmglp14(i,:) = NaN;
    end
    
    
    if ~isempty(sithsmlp14(i).active)
        actsithsmlp14(i,:) = interp1(sithsmlp14(i).normT,sithsmlp14(i).activetime./sithsmlp14(i).activetime(end),[.01:.01:1]);
    else
        actsithsmlp14(i,:) = NaN;
    end
    
    if ~isempty(sithsmglp14(i).active)
        actsithsmglp14(i,:) = interp1(sithsmglp14(i).normT,sithsmglp14(i).activetime./sithsmglp14(i).activetime(end),[.01:.01:1]);
    else
        actsithsmglp14(i,:) = NaN;
    end
    
    
    if ~isempty(hbmlp14(i).active)
        acthbmlp14(i,:) = interp1(hbmlp14(i).normT,hbmlp14(i).activetime./hbmlp14(i).activetime(end),[.01:.01:1]);
    else
       acthbmlp14(i,:) = NaN;
    end
    
    if ~isempty(spmlp14(i).active)
        actspmlp14(i,:) = interp1(spmlp14(i).normT,spmlp14(i).activetime./spmlp14(i).activetime(end),[.01:.01:1]);
    else
       actspmlp14(i,:) = NaN;
    end
    
    if ~isempty(ssmlp14(i).active)
        actssmlp14(i,:) = interp1(ssmlp14(i).normT,ssmlp14(i).activetime./ssmlp14(i).activetime(end),[.01:.01:1]);
    else
        actssmlp14(i,:) = NaN;
    end
    
    if ~isempty(sshsogmlp14(i).active)
        actsshsogmlp14(i,:) = interp1(sshsogmlp14(i).normT,sshsogmlp14(i).activetime./sshsogmlp14(i).activetime(end),[.01:.01:1]);
    else
        actsshsogmlp14(i,:) = NaN;
    end
    
    if ~isempty(sshthsmlp14(i).active)
        actsshthsmlp14(i,:) = interp1(sshthsmlp14(i).normT,sshthsmlp14(i).activetime./sshthsmlp14(i).activetime(end),[.01:.01:1]);
    else
        actsshthsmlp14(i,:) = NaN;
    end
end

%take average and std of activity
aveactsdsnamlp14 = nanmean(actsdsnamlp14);
aveactsdsnamglp14 = nanmean(actsdsnamglp14);
aveactsdsnamglgp14 = nanmean(actsdsnamglgp14);

aveactsdsogmlp14 = nanmean(actsdsogmlp14);
aveactsdsogmglp14 = nanmean(actsdsogmglp14);
aveactsdsogmglgp14 = nanmean(actsdsogmglgp14);

aveactsdthsmlp14 = nanmean(actsdthsmlp14);
aveactsdthsmglp14 = nanmean(actsdthsmglp14);
aveactsdthsmglgp14 = nanmean(actsdthsmglgp14);

aveactsisnamlp14 = nanmean(actsisnamlp14);
aveactsisnamglp14 = nanmean(actsisnamglp14);

aveactsisogmlp14 = nanmean(actsisogmlp14);
aveactsisogmglp14 = nanmean(actsisogmglp14);

aveactsithsmlp14 = nanmean(actsithsmlp14);
aveactsithsmglp14 = nanmean(actsithsmglp14);


aveacthbmlp14 = nanmean(acthbmlp14);
aveactspmlp14 = nanmean(actspmlp14);
aveactssmlp14 = nanmean(actssmlp14);
aveactsshsogmlp14 = nanmean(actsshsogmlp14);
aveactsshthsmlp14 = nanmean(actsshthsmlp14);


stdactsdsnamlp14 = nanstd(actsdsnamlp14);
stdactsdsnamglp14 = nanstd(actsdsnamglp14);
stdactsdsnamglgp14 = nanstd(actsdsnamglgp14);

stdactsdsogmlp14 = nanstd(actsdsogmlp14);
stdactsdsogmglp14 = nanstd(actsdsogmglp14);
stdactsdsogmglgp14 = nanstd(actsdsogmglgp14);

stdactsdthsmlp14 = nanstd(actsdthsmlp14);
stdactsdthsmglp14 = nanstd(actsdthsmglp14);
stdactsdthsmglgp14 = nanstd(actsdthsmglgp14);

stdactsisnamlp14 = nanstd(actsisnamlp14);
stdactsisnamglp14 = nanstd(actsisnamglp14);

stdactsisogmlp14 = nanstd(actsisogmlp14);
stdactsisogmglp14 = nanstd(actsisogmglp14);

stdactsithsmlp14 = nanstd(actsithsmlp14);
stdactsithsmglp14 = nanstd(actsithsmglp14);

stdacthbmlp14 = nanstd(acthbmlp14);
stdactspmlp14 = nanstd(actspmlp14);
stdactssmlp14 = nanstd(actssmlp14);
stdactsshsogmlp14 = nanstd(actsshsogmlp14);
stdactsshthsmlp14 = nanstd(actsshthsmlp14);

figure;
shadedErrorBar([.01:.01:1],aveactsdsogmlp14,stdactsdsogmlp14,{'color',rgb('black'),'LineWidth',2},1); hold on
shadedErrorBar([.01:.01:1],aveactsdsogmglp14,stdactsdsogmglp14,{'color',rgb('red'),'LineWidth',2},1); 
shadedErrorBar([.01:.01:1],aveactsdsogmglgp14,stdactsdsogmglgp14,{'color',rgb('blue'),'LineWidth',2},1); 
legend('SD sog mlp','SD sog mglp','SD sog mglgp')
ylim([0 1])

figure;
shadedErrorBar([.01:.01:1],aveactsisnamglp14,stdactsisnamglp14,{'color',rgb('black'),'LineWidth',2},1); hold on
shadedErrorBar([.01:.01:1],aveactsisogmglp14,stdactsisogmglp14,{'color',rgb('red'),'LineWidth',2},1); 
shadedErrorBar([.01:.01:1],aveactsithsmglp14,stdactsithsmglp14,{'color',rgb('blue'),'LineWidth',2},1); 
legend('SD sna mlp','SD sog mlp','SD ths mlp')
ylim([0 1])

figure;
shadedErrorBar([.01:.01:1],aveactsdthsmglp14,stdactsdthsmglp14,{'color',rgb('black'),'LineWidth',2},1); hold on
shadedErrorBar([.01:.01:1],aveactsithsmglp14,stdactsithsmglp14,{'color',rgb('red'),'LineWidth',2},1); 
legend('SD ths mglp','SI ths mglp')
ylim([0 1])

figure;
shadedErrorBar([.01:.01:1],aveacthbmlp14,stdacthbmlp14,{'color',rgb('black'),'LineWidth',2},1); hold on
shadedErrorBar([.01:.01:1],aveactspmlp14,stdactspmlp14,{'color',rgb('red'),'LineWidth',2},1); hold on
shadedErrorBar([.01:.01:1],aveactssmlp14,stdactssmlp14,{'color',rgb('blue'),'LineWidth',2},1); hold on
shadedErrorBar([.01:.01:1],aveactsdsnamlp14,stdactsdsnamlp14,{'color',rgb('green'),'LineWidth',2},1); hold on
% shadedErrorBar([.01:.01:1],aveactsisnamlp14,stdactsisnamlp14,{'color',rgb('orange'),'LineWidth',2},1); hold on
legend('hb mlp','sna prox mlp','sna sh sna mlp','sd sna mlp')
ylim([0 1])

%% False Color Embryos

path = 'G:\Shared drives\Lim_Lab\Sam\Elongation Project\sna_sh_ms2_lacz_pp7\ssmlpn';


% load([path 'trajectories_ms2.mat']);

marker = ssmlp14;
for q = 1:4
    pathfull = sprintf('%d',q);
    fullpath = [path,pathfull,'\'];
    load([fullpath 'segmentation_lineage14.mat'], 'convex_image','nuc_lineage','images_seg','lineage_cx','lineage_cy','images_ms2', ...
        'images_nuc');
    
for i=round(size(marker(q).M,1)/2) %middle time point
% for i = 50:50
    F = zeros(size(images_seg,1),size(images_seg,2),3,'uint16');
    
    for j=1:size(marker(q).M,2) % all nuclei
        c_image = convex_image{i,nuc_lineage(i,(j))};
        c_image = imresize(c_image,0.8);
        r1 = regionprops(c_image,'PixelList','centroid');
        
        dx = lineage_cx(i,(j))-r1.Centroid(1); dx = round(dx);
        dy = lineage_cy(i,(j))-r1.Centroid(2); dy = round(dy);
        r1.PixelList(:,1) = r1.PixelList(:,1)+dx;
        r1.PixelList(:,2) = r1.PixelList(:,2)+dy;

        ddum = find(r1.PixelList(:,1)<1);
        r1.PixelList(ddum,1) = 1;
        ddum= find(r1.PixelList(:,2)<1);
        r1.PixelList(ddum,2)=1;
        ddum = find(r1.PixelList(:,1)>size(images_ms2,2));
        r1.PixelList(ddum,1) = size(images_ms2,2);
        ddum = find(r1.PixelList(:,2)>size(images_seg,1));
        r1.PixelList(ddum,2) = size(images_ms2,1);
        
        xx = r1.PixelList(:,1);
        yy = r1.PixelList(:,2);
        
        for k=1:length(xx)
            F(yy(k),xx(k),2) = 10000;
%             if M1(i,j) > 100
%                 F(yy(k),xx(k),1) = 35000;
%             end
            if ismember(j,marker(q).active)
                F(yy(k),xx(k),1) = 70000;%M1(i,j)*10;
            end
%             F(yy(k),xx(k),1) = outM1(j);
        end        
    end
    F(:,:,3) = images_nuc(:,:,i)*12;
%         imwrite(F,sprintf('%s%03d.png',[path 'falsecolor\'],i));
end
figure,imshow(F)

clear c_image convex_image nuc_lineage images_seg
end

%% Output Heatmap

% path = 'G:\Shared drives\Lim_Lab\Sam\Elongation Project\hb_ms2_lacz_pp7\hbmlpn'; %n5
% path = 'G:\Shared drives\Lim_Lab\Sam\Elongation Project\sna_prox_ms2_lacz_pp7\spmlpn'; %n2
% path = 'G:\Shared drives\Lim_Lab\Sam\Elongation Project\sna_sh_ms2_lacz_pp7\ssmlpn'; %n4
path = 'G:\Shared drives\Lim_Lab\Sam\Elongation Project\sogdist_sna_MLP\sdsnamlpn'; %n2

q = 2; %embryo number

% load([path 'trajectories_ms2.mat']);

marker = sdsnamlp14;

dX = marker(q).intT(2)-marker(q).intT(1);

% %Output for all nuclei
% for i = 1:size(marker(q).M,2)
%     Output(i) = trapz(marker(q).intT,marker(q).intM(:,i))/(dX*100);
% end
    
%Output looking only at first 10 min
dX = marker(q).intT(2)-marker(q).intT(1);
%find where intT is 10 min
min5 = find(marker(q).intT > 10);
for i = 1:size(marker(q).M,2)
    Output(i) = trapz(marker(q).intT(1:min5(1)),marker(q).intM(1:min5(1),i))/(dX*min5(1));
end


% pathfull = sprintf('%d',q);
fullpath = [path,'2\'];
load([fullpath 'segmentation_lineage14.mat'], 'convex_image','nuc_lineage','images_seg','lineage_cx','lineage_cy','images_ms2', ...
    'images_nuc');
    
for i=round(size(marker(q).M,1)/2) %middle time point
% for i = 50:50
    F = zeros(size(images_seg,1),size(images_seg,2));
    
    for j=1:size(marker(q).M,2) % all nuclei
        c_image = convex_image{i,nuc_lineage(i,(j))};
        c_image = imresize(c_image,0.85);
        r1 = regionprops(c_image,'PixelList','centroid');
        
        dx = lineage_cx(i,(j))-r1.Centroid(1); dx = round(dx);
        dy = lineage_cy(i,(j))-r1.Centroid(2); dy = round(dy);
        r1.PixelList(:,1) = r1.PixelList(:,1)+dx;
        r1.PixelList(:,2) = r1.PixelList(:,2)+dy;

        ddum = find(r1.PixelList(:,1)<1);
        r1.PixelList(ddum,1) = 1;
        ddum= find(r1.PixelList(:,2)<1);
        r1.PixelList(ddum,2)=1;
        ddum = find(r1.PixelList(:,1)>size(images_ms2,2));
        r1.PixelList(ddum,1) = size(images_ms2,2);
        ddum = find(r1.PixelList(:,2)>size(images_seg,1));
        r1.PixelList(ddum,2) = size(images_ms2,1);
        
        xx = r1.PixelList(:,1);
        yy = r1.PixelList(:,2);
        
        for k=1:length(xx)
            F(yy(k),xx(k)) = Output(j);
            
%             if M1(i,j) > 100
%                 F(yy(k),xx(k),1) = 35000;
%             end
%             if ismember(j,marker(q).active)
%                 F(yy(k),xx(k),1) = 70000;%M1(i,j)*10;
%             end
%             F(yy(k),xx(k),1) = outM1(j);
        end        
    end
%     F(:,:,3) = images_nuc(:,:,i)*12;
%         imwrite(F,sprintf('%s%03d.png',[path 'falsecolor\'],i));
end

%create transparancy
transp = zeros(length(F));
transp(F ~= 0) = 1;
figure;
ax1 = axes;
imshow(imread([fullpath 'sdsnamlpn2 2proj.tif'],380)*15); hold on; %hbmlpn5 332; spmlpn2 269; ssmlpn4 335 sdsnamlpn2 380
colormap(ax1,'gray');
ax2 = axes;
imagesc(F,'Alphadata',transp);
colormap(ax2,parula)
caxis(ax2,[min(nonzeros(F)) max(nonzeros(F))]);
ax2.Visible = 'off';
linkprop([ax1,ax2],'Position');
colorbar
caxis([0 8e3])
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);

clear c_image convex_image nuc_lineage images_seg


%% ANOVA and multi comparison test
% marker = sisnamlp14;
% marker2 = sisnamlprate14;
% for i = 1:4
%     split(i) = length(marker(i).active);
% end

% Cell = {marker2(1:split(1));marker2(split(1)+1:sum(split(1:2)));marker2(sum(split(1:2))+1:sum(split(1:3)));marker2(sum(split(1:3))+1:end)};
% Cell = {[hbmlp14.OutputWindowM];[spmlp14.OutputWindowM]; [ssmlp14.OutputWindowM];[sdsnamlp14.OutputWindowM]};
% Cell = {sdsnamlprate14; sdsogmlprate14; sdthsmlprate14};
Cell = {sdsogmglprate14; sdsogmyeprate14; sdsogmmweprate14; sdsogmyprate14; sdsogmmwprate14};

% Cell = {hbrate14;sprate14;ssrate14;sdsnamlprate14};
dataCell = [];
grp = [];



for i = 1:length(Cell)
    dataCell=[dataCell, Cell{i}];
    grp=[grp, ones(1, length(Cell{i}))*i];
end

[p,tbl,stats] = anova1(dataCell,grp);
xticklabels({'hb','sna proximal','sna distal','SD sna'})
% set(gca,'XTickLabelRotation',45)
ylabel('Elongation Rate (kb/min)')

figure; c = multcompare(stats,'alpha',0.01)
yticklabels(flip({'mglp','myp','myep','mw','mwe','ke'}))


%% Proving xcorr vs leading edge
edges = linspace(0,10,50);
figure; histogram(sdsogmyprate14,edges); hold on
histogram(sdsogmyeprate14,edges);

%% Transcription Acceleration
hbacc = ([hbmlp14.xcorr] - hbrate14);%./hbrate14;
spacc = ([spmlp14.xcorr] - sprate14);%./sprate14;
ssacc = ([ssmlp14.xcorr] - ssrate14);%./ssrate14;
sshsogacc = ([sshsogmlp14.xcorr] - sshsograte14);%./sshsograte14;
sshthsacc = ([sshthsmlp14.xcorr] - sshthsrate14);%./sshthsrate14;
sdsnamlpacc = ([sdsnamlp14.xcorr] - sdsnamlprate14);%./sdsnamlprate14;
sdsnamglpacc = ([sdsnamglp14.xcorr] - sdsnamglprate14);%./sdsnamglprate14;
sdsnamglgpacc = ([sdsnamglgp14.xcorr] - sdsnamglgprate14);%./sdsnamglgprate14;
sdsogmlpacc = ([sdsogmlp14.xcorr] - sdsogmlprate14);%./sdsogmlprate14;
sdsogmglpacc = ([sdsogmglp14.xcorr] - sdsogmglprate14);%./sdsogmglprate14;
sdsogmglgpacc = ([sdsogmglgp14.xcorr] - sdsogmglgprate14);%./sdsogmglgprate14;
sdthsmlpacc = ([sdthsmlp14.xcorr] - sdthsmlprate14);%./sdthsmlprate14;
sdthsmglpacc = ([sdthsmglp14.xcorr] - sdthsmglprate14);%./sdthsmglprate14;
sdthsmglgpacc = ([sdthsmglgp14.xcorr] - sdthsmglgprate14);%./sdthsmglgprate14;
sisnamlpacc = ([sisnamlp14.xcorr] - sisnamlprate14);%./sisnamlprate14;
sisnamglpacc = ([sisnamglp14.xcorr] - sisnamglprate14);%./sisnamglprate14;
sisogmlpacc = ([sisogmlp14.xcorr] - sisogmlprate14);%./sisogmlprate14;
sisogmglpacc = ([sisogmglp14.xcorr] - sisogmglprate14);%./sisogmglprate14;
sithsmlpacc = ([sithsmlp14.xcorr] - sithsmlprate14);%./sithsmlprate14;
sithsmglpacc = ([sithsmglp14.xcorr] - sithsmglprate14);%./sithsmglprate14;
hbmypacc = ([hbmyp14.xcorr] - hbmyprate14);%./hbmyprate14;
sdsogmypacc = ([sdsogmyp14.xcorr] - sdsogmyprate14);%./sdsogmyprate14;
sdsogmyepacc = ([sdsogmyep14.xcorr] - sdsogmyeprate14);%./sdsogmyeprate14;
sdsogmmwpacc = ([sdsogmmwp14.xcorr] - sdsogmmwprate14);%./sdsogmmwprate14;
sdsogmmwepacc = ([sdsogmmwep14.xcorr] - sdsogmmweprate14);%./sdsogmmweprate14;
shsmyepacc = ([shsmyep14.xcorr] - shsmyeprate14);%./shsmyeprate14;
sdsogmkepacc = ([sdsogmkep14.xcorr] - sdsogmkeprate14);%./sdsogmkeprate14;


Cell = {hbacc; spacc;ssacc;sdsnamlpacc;sisnamlpacc};
dataCell = [];
grp = [];

for i = 1:length(Cell)
    dataCell=[dataCell, Cell{i}];
    grp=[grp, ones(1, length(Cell{i}))*i];
end
figure;
boxplot(dataCell,grp,'notch','on'); hold on
xticklabels({'hb mlp','sna prox mlp','sna sh sna mlp','SD sna mlp','SI sna mlp'})
set(gca,'XTickLabelRotation',45)
ylabel('Pol II Acceleration')
hold on; vline(0,'k:');
ylim([-1 2])


%get 10 - 90 percentile

% tmp = prctile(sdsogmmwpacc,[10,90]);
% sdsogmmwpacc(sdsogmmwpacc < tmp(1)) = [];
% sdsogmmwpacc(sdsogmmwpacc > tmp(2)) = [];


%Try violin plot
figure; violin(1,hbacc,'rotation','horizontal','withmdn',1); hold on
violin(2,spacc,'rotation','horizontal','withmdn',1);
violin(3,ssacc,'rotation','horizontal','withmdn',1);
violin(4,sdsnamlpacc,'rotation','horizontal','withmdn',1);
violin(5,sdsogmyepacc,'rotation','horizontal','withmdn',1);
violin(6,sdsogmypacc,'rotation','horizontal','withmdn',1);
violin(7,sdsogmmwepacc,'rotation','horizontal','withmdn',1);
violin(8,sdsogmmwpacc,'rotation','horizontal','withmdn',1);
violin(9,sdsogmlpacc,'rotation','horizontal','withmdn',1);
violin(10,sdsogmglpacc,'rotation','horizontal','withmdn',1);
violin(11,sdsogmglgpacc,'rotation','horizontal','withmdn',1);
% violin(5,sisnamlpacc,'rotation','horizontal','withmdn',1);
% violin(10,shsmyepacc,'rotation','horizontal','withmdn',1);
% violin(11,sdsogmkepacc,'rotation','horizontal','withmdn',1);

vline(0,'r')
xlim([-8 8])
set(gca,'XMinorTIck','on')
yticks([1:11])
xlabel('Elongation Difference (kb/min)')
yticklabels({'hb mlp','sna prox mlp','sna sh sna mlp','sd sna mlp','sd sog my(e)p','sd sog myp','sd sog mmw(e)p','sd sog mmwp'})
