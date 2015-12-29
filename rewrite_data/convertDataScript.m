close all; clc;
clear all;

load svmdata;
uwSvmWrite('svm.dat', svm.SupportVectors, svm.Alpha, svm.Bias, svm.ScaleData.shift, svm.ScaleData.scaleFactor, b);

load scriptdata;
for i=1:8
    baseName = Filenames(i,1:end-8);
    baseName = ['Prm_', baseName];
    baseName = ['./', Pathnames(i,:), baseName];
    eval(['load ', baseName]);
    baseWriteName = [baseName, '.dat'];
    uwBaseWrite(baseWriteName, Ang_endpt, MagMax, MagMin, fft_dR);
    disp('writing base done.');
    
    for j=10001:10000+tns(i)
        frameName = ['./', Pathnames(i,:), Filenames(i,1:end-9)];
        frameName = sprintf('%s%d', frameName, j);
        frameWriteName = [frameName, '.dat'];
        eval(['load ', frameName]);
        uwEachWrite(frameWriteName, Ang_plot, RP_dBm);
        disp(['writing frame ', frameName, ' done.']);
    end
end
