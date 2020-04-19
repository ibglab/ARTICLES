function []= fig1BCD ()

SP=[2,5,1];
counter=1;
figure; hold on;

% Sampling rate or rat data
Fs=44000;

% Medium spiny neuron (MSN)
type='MSN';
sortedTxtFile='F130623-0012.matCRAW_013FiltFilt.txt';
rawBinFile='F130623-0012.matCRAW_013.bin';
rawFiltFiltBinFile='F130623-0012.matCRAW_013FiltFilt.bin';
rawFilterBinFile='F130623-0012.matCRAW_013Filter.bin';
[spkTime, rawSig,LPsig, NLPsig]=openFiles (sortedTxtFile,rawBinFile,rawFiltFiltBinFile, rawFilterBinFile,Fs);
[dLP(counter), dNLP(counter), snrLP(counter), snrNLP(counter)] = meanWF (spkTime, rawSig,LPsig, NLPsig,Fs, type, SP);
counter=counter+1;
SP(3)=counter;

% Tonically active neuron (TAN)
type='TAN';
sortedTxtFile='F130702-0007.matCRAW_004FiltFilt.txt';
rawBinFile='F130702-0007.matCRAW_004.bin';
rawFiltFiltBinFile='F130702-0007.matCRAW_004FiltFilt.bin';
rawFilterBinFile='F130702-0007.matCRAW_004Filter.bin';
[spkTime, rawSig,LPsig, NLPsig]=openFiles (sortedTxtFile,rawBinFile,rawFiltFiltBinFile, rawFilterBinFile,Fs);
[dLP(counter), dNLP(counter), snrLP(counter), snrNLP(counter)] = meanWF (spkTime, rawSig,LPsig, NLPsig,Fs, type, SP);
counter=counter+1;
SP(3)=counter;

% Fast spiking interneuron (FSI)
type='FSI';
sortedTxtFile='F130623-0012.matCRAW_004FiltFilt.txt';
rawBinFile='F130623-0012.matCRAW_004.bin';
rawFiltFiltBinFile='F130623-0012.matCRAW_004FiltFilt.bin';
rawFilterBinFile='F130623-0012.matCRAW_004Filter.bin';
[spkTime, rawSig,LPsig, NLPsig]=openFiles (sortedTxtFile,rawBinFile,rawFiltFiltBinFile, rawFilterBinFile,Fs);
[dLP(counter), dNLP(counter), snrLP(counter), snrNLP(counter)] = meanWF (spkTime, rawSig,LPsig, NLPsig,Fs, type, SP);
counter=counter+1;
SP(3)=counter;

%primates data y axis: AU.
Fs=40000;

% Low frequency burster (LFB)
type='LFB';
sortedTxtFile='A0801068.matCSPK16FiltFilt.txt';
rawBinFile='A0801068.matCSPK16.bin';
rawFiltFiltBinFile='A0801068.matCSPK16FiltFilt.bin';
rawFilterBinFile='A0801068.matCSPK16Filter.bin';
[spkTime, rawSig,LPsig, NLPsig]=openFiles (sortedTxtFile,rawBinFile,rawFiltFiltBinFile, rawFilterBinFile,Fs);
[dLP(counter), dNLP(counter), snrLP(counter), snrNLP(counter)] = meanWF (spkTime, rawSig,LPsig, NLPsig,Fs, type, SP);
counter=counter+1;
SP(3)=counter;

% High frequency pauser (HFP)
type='HFP';
sortedTxtFile='A0801064.matCSPK11FiltFilt.txt';
rawBinFile='A0801064.matCSPK11.bin';
rawFiltFiltBinFile='A0801064.matCSPK11FiltFilt.bin';
rawFilterBinFile='A0801064.matCSPK11Filter.bin';
[spkTime, rawSig,LPsig, NLPsig]=openFiles (sortedTxtFile,rawBinFile,rawFiltFiltBinFile, rawFilterBinFile,Fs);
[dLP(counter), dNLP(counter), snrLP(counter), snrNLP(counter)] = meanWF (spkTime, rawSig,LPsig, NLPsig,Fs, type, SP);
counter=counter+1;
SP(3)=counter;

subplot (2,5,[6:7]);
bar ([dLP;dNLP]');
ylabel ('D');
x=({'MSN','TAN','FSI','LFB','HFP'});
set(gca,'XTickLabel',x);

subplot (2,5,[9:10]);
bar ([snrNLP./snrLP]);
ylabel ('%SNR');
v=axis;
axis([v(1) v(2) 0 1]);
set(gca,'XTickLabel',x);

end


function [spkTime, rawSig,LPsig, NLPsig] = openFiles (sortedTxtFile,rawBinFile,rawFiltFiltBinFile, rawFilterBinFile,Fs)

micVCorr = 200/38.147;          % Correction to microVolt units
fid=fopen (rawBinFile);
rawSig=fread(fid,'int16');
rawSig=rawSig/micVCorr;
fclose (fid);
fid=fopen (rawFiltFiltBinFile);
LPsig=fread(fid,'int16');
LPsig=LPsig/micVCorr;
fclose (fid);
fid=fopen (rawFilterBinFile);
NLPsig=fread(fid,'int16');
NLPsig=NLPsig/micVCorr;
fclose (fid);

%spkTime
x=csvread(sortedTxtFile,1);
spkTime=x(find(x(:,2)==1),3);
spkTime=round((spkTime)*Fs);

end

function [dLP, dNLP, snrLP, snrNLP] = meanWF (spkTime, rawSig,LPsig, NLPsig,Fs, type, SP) %SP=subplot

timeWinBefore=0.5*Fs/1000;
timeWinAfter=1*Fs/1000;

rawWF=zeros(length(spkTime),timeWinBefore+timeWinAfter+1);
LPWF=zeros(length(spkTime),timeWinBefore+timeWinAfter+1);
NLPWF=zeros(length(spkTime),timeWinBefore+timeWinAfter+1);

for i=1:length (spkTime)
    rawWF(i,:)=rawSig((spkTime(i)-timeWinBefore):(spkTime(i)+timeWinAfter));
    LPWF(i,:)=LPsig((spkTime(i)-timeWinBefore):(spkTime(i)+timeWinAfter));
    NLPWF(i,:)=NLPsig((spkTime(i)-timeWinBefore):(spkTime(i)+timeWinAfter));
end


meanRawWF=(mean(rawWF));
meanRawWF=meanRawWF-meanRawWF(10);
meanLPWF=(mean(LPWF));
meanLPWF=meanLPWF-meanLPWF(10);
meanNLPWF=mean(NLPWF);

%align WFs
if (type =='TAN')
    [vRaw,iRaw]=max(meanRawWF);
    [vNLP,iNLP]=max(meanNLPWF);
    [vLP,iLP]=max(meanLPWF);
else
    [vRaw,iRaw]=min(meanRawWF);
    [vNLP,iNLP]=min(meanNLPWF);
    [vLP,iLP]=min(meanLPWF);
end
diffNLPWF=iNLP-iRaw;
diffLPWF=iLP-iRaw;
if (type =='MSN')
    diffLPWF=0;
elseif (type =='FSI')
    diffLPWF=0;
end

for i=1:length (spkTime)
    NLPWF(i,:)=NLPsig((spkTime(i)+diffNLPWF-timeWinBefore):(spkTime(i)+diffNLPWF+timeWinAfter));
    LPWF(i,:)=LPsig((spkTime(i)+diffLPWF-timeWinBefore):(spkTime(i)+diffLPWF+timeWinAfter));
end

meanNLPWF=mean(NLPWF);
meanLPWF=mean(LPWF);
meanLPWF=meanLPWF-meanLPWF(10);

subplot (SP(1), SP(2), SP(3)); hold on;
plot([0:size(meanRawWF,2)-1]/(Fs/1000),meanRawWF,'k');
plot([0:size(meanLPWF,2)-1]/(Fs/1000),(meanLPWF),'r');
plot([0:size(meanNLPWF,2)-1]/(Fs/1000),(meanNLPWF),'b');


title (type);
xlabel('Time (ms)');
ylabel ('V (µV)');


% Distance
dLP=sqrt(sum((meanLPWF-meanRawWF).^2))./sqrt(sum((meanRawWF).^2));
dNLP=sqrt(sum((meanNLPWF-meanRawWF).^2))./sqrt(sum((meanRawWF).^2));

% Signal to noise
if (type =='TAN')
    snrLP=(max(meanLPWF)./std(LPsig));
    snrNLP=(max(meanNLPWF)./std(NLPsig));
else
    snrLP=abs(min(meanLPWF)./std(LPsig));
    snrNLP=abs(min(meanNLPWF)./std(NLPsig));
end

end

