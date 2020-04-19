function []= figure5BC ()

sortedTxtFile='F130623-0012.matCRAW_013FiltFilt.txt';
rawBinFile='F130623-0012.matCRAW_013.bin';
rawFilterBinFile='F130623-0012.matCRAW_013Filter.bin';
Fs=44000;

[spkTime, rawSig, NLPsig]=openFiles (sortedTxtFile,rawBinFile, rawFilterBinFile,Fs);

[rawWF]=meanWF (spkTime, rawSig,Fs);
[NLPWF] = meanWF (spkTime, NLPsig,Fs);


%IBG%%%%%%%%%%%%%%%%%%%%%%

opt = 4;
l=size(NLPWF,2);
for k=1:size(NLPWF,1)
    switch (opt)
        case 1 %linear decay
            sl = (NLPWF(k,l)-NLPWF(k,l-2))/2;
            NLPWF(k,l+[1:10]) = NLPWF(k,l)+sl*[1:10];
        case 2 %constant value
            NLPWF(k,l+[1:10]) = NLPWF(k,l);
        case 3 %zero padding
            NLPWF(k,l+[1:10]) = 0;
        case 4%spline fit 
            NLPWF(k,l+[1:10])  = interp1([1:l l+11:l+20],[NLPWF(k,1:l) zeros(1,10)],[l+1:l+10],'spline');
    end
end

%IBG%%%%%%%%%%%%%%%%%%%%%%

%correct concatenated segments
concatSeg=reshape (NLPWF', ([1, size(NLPWF,1)*size(NLPWF,2)]));
correctedSig= correction (concatSeg, Fs);
[correctedWF] = (reshape(correctedSig,[size(NLPWF,2),length(spkTime)]))';

%correct segmented waveforms matrix
[correctedWFmat]=correctedMat (NLPWF, Fs);

%elongation (padding) removal
correctedWF(:,(l+1:end))=[];
correctedWFmat(:,(l+1:end))=[];

%mean WFs
meanRawWF=mean(rawWF);
meanNLPWF=mean(NLPWF);
meanCWF=mean(correctedWF);
meanCorrectedWF=meanCWF-meanCWF(5);
meanCorrectedMatWF=mean(correctedWFmat);
meanCorrectedMatWF=meanCorrectedMatWF-meanCorrectedMatWF(5);

%plot
figure; hold on;
plot([0:size(meanRawWF,2)-1]/(Fs/1000),meanRawWF,'k');
plot([0:size(meanNLPWF,2)-1]/(Fs/1000),(meanNLPWF),'b');
plot([0:size(meanCorrectedWF,2)-1]/(Fs/1000),(meanCorrectedWF),'g');
plot([0:size(meanCorrectedMatWF,2)-1]/(Fs/1000),(meanCorrectedMatWF),'m');



function [spkTime, rawSig, NLPsig] = openFiles (sortedTxtFile,rawBinFile, rawFilterBinFile,Fs)

fid=fopen (rawBinFile);
rawSig=fread(fid,'int16');
rawSig=(rawSig*38.147)/200;
fclose (fid);
fid=fopen (rawFilterBinFile);
NLPsig=fread(fid,'int16');
NLPsig=(NLPsig*38.147)/200;
fclose (fid);

%spkTime
x=csvread(sortedTxtFile,1);
spkTime=x(find(x(:,2)==1),3);
spkTime=round((spkTime)*Fs);




function [WFmat] = meanWF (spkTime, signal, Fs)

timeWinBefore=round(0.3*Fs/1000);
timeWinAfter=round(1.1*Fs/1000);

if length(spkTime)==1
   spkTime= [0:(spkTime)-1]*(timeWinBefore+timeWinAfter+1)+1; 
   timeWinAfter=timeWinAfter+timeWinBefore;
   timeWinBefore=0;

end

WFmat=zeros(length(spkTime),timeWinBefore+timeWinAfter+1);

for i=1:length (spkTime)
    WFmat(i,:)=signal((spkTime(i)-timeWinBefore):(spkTime(i)+timeWinAfter));
end




function [correctedSig]= correction (sig, Fs)

HP=300;
LP=6000;
[b, a] = butter(4, [HP LP]/(Fs/2), 'bandpass');

correctedSig=filter(b,a, fliplr(sig));
correctedSig= fliplr(correctedSig);


function [correctedWFmat]=correctedMat (WFmat, Fs)

HP=300;
LP=6000;
[b, a] = butter(4, [HP LP]/(Fs/2), 'bandpass');

% %zero padding
% zeroNum=10;
% WFmat=[zeros(size(WFmat,1),zeroNum), WFmat, zeros(size(WFmat,1),zeroNum) ];

correctedWFmat=filter(b,a, fliplr(WFmat)');
correctedWFmat= fliplr(correctedWFmat');

% %zeros removal
% correctedWFmat (:,1:zeroNum)=[];
% correctedWFmat(:,(size(correctedWFmat,2)-zeroNum+1):size(correctedWFmat,2))=[];

a=1;