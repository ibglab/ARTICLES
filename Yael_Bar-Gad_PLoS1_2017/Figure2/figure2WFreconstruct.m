function []= figure2WFreconstruct (meanWF)
%run figure2meanWF

Fs=44000;
mFs2 = 2*Fs/1000;

numberOfSines=9; %for summation into spike shape, including DC
sinesForPresentation=[1,3,5]+1;%3 sine functions

%WAVEFORM RECONSTRUCTION
%reconstruct spike shape by a summation of sinusoidal functions
dft=fft(meanWF(2:end));
window=length(dft); %ms
t=[0:1/mFs2:3-(3/mFs2)];
sum=zeros(1,3*mFs2-2)+real(dft(1)/(window));
for i=2:numberOfSines
    sum=(sum+(1/(mFs2/2))*abs(dft(i)).*sin([2*pi*(i-1)*t]+pi/2+unwrap(angle(dft(i)))));
end

sum=([zeros(1,mFs2), sum, zeros(1,mFs2)]); %reconstructed spike shape


%Sinusoidal functions for presentation
zPad = zeros(1,88);
f1=((1/(mFs2/2))*abs(dft(sinesForPresentation(1))).*sin([2*pi*(sinesForPresentation(1)-1)*t]+pi/2+unwrap(angle(dft((sinesForPresentation(1)))))));
f1=([zPad f1 zPad]);
f2=((1/(88/2))*abs(dft(sinesForPresentation(2))).*sin([2*pi*(sinesForPresentation(2)-1)*t]+pi/2+unwrap(angle(dft((sinesForPresentation(2)))))));
f2=([zPad f2 zPad]);
f3=((1/(88/2))*abs(dft(sinesForPresentation(3))).*sin([2*pi*(sinesForPresentation(3)-1)*t]+pi/2+unwrap(angle(dft((sinesForPresentation(3)))))));
f3=([zPad f3 zPad]);

%FILTERS
%NLP filter
HP =300; LP =  6000;
[b, a] = butter(4, [HP LP]/(Fs/2), 'bandpass');
f1NonLinear=filter(b, a, f1);
f2NonLinear=filter(b, a, f2);
f3NonLinear=filter(b, a, f3);
fsumNonLinear=filter(b, a, sum);

%ZP filter
[bbp, abp] = butter(4, [HP LP]/(Fs/2), 'bandpass');
f1Linear=filtfilt(bbp, abp, f1);
f2Linear=filtfilt(bbp, abp, f2);
f3Linear=filtfilt(bbp, abp, f3);
fsumLinear=filtfilt(bbp, abp, sum);

%FIR LP filter
bfir=fir1(20, [300 6000]/(Fs/2), 'bandpass');
f1FIR=filter(bfir,1,f1);
f2FIR=filter(bfir,1,f2);
f3FIR=filter(bfir,1,f3);
fsumFIR=filter(bfir,1,sum);

freq1=(sinesForPresentation(1)-1)*((Fs/window));
freq2=(sinesForPresentation(2)-1)*((Fs/window));
freq3=(sinesForPresentation(3)-1)*((Fs/window));


% PLOT FIGURE 
xmin=2*window+1; xmax=xmin+window-1;
ymin=-150; ymax=75;

figure; hold on;
subplot (4,3,1); hold on;
plot (f1,'k');
plot (f1Linear ,'r');
axis([xmin xmax ymin ymax]);
title ('ZP');
ylabel ([num2str(freq1) 'Hz']);
set(gca,'xtick',[])

subplot (4,3,2); hold on;
plot (f1,'k');
plot (f1FIR ,'c');
axis([xmin xmax ymin ymax]);
title ('LP');
set(gca,'xtick',[])

subplot (4,3,3); hold on;
plot (f1,'k');
plot (f1NonLinear,'b');
axis([xmin xmax ymin ymax]);
title ('NLP');
set(gca,'xtick',[])

subplot (4,3,4); hold on;
plot (f2,'k');
plot (f2Linear ,'r');
axis([xmin xmax ymin ymax]);
ylabel ([num2str(freq2) 'Hz']);
set(gca,'xtick',[])

subplot (4,3,5); hold on;
plot (f2,'k');
plot (f2FIR,'c');
axis([xmin xmax ymin ymax]);
set(gca,'xtick',[])

subplot (4,3,6); hold on;
plot (f2,'k');
plot (f2NonLinear ,'b');
axis([xmin xmax ymin ymax]);
set(gca,'xtick',[])

subplot (4,3,7); hold on;
plot (f3,'k');
plot (f3Linear ,'r');
axis([xmin xmax ymin ymax]);
ylabel ([num2str(freq3) 'Hz']);
set(gca,'xtick',[])

subplot (4,3,8); hold on;
plot (f3,'k');
plot (f3FIR ,'c');
axis([xmin xmax ymin ymax]);
set(gca,'xtick',[])

subplot (4,3,9); hold on;
plot (f3,'k');
plot (f3NonLinear ,'b');
axis([xmin xmax ymin ymax]);
set(gca,'xtick',[])

subplot (4,3,10); hold on;
plot (sum,'k');
plot ( fsumLinear,'r');
axis([xmin xmax ymin ymax]);
ylabel ('Sum');
set(gca,'xtick',[])
xlabel('T=2ms');


subplot (4,3,11); hold on;
plot (sum,'k');
plot (fsumFIR ,'c');
axis([xmin xmax ymin ymax]);
set(gca,'xtick',[])
xlabel('T=2ms');

subplot (4,3,12); hold on;
plot (sum,'k');
plot (fsumNonLinear ,'b');
axis([xmin xmax ymin ymax]);
set(gca,'xtick',[])
xlabel('T=2ms');


