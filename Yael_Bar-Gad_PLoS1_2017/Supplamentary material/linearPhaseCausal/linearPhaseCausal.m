function linearPhaseCausal()

range = [-10:10];
k = 7;

figure;
subplot(2,2,1)
stem(range, LPC(range,k));
title('Linear phase - Causal');

subplot(2,2,2)
stem(range, LPNC(range,k));
title('Linear phase - Non causal');

subplot(2,2,3)
stem(range, NLPC(range,k));
title('Non linear phase - Causal');

subplot(2,2,4)
stem(range, NLPNC(range,k));
title('Non linear phase - Non causal');

IBG = 1;

function y = LPC(x, n)
    y=zeros(1,length(x));
    y (x>=0 & x<=n-1) = 1/n;

function y = LPNC(x, n)
    y=zeros(1,length(x));
    y(x>=(-n-1)/2 & x<=(n-1)/2) = 1/n;

function y = NLPC(x, n)
    y=zeros(1,length(x));
   for i=1:length(x)
       if  (x(i)>=0 && x(i)<=n-2)
           y(i) = 2^(-x(i)-1);
       elseif (x(i)==n-1)
            y(i) = 2^(-x(i));
       end
   end
   
function y = NLPNC(x, n)
    y=zeros(1,length(x));
    o = round(n/2);
   for i=1:length(x)
       if  (x(i)>=-o && x(i)<=n-2-o)
           y(i) = 2^(-x(i)-1+o);
       elseif (x(i)==n-1-o)
            y(i) = 2^(-x(i)+o);
       end
   end
    
    