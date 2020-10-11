clear;
%% create the carrier
n = 5; resolution_n = 11;
fs = 50000000;
fc = fs/(2^n);
x=linspace(0,2*pi,2^n);		   %32 point, carrier frequency is 1Mhz
y1=cos(x)+1;					   
carrier_cos=ceil(y1*(2^resolution_n-1));				

y2=sin(x)+1;					   
carrier_sin=ceil(y2*(2^resolution_n-1));				

% there is a straight signal in the carrier
%% create the source code
N = 10;    %% the number of QAM symbols
Code_N = N*10;
info=random_binary(Code_N);
info = reshape(info, [N, 10]);
I_info = info(:,1:5);
Q_info = info(:,6:10);

%% transmit:Differential coding  for I_info
I_delay_info = I_info(:,1:4);
I_diff_info = I_info;
I_diff_info(:, 2:5) =  xor(I_delay_info, I_info(:, 2:5));

%% transmit:Gray code 2 binary
Transmit_I_Gray_bits = I_diff_info;
Transmit_I_two_bits(:,1) = Transmit_I_Gray_bits(:,1);
Transmit_Q_Gray_bits = Q_info;
Transmit_Q_two_bits(:,1) = Transmit_Q_Gray_bits(:,1);
for i = 2:5
    Transmit_I_two_bits(:,i) = xor(Transmit_I_Gray_bits(:,i), Transmit_I_two_bits(:, i-1));
    Transmit_Q_two_bits(:,i) = xor(Transmit_Q_Gray_bits(:,i), Transmit_Q_two_bits(:, i-1));
end


%% transmit:transform to symbol
I_symbol = zeros(N,1);
Q_symbol = zeros(N,1);
for i = 1:5
    I_symbol(:,1) = I_symbol(:,1) + 2^(5-i) * Transmit_I_two_bits(:,i);
    Q_symbol(:,1) = Q_symbol(:,1) + 2^(5-i) * Transmit_Q_two_bits(:,i);
end 

%%  transmit:multiple the carrier and add the I and Q
% code rate: 100kHz
f_k = 5;    % the frequency k between fc and fb
fb = fc/f_k;
I_symbol_repeat = repmat(I_symbol,1,f_k*2^n);
I_symbol_all = reshape(I_symbol_repeat', 1, f_k*2^n*N);
Q_symbol_repeat = repmat(Q_symbol,1,f_k*2^n);
Q_symbol_all = reshape(Q_symbol_repeat', 1, f_k*2^n*N);


%% transmit:transmit filter
[b,a]=butter(2,4*fb/fs);  
I_symbol_bshape=filtfilt(b,a,I_symbol_all);	
Q_symbol_bshape=filtfilt(b,a,Q_symbol_all);

%% transmit:finish the modulation
carrier_sin = repmat(carrier_sin,1,f_k*N);
carrier_cos = repmat(carrier_cos,1,f_k*N);

I_fpga_signal = (I_symbol_bshape-15.5).*(carrier_cos-(2^(resolution_n)-1));
Q_fpga_signal = (Q_symbol_bshape-15.5).*(carrier_sin-(2^(resolution_n)-1));

fpga_signal = I_fpga_signal + Q_fpga_signal;

subplot(3,1,1)
plot(I_fpga_signal)
legend('I signal')
subplot(3,1,2)
plot(Q_fpga_signal)
legend('Q signal')
subplot(3,1,3)
plot(fpga_signal)
legend('QAM signal')

%% channel:now is empty

snr = 40;
noise = max(abs(fpga_signal)) * wgn(1,length(fpga_signal),-snr);
Receive_signals = fpga_signal + noise; 

%% receive:multiple the carrier
Receive_I_signals = Receive_signals.*(carrier_cos-(2^(resolution_n)-1));
Receive_Q_signals = Receive_signals.*(carrier_sin-(2^(resolution_n)-1));

%% receive:LPF
[b,a]=butter(2,4*fb/fs);  
Receive_I_symbols=filtfilt(b,a,Receive_I_signals)*0.95;	
Receive_Q_symbols=filtfilt(b,a,Receive_Q_signals);	

%% receive:judge
judge_n=fix((.5:1:N)*f_k*(2^n));

Receive_I_codes1 = Receive_I_symbols/(0.5*((2^resolution_n-1)^2))+15.5;
Receive_Q_codes1 = Receive_Q_symbols/(0.5*((2^resolution_n-1)^2))+15.5;

Receive_I_codes2 = Receive_I_codes1(judge_n);
Receive_Q_codes2 = Receive_Q_codes1(judge_n);

hold on;

%% draw the constellation symbols
% scatter(Receive_I_codes2-15.5,Receive_Q_codes2-15.5,'.')
% hold off;  box on
% hold on;
% DrawCircle(0,0,15.5*sqrt(2),100,'r')
% title('1024QAM:constellation symbols');
%% receive:32--->2

Receive_I_bits = thirty_two2two(Receive_I_codes2);
Receive_Q_bits = thirty_two2two(Receive_Q_codes2);

Receive_I_two_bits = reshape(Receive_I_bits,[5,numel(Receive_I_bits)/5])';
Receive_Q_two_bits = reshape(Receive_Q_bits,[5,numel(Receive_Q_bits)/5])';
%% two2gray

Receive_I_gray_bits = Receive_I_two_bits;
Receive_I_gray_bits(:,1) = Receive_I_two_bits(:,1);
Receive_Q_gray_bits = Receive_Q_two_bits;
Receive_Q_gray_bits(:,1) = Receive_Q_two_bits(:,1);
for i = 2:5
    Receive_I_gray_bits(:,i) = xor(Receive_I_two_bits(:,i-1), Receive_I_two_bits(:, i));
    Receive_Q_gray_bits(:,i) = xor(Receive_Q_two_bits(:,i-1), Receive_Q_two_bits(:, i));
end

%% I--decoding diffiential
Receive_I_dediff_bits  = Receive_I_gray_bits;
Receive_I_dediff_bits(:,1) = Receive_I_gray_bits(:,1);
for i=2:5
    Receive_I_dediff_bits(:,i) = xor(Receive_I_gray_bits(:,i),Receive_I_dediff_bits(:,i-1));
end

%% compose the all receive info
Receive_info = zeros(N,10);
Receive_info(:,1:5) = Receive_I_dediff_bits;
Receive_info(:,6:10) = Transmit_Q_Gray_bits;

%% calculate the BER
disp(sum(sum(abs(Receive_info-info)))/Code_N)

% subplot(3,1,1)
% plot(Receive_signals/(0.5*((2^resolution_n-1)^1))+15.5)
% hold on;plot(Receive_I_symbols/(0.5*((2^resolution_n-1)^2))+15.5);
% hold on;plot(Receive_Q_symbols/(0.5*((2^resolution_n-1)^2))+15.5);
% legend('QAM receice signal','receive I symbols','receive Q symbols');
% subplot(3,1,2)
% plot(Receive_I_signals/(0.5*((2^resolution_n-1)^2))+15.5);hold on;
% plot(Receive_I_symbols/(0.5*((2^resolution_n-1)^2))+15.5)
% legend('I receice signal','receive I symbols');
% subplot(3,1,3)
% plot(Receive_Q_signals/(0.5*((2^resolution_n-1)^2))+15.5);hold on;
% plot(Receive_Q_symbols/(0.5*((2^resolution_n-1)^2))+15.5)
% legend('Q receice signal','receive Q symbols');

y = Receive_signals/(0.5*((2^resolution_n-1)^1))+15.5;
n=length(y);     y=fft(y)/n;    y=abs(y(1:fix(n/2)))*2;
q=find(y<1e-04); y(q)=1e-04;    y=20*log10(y);
m=fs/fb;
f1=m/n;          f=0:f1:(length(y)-1)*f1;
subplot(3,1,1)
plot(f,y,'r');	
grid on;  
title('receive signal'); 	xlabel('f/fb'); 

y = Receive_I_signals/(0.5*((2^resolution_n-1)^2))+15.5;
n=length(y);     y=fft(y)/n;    y=abs(y(1:fix(n/2)))*2;
q=find(y<1e-04); y(q)=1e-04;    y=20*log10(y);
m=fs/fb;
f1=m/n;          f=0:f1:(length(y)-1)*f1;
subplot(3,1,2)
plot(f,y,'r');	
grid on;  
title('after multiple: receive I signal'); 	xlabel('f/fb'); 

y = Receive_I_symbols/(0.5*((2^resolution_n-1)^2))+15.5;
n=length(y);     y=fft(y)/n;    y=abs(y(1:fix(n/2)))*2;
q=find(y<1e-04); y(q)=1e-04;    y=20*log10(y);
m=fs/fb;
f1=m/n;          f=0:f1:(length(y)-1)*f1;
subplot(3,1,3)
plot(f,y,'r');	
grid on;  
title('after LPF: receive I symbols'); 	xlabel('f/fb'); 
