file = 'C:\Users\13599\Downloads\202100172010.wav';
[y, Fs] = audioread(file);
% 播放音频
sound(y, Fs);
t = (0:length(y)-1)/Fs;	% 时间轴
figure;
subplot(3,1,1);
plot(t, y);
xlabel('时间/s');
ylabel('幅度');	%归一化幅度值,在-1到1之间
title('原始音频信号的时域波形');
% 绘制原始信号的频谱
N=length(y); 
Yk=fft(y,N);
Yk_abs=abs(Yk);
Yk_ang=angle(Yk); 
subplot(3,1,2); 
stem(Yk_abs);
xlabel('k');
ylabel('幅度');
title('原始音频信号FFT的幅度谱');
subplot(3,1,3);
stem(Yk_ang); 
xlabel('k');
ylabel('相位');
title('原始音频信号FFT的相位谱');
%对原始信号进行延时
% 读入音频信号和采样率
[x, fs] = audioread(file);
% 延时一定时间（3s）
delay_sec =3;
delay_samples = round(delay_sec * fs);
% 将整个信号插值成一个函数
t1 = (0:length(x)-1)'/fs;
f = griddedInterpolant(t1, x);
t1_new = (0:length(x)+delay_samples-1)'/fs;
% 将信号延时
x_new = f(t1_new-delay_sec); 
x_new(1:delay_samples) = 0;	% 将前面的填充为0
% 播放和保存
sound(x_new, fs); 
audiowrite('video_delay.wav', x_new, fs);
% 绘制延时信号时域波形
t2 = (0:length(x_new)-1)/Fs;	% 时间轴
figure; 
subplot(3,1,1);
plot(t2, x_new);
xlabel('时间/s');
ylabel('幅度');	%归一化幅度值,在-1到1之间title('延时音频信号的时域波形');
% 绘制延时信号的频谱
N1=length(x_new); 
Xk=fft(x_new,N1); 
Xk_abs=abs(Xk);
Xk_ang=angle(Xk); 
subplot(3,1,2);
stem(Xk_abs); 
xlabel('k');
ylabel('幅度');
title('音频信号FFT的幅度谱');
subplot(3,1,3); 
stem(Xk_ang);
xlabel('k');
ylabel('相位');
title('音频信号FFT的相位谱');
% 获取回声长度
echoLength = round(0.2 * Fs); % 假设回声长度为 0.2s
% 设置延迟时间和增益
delay = 0.1
gain = 0.8;
% 构造一个单回声滤波器
h = zeros(echoLength,1); 
h(1) = 1; % 单位脉冲响应
h(round(delay*Fs)+1) = gain; % 添加回声
% 对输入信号进行滤波
y_echo = filter(h,1,y);
sound(y_echo, Fs);
t_echo = (0:length(y_echo)-1)/Fs;	% 时间轴
figure; subplot(3,1,1);
plot(t_echo, y_echo);
xlabel('时间/s');
ylabel('幅度');%归一化幅度值,在-1到1之间
title('单回声信号的时域波形');
% 绘制单回声信号的频谱
N=length(y_echo);
Yk_echo=fft(y_echo,N); 
Yk_echo_abs=abs(Yk_echo);
Yk_echo_ang=angle(Yk_echo);
subplot(3,1,2); 
stem(Yk_echo_abs);
xlabel('k');
ylabel('幅度');
title('单回声信号FFT的幅度谱');
subplot(3,1,3);
stem(Yk_echo_ang);
xlabel('k');
ylabel('相位');
title('单回声信号FFT的相位谱');
%多重回声滤波器
% 定义滤波器的参数
delays = [0.05 0.1 0.15 0.2 0.25]; % 回声延迟时间
gains = [0.8 0.6 0.3 0.1 0.03]; % 回声增益因子

% 循环创建多个单回声滤波器
echoFilters = cell(1, length(delays)); 
for n = 1:length(delays)
    disp(4)
% 第 n 个单回声滤波器的延迟和增益
d = delays(n); 
g = gains(n);
% 创建一个一阶 FIR 滤波器
H = [1; zeros(round(d*Fs), 1); 
g*ones(length(y)-round(d*Fs)-1, 1)];
echoFilters{n} = dsp.FIRFilter('Numerator', H');
end
% 对输入信号应用多重回声滤波器
y_duo = y;
for n = 1:length(echoFilters) 
    disp(2);
    disp(length(echoFilters));
    y_duo = echoFilters{n}(y_duo);
end
disp(7)
% 播放加有多重回声的音频信号
sound(y_duo, Fs);
% 绘制多重回声信号的频谱
N_duo=length(y_duo);
Yk_duo=fft(y_duo,N_duo);
Yk_duo_abs=abs(Yk_duo);
Yk_duo_ang=angle(Yk_duo); 
figure;
subplot(2,1,1); 
stem(Yk_duo_abs);
xlabel('k');
ylabel('幅度');
title('多重回声信号FFT的幅度谱'); 
subplot(2,1,2); 
stem(Yk_duo_ang);
xlabel('k');
ylabel('相位');
title('多重回声信号FFT的相位谱');
disp(3)
%无限回声滤波器
% 定义滤波器的参数
delay1 = 0.2;	% 回声延迟时间
gain1 = 0.8;	% 回声增益因子

% 创建一个一阶 FIR 滤波器
fil = [1; zeros(round(delay1*Fs),1);
gain1;
zeros(length(y),1)];
hin = dsp.FIRFilter('Numerator',fil');

% 对输入信号应用无限回声滤波器
y_in = hin(y);

% 播放加有无限回声的音频信号
sound(y_in, Fs);

% 绘制单回声信号时域波形
t_in = (0:length(y_in)-1)/Fs;	% 时间轴
figure; subplot(3,1,1);
plot(t_in, y_in);
xlabel('时间/s');
ylabel('幅度');	%归一化幅度值,在-1到1之间title('无限回声信号的时域波形');

% 绘制无限回声信号的频谱N_in=length(y_in); Yk_in=fft(y_in,N_in);
Yk_in_abs=abs(Yk_in); 
Yk_in_ang=angle(Yk_in);
subplot(3,1,2); 
stem(Yk_in_abs);
xlabel('k');
ylabel('幅度');
title('无限回声信号FFT的幅度谱');
subplot(3,1,3); 
stem(Yk_in_ang);
xlabel('k');
ylabel('相位');
title('无限回声信号FFT的相位谱');
% 指定低通滤波器的通带截止频率、阻带截止频率和衰减特性
fp = 1000;  % 通带截止频率fc = 1200;  % 阻带截止频率As = 100;   % 阻带最小衰减Ap = 1; % 通 带 最 大 衰 减
% 创建滤波器规格对象
lowpassSpecs = fdesign.lowpass('Fp,Fst,Ap,Ast',fp,fc,Ap,As,5000);
% 根据滤波器规格对象设计滤波器
lowpassFilter = design(lowpassSpecs);
% 频率响应分析
freqz(lowpassFilter);
% 滤波处理
y_low = filter(lowpassFilter, y);
% 播放滤波后的音频信号
%sound(y_low,Fs);

% 绘制低通信号时域波形
t_low = (0:length(y_low)-1)/Fs;	% 时间轴
figure; subplot(3,1,1); plot(t_low, y_low);
xlabel('时间/s');
ylabel('幅度');%归一化幅度值,在-1到1之间
title('低通信号的时域波形');

% 绘制无限回声信号的频谱N_low=length(y_low); Yk_low=fft(y_low,N_low);
Yk_low_abs=abs(Yk_low); 
Yk_low_ang=angle(Yk_low); 
subplot(3,1,2); 
stem(Yk_low_abs);
xlabel('k');
ylabel('幅度');
title('低通信号FFT的幅度谱'); subplot(3,1,3);
stem(Yk_low_ang); xlabel('k');
ylabel('相位');
title('低通信号FFT的相位谱');

% 指定高通滤波器的通带截止频率、阻带截止频率和衰减特性
fc1 = 4000; % 通带截止频率fp1 = 3500; % 阻带截止频率As1 = 100; % 阻带最小衰减Ap1 = 1; % 通带最大衰减

% 计算滤波器的阶数和截止频率
[n, Wn] = buttord(2*fp1/10000, 2*fc1/10000, Ap1, As1);
disp(1)

% 根据阶数和截止频率设计滤波器
[b, a] = butter(n, Wn, 'high');
% 绘制高通滤波器的幅度响应和相位响应
freqz(b, a);
% 绘制振幅响应曲线和相位响应曲线
[Hh, f] = freqz(b, a, 512, 10000);
amplitudeResponse = 20 * log10(abs(Hh));
phaseResponse = unwrap(angle(Hh)) .* 180 ./ pi;

figure; subplot(2,1,1);
plot(f, amplitudeResponse); 
xlabel('Frequency (Hz)'); 
ylabel('Amplitude (dB)'); 
title('Amplitude Response');

subplot(2,1,2); 
plot(f, phaseResponse);
xlabel('Frequency (Hz)'); 
ylabel('Phase (degrees)'); 
title('Phase Response');

% 滤波处理
y_high= filter(b,a,y);
% 播放滤波后的音频信号
sound(y_high,Fs);
% 绘制高通信号时域波形
t_high = (0:length(y_high)-1)/Fs;	% 时间轴
figure; subplot(3,1,1); 
plot(t_high, y_high);
xlabel('时间/s');
ylabel('幅度');	%归一化幅度值,在-1到1之间title('高通信号的时域波形');

% 绘制高通信号的频谱N_high=length(y_high); Yk_high=fft(y_high,N_high);
Yk_high_abs=abs(Yk_high); 
Yk_high_ang=angle(Yk_high); 
subplot(3,1,2);
stem(Yk_high_abs); 
xlabel('k');
ylabel('幅度');
title('高通信号FFT的幅度谱'); 
subplot(3,1,3);
stem(Yk_high_ang); 
xlabel('k');
ylabel('相位');
title('高通信号FFT的相位谱');

%设计带通滤波器
%fp1=1200Hz，fp2=3000hZ，fc1=1000Hz，fc2=3200Hz，As=100dB，Ap=1dB,取采样频率fs=8000Hz
wp=[0.3,0.75];
ws=[0.25,0.8];
ap=1; 
as=100;
[n1,wc1]=buttord(wp,ws,ap,as); 
[b1,a1]=butter(n1,wc1);
figure;
[h1,w1]=freqz(b1,a1);
db1=20*log10(abs(h1)/max(abs(h1)));
subplot(3,1,1);
plot(w1/pi,abs(h1)); 
xlabel('pi');
ylabel('幅度谱');
axis([0,1,0,1.1]);
subplot(3,1,2);
plot(w1/pi,angle(h1));
xlabel('pi');
ylabel(' 相 位 谱 ');
subplot(3,1,3);
plot(w1/pi,db1);
xlabel('pi');
ylabel('dB');
% 滤波处理
y_dai= filter(b1,a1,y);
% 播放滤波后的音频信号
sound(y_dai,Fs);

% 绘制带通信号时域波形
t_dai = (0:length(y_dai)-1)/Fs;	% 时间轴
figure; subplot(3,1,1);
plot(t_dai, y_dai);
xlabel('时间/s');
ylabel('幅度');%归一化幅度值,在-1到1之间
title('带通信号的时域波形');

% 绘制高通信号的频谱
N_dai=length(y_dai); 
Yk_dai=fft(y_dai,N_dai);
Yk_dai_abs=abs(Yk_dai); 
Yk_dai_ang=angle(Yk_dai); 
subplot(3,1,2); 
stem(Yk_dai_abs); 
xlabel('k');
ylabel('幅度');
title('带通信号FFT的幅度谱');
subplot(3,1,3);
stem(Yk_dai_ang);
xlabel('k');
ylabel('相位');
title('带通信号FFT的相位谱');