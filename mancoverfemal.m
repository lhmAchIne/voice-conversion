clear all, close all, clc;

framelength = 80;   %帧长
windowlength = 240; %窗长
P = 10;             %LPC预测阶数

[signal, fs] = wavread('man.wav');      %读取一段语音数据
L = length(signal);                     %语音数据长度
FrameNumber = floor(L/framelength)-2;   %帧数

%初始化
exc = zeros(L, 1);          %激励信号（预测误差）
zi_pre = zeros(P, 1);       %预测滤波器的状态
s_rec = zeros(L, 1);        %重建语音信号
exc_syn_t = zeros(L, 1);    %合成激励信号
s_syn_t = zeros(L, 1);      %合成语音
last_syn_t = 0;             %存储上一个（或多个）段的最后一个脉冲的下场
zi_syn_t = zeros(P, 1);     %合成滤波器状态
hw = hamming(windowlength); %汉明窗

%依次处理每帧信号
for n = 3:FrameNumber
    %计算预测系数
    s_w = signal(n*framelength - windowlength + 1:n*framelength).*hw;   %用汉明窗加权信号
    [A E] = lpc(s_w, P);                %计算P个预测系数
                                        %A是预测系数，E是用来计算合成激励能量
    s_f = signal((n-1)*framelength + 1:n*framelength);                  %取当前帧数据
    [exc1, zi_pre] = filter(A, 1, s_f, zi_pre);
    exc((n-1)*framelength + 1:n*framelength) = exc1;                    %获取激励信号
    s_Pitch = exc(n*framelength - 222:n*framelength);
    PT = findpitch(s_Pitch);            %计算基音周期
    G = sqrt(E*PT);                     %计算增益
    
    %将基因周期减少一半，将共振峰频率增加150Hz，重新合成语音
    PT1 = floor(PT/2);                  %减少基因周期
    poles = roots(A);
    deltaOMG = 150*2*pi/8000;
    
    %增加共振峰频率，实轴上方的极点逆时针转，下方顺时针转
    for p = 1:10
        if imag(poles(p))>0
            poles(p) = poles(p)*exp(j*deltaOMG);
        elseif imag(poles(p))<0
            poles(p) = poles(p)*exp(-j*deltaOMG);
        end
    end
    A1 = poly(poles);
    tempn_syn_t = [1:n*framelength - last_syn_t]';
    exc_syn1_t = zeros(length(tempn_syn_t), 1);
    exc_syn1_t(mod(tempn_syn_t, PT1) == 0) = G;         %计算该段脉冲
    exc_syn1_t = exc_syn1_t((n-1)*framelength - last_syn_t + 1:n*framelength - last_syn_t);
    [s_syn1_t, zi_syn_t] = filter(1, A1, exc_syn1_t, zi_syn_t);
    exc_syn_t((n-1)*framelength + 1:n*framelength) = exc_syn1_t;        %得到合成激励
    s_syn_t((n-1)*framelength + 1:n*framelength) = s_syn1_t;             %得到合成语音
    last_syn_t = last_syn_t + PT1*floor((n*framelength - last_syn_t)/PT1);
end

sound(signal, fs);
pause(2);
sound(s_syn_t, fs);
