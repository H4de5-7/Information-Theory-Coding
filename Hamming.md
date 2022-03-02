# MATLAB实现7,4汉明码的编码解码纠错及BER的分析

本科时信息论与编码的作业

课程为

```
Information Theory & Coding
Vaibhav Kumar, PhD
School of Electrical & Electronic Engineering
University College Dublin – The Republic of Ireland
```

下述作业也是老师布置的，如果涉及到版权问题我会删掉该博客。

流程逻辑如下：

数据->汉明编码->BPSK调制->加入AWGN（模拟传输时的噪声）->BPSK解调->汉明解码->输出数据

给定生成矩阵为

```
G=[1 1 0 1 0 0 0;
   0 1 1 0 1 0 0;
   1 1 1 0 0 1 0;
   1 0 1 0 0 0 1]; %generator matrix
```

给定的syndrome和coset leader的对应关系为

syndrome     | coset leader
-------- | -----
000  | 0000000
001  | 0010000
010  | 0100000
011  | 0000100
100  | 1000000
101  | 0000001
110  | 0001000
111  | 0000010

</br>

## 代码

```matlab
clc,clear
Eb_N0_log=0:0.5:10; %dB
Eb_N0=10.^(Eb_N0_log/10);
R_coded=4/7; %R
BER_coded=zeros(1,21); %store BER in different Eb/N0
G=[1 1 0 1 0 0 0;
   0 1 1 0 1 0 0;
   1 1 1 0 0 1 0;
   1 0 1 0 0 0 1]; %generator matrix
H=[1 0 0 1 0 1 1;
   0 1 0 1 1 1 0;
   0 0 1 0 1 1 1]; %parity-check matrix
for i=1:21 %calculate BER in different Eb/N0
    SNR_coded=2*R_coded*Eb_N0(i); %calculate SNR
    message=round(rand(1,1000000)); %generate original data
    px_w=1; %signal's power
    pn_w=px_w/SNR_coded; %noise's power
    error_number=0; %error bit number
    for index=1:250000 %250000 blocks
        m=[message(4*index-3),message(4*index-2),message(4*index-1),message(4*index)]; %one block data
        v=mod(m*G,2); %generate v
        message_modul=1-v*2; %0->1,1->-1 modulation
        r=message_modul+sqrt(pn_w)*randn(1,7); %add white Gaussian noise
        r=r<0; %demodulation
        s=mod(r*H',2); %generate syndrome
        syndrome=char(s+'0'); %find correspond error pattern
        e=zeros(1,7);
        switch syndrome
            case '000'
                e=[0 0 0 0 0 0 0];
            case '001'
                e=[0 0 1 0 0 0 0];
            case '010'
                e=[0 1 0 0 0 0 0];
            case '011'
                e=[0 0 0 0 1 0 0];
            case '100'
                e=[1 0 0 0 0 0 0];
            case '101'
                e=[0 0 0 0 0 0 1];
            case '110'
                e=[0 0 0 1 0 0 0];
            case '111'
                e=[0 0 0 0 0 1 0];
            otherwise
                disp('wrong syndrome')
        end
        signal=mod(r+e,2); %correct recieved data
        error_number=error_number+sum(sum(signal(1,4:7)~=m)); %find the number of error bits
    end
    BER_coded(i)=error_number/1000000;
end
semilogy(Eb_N0_log,BER_coded);
xlabel('Eb/N0 [dB]');
ylabel('BER');
title('BER versus Eb/N0');

BER =0.5*erfc(sqrt(2*Eb_N0)/sqrt(2)); %Q function
hold on;
semilogy(Eb_N0_log,BER);
xlabel('Eb/N0 [dB]');
ylabel('BER');
title('BER versus Eb/N0');

R_uncoded=1; %R
BER_uncoded=zeros(1,18);
for i=1:21 %calculate BER in different Eb/N0
    SNR_uncoded=2*R_uncoded*Eb_N0(i); %SNR
    message=round(rand(1,1000000)); %generate original data
    message_modul=1-message.*2; %0->1,1->-1 modulation
    px_w=1; %signal's power
    pn_w=px_w/SNR_uncoded; %noise's power
    r=message_modul+sqrt(pn_w)*randn(1,1000000); %add white Gaussian noise
    error_number=0;
    r=r<0; %demodulation
    error_number=error_number+sum(sum(r~=message)); %find the number of error bits
    BER_uncoded(i)=error_number/1000000;
end
hold on;
semilogy(Eb_N0_log,BER_uncoded);
xlabel('Eb/N0 [dB]');
ylabel('BER ');
title('BER versus Eb/N0');
legend('Coded(simulation)','Uncoded(theory)','Uncoded(simulation)');
```


## 计算BER

1、随机生成1000000比特，分为250000 block，每个block称为m。

2、进行汉明编码。

3、进行BPSK调制。

4、添加高斯白噪声

![高斯白噪声](https://img-blog.csdnimg.cn/a6bb64961b9d4cb08f6e9da0f964e2c6.png?x-oss-process=image/watermark,type_d3F5LXplbmhlaQ,shadow_50,text_Q1NETiBASDRkZTU=,size_20,color_FFFFFF,t_70,g_se,x_16)

5、BPSK解码。

6、选择coset leader e并计算得出接收到的数据为r+e。计算错误比特的个数。

7、计算BER。

结果：

![结果1](https://img-blog.csdnimg.cn/92a0f448089d40ad8669a7dd068b650f.png?x-oss-process=image/watermark,type_d3F5LXplbmhlaQ,shadow_50,text_Q1NETiBASDRkZTU=,size_20,color_FFFFFF,t_70,g_se,x_16)
</br>

## 计算不使用汉明编码情况时的BER（模拟环境与理论情况）

1、随机生成1000000比特。

2、进行BPSK调制。

3、添加高斯白噪声。

4、BPSK解码。

5、计算错误比特的个数。

6、利用第五部分结果计算模拟环境下的BER。

6、使用erfc函数计算理论情况下的BER。

结果：

![在这里插入图片描述](https://img-blog.csdnimg.cn/4fdd5825843e496a976924007c5da0cd.png?x-oss-process=image/watermark,type_d3F5LXplbmhlaQ,shadow_50,text_Q1NETiBASDRkZTU=,size_20,color_FFFFFF,t_70,g_se,x_16)
</br>

## 总结

可以看到当Eb/N0大于5.8dB的时候汉明编码系统比未编码系统的BER低




