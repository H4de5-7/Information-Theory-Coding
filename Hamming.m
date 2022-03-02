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
