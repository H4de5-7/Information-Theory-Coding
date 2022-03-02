clc,clear
Eb_N0_log=0:0.5:10; %dB
Eb_N0=10.^(Eb_N0_log/10);
R_coded=11/16; %R
BER_coded=zeros(1,21); %store BER in different Eb/N0
G=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; %v0
   0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1; %v4
   0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1; %v3
   0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1; %v2
   0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1; %v1
   0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1; %v3v4
   0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 1; %v2v4
   0 0 0 0 0 0 0 0 0 1 0 1 0 1 0 1; %v1v4
   0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1; %v2v3
   0 0 0 0 0 1 0 1 0 0 0 0 0 1 0 1; %v1v3
   0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1; %v1v2
   ]; %generator matrix
for i=1:21 %calculate BER in different Eb/N0
    SNR_coded=2*R_coded*Eb_N0(i); %calculate SNR
    message=round(rand(1,1980000)); %generate original data u = (a0 a4 a3 a2 a1 a34 a24 a14 a23 a13 a12),
    px_w=1; %signal's power
    pn_w=px_w/SNR_coded; %noise's power
    error_number=0; %error bit number
    for index=1:180000 %90000 blocks
        m=message(1,11*index-10:11*index);
        %m=[message(11*index-10),message(11*index-9),message(11*index-8),message(11*index-7),message(11*index-6),message(11*index-5),message(11*index-4),message(11*index-3),message(11*index-2),message(11*index-1),message(11*index)]; %one blobk data
        v=mod(m*G,2); %generate v
        message_modul=1-v*2; %0->1,1->-1 modulation
        r=message_modul+sqrt(pn_w)*randn(1,16); %add white Gaussian noise
        r=r<0; %demodulation
        a12=sum([mod(r(1)+r(2)+r(3)+r(4),2),mod(r(5)+r(6)+r(7)+r(8),2),mod(r(9)+r(10)+r(11)+r(12),2),mod(r(13)+r(14)+r(15)+r(16),2)]==1);
        if a12>2
            a12=1;
        elseif a12<2
            a12=0;
        else
            a12=round(rand());
        end %check-sums for a12
        a13=sum([mod(r(1)+r(2)+r(5)+r(6),2),mod(r(3)+r(4)+r(7)+r(8),2),mod(r(9)+r(10)+r(13)+r(14),2),mod(r(11)+r(12)+r(15)+r(16),2)]==1);
        if a13>2
            a13=1;
        elseif a13<2
            a13=0;
        else
            a13=round(rand());
        end %check-sums for a13
        a23=sum([mod(r(1)+r(3)+r(5)+r(7),2),mod(r(2)+r(4)+r(6)+r(8),2),mod(r(9)+r(11)+r(13)+r(15),2),mod(r(10)+r(12)+r(14)+r(16),2)]==1);
        if a23>2
            a23=1;
        elseif a23<2
            a23=0;
        else
            a23=round(rand());
        end %check-sums for a23
        a14=sum([mod(r(1)+r(2)+r(9)+r(10),2),mod(r(3)+r(4)+r(11)+r(12),2),mod(r(5)+r(6)+r(13)+r(14),2),mod(r(7)+r(8)+r(15)+r(16),2)]==1);
        if a14>2
            a14=1;
        elseif a14<2
            a14=0;
        else
            a14=round(rand());
        end %check-sums for a14
        a24=sum([mod(r(1)+r(3)+r(9)+r(11),2),mod(r(2)+r(4)+r(10)+r(12),2),mod(r(5)+r(7)+r(13)+r(15),2),mod(r(6)+r(8)+r(14)+r(16),2)]==1);
        if a24>2
            a24=1;
        elseif a24<2
            a24=0;
        else
            a24=round(rand());
        end %check-sums for a24
        a34=sum([mod(r(1)+r(5)+r(9)+r(13),2),mod(r(2)+r(6)+r(10)+r(14),2),mod(r(3)+r(7)+r(11)+r(15),2),mod(r(4)+r(8)+r(12)+r(16),2)]==1);
        if a34>2
            a34=1;
        elseif a34<2
            a34=0;
        else
            a34=round(rand());
        end %check-sums for a34
        
        r=mod(r-a12*G(11,:)-a13*G(10,:)-a23*G(9,:)-a14*G(8,:)-a24*G(7,:)-a34*G(6,:),2);
        a1=sum([mod(r(1)+r(2),2),mod(r(3)+r(4),2),mod(r(5)+r(6),2),mod(r(7)+r(8),2),mod(r(9)+r(10),2),mod(r(11)+r(12),2),mod(r(13)+r(14),2),mod(r(15)+r(16),2),]==1);
        if a1>4
            a1=1;
        elseif a1<4
            a1=0;
        else
            a1=round(rand());
        end %check-sums for a1
        a2=sum([mod(r(1)+r(3),2),mod(r(2)+r(4),2),mod(r(5)+r(7),2),mod(r(6)+r(8),2),mod(r(9)+r(11),2),mod(r(10)+r(12),2),mod(r(13)+r(15),2),mod(r(14)+r(16),2),]==1);
        if a2>4
            a2=1;
        elseif a2<4
            a2=0;
        else
            a2=round(rand());
        end %check-sums for a2
        a3=sum([mod(r(1)+r(5),2),mod(r(2)+r(6),2),mod(r(3)+r(7),2),mod(r(4)+r(8),2),mod(r(9)+r(13),2),mod(r(10)+r(14),2),mod(r(11)+r(15),2),mod(r(12)+r(16),2),]==1);
        if a3>4
            a3=1;
        elseif a3<4
            a3=0;
        else
            a3=round(rand());
        end %check-sums for a3
        a4=sum([mod(r(1)+r(9),2),mod(r(2)+r(10),2),mod(r(3)+r(11),2),mod(r(4)+r(12),2),mod(r(5)+r(13),2),mod(r(6)+r(14),2),mod(r(7)+r(15),2),mod(r(8)+r(16),2),]==1);
        if a4>4
            a4=1;
        elseif a4<4
            a4=0;
        else
            a4=round(rand());
        end %check-sums for a4
        
        r=mod(r-a1*G(5,:)-a2*G(4,:)-a3*G(3,:)-a4*G(2,:),2);
        a0=sum(r==1);
        if a0>8
            a0=1;
        elseif a0<8
            a0=0;
        else
            a0=round(rand());
        end %check-sums for a0
        signal=[a0,a4,a3,a2,a1,a34,a24,a14,a23,a13,a12];
        error_number=error_number+sum(sum(signal~=m)); %find the number of error bits
    end
    BER_coded(i)=error_number/1980000;
end
semilogy(Eb_N0_log,BER_coded);
xlabel('Eb/N0 [dB]');
ylabel('BER');
title('BER versus Eb/N0');

Eb_N0_log=0:0.5:10; %dB
Eb_N0=10.^(Eb_N0_log/10);
BER =0.5*erfc(sqrt(2*Eb_N0)/sqrt(2)); %Q function
hold on;
semilogy(Eb_N0_log,BER);
xlabel('Eb/N0 [dB]');
ylabel('BER');
title('BER versus Eb/N0');

R_uncoded=1; %R
BER_uncoded=zeros(1,21);
for i=1:21 %calculate BER in different Eb/N0
    SNR_uncoded=2*R_uncoded*Eb_N0(i); %SNR
    message=round(rand(1,1980000)); %generate original data
    message_modul=1-message.*2; %0->1,1->-1 modulation
    px_w=1; %signal's power
    pn_w=px_w/SNR_uncoded; %noise's power
    r=message_modul+sqrt(pn_w)*randn(1,1980000); %add white Gaussian noise
    error_number=0;
    r=r<0; %demodulation
    error_number=error_number+sum(sum(r~=message)); %find the number of error bits
    BER_uncoded(i)=error_number/1980000;
end
hold on;
semilogy(Eb_N0_log,BER_uncoded);
xlabel('Eb/N0 [dB]');
ylabel('BER ');
title('BER versus Eb/N0');
legend('RM-coded(simulation)','Uncoded(theory)','Uncoded(simulation)');

%I extend the horizontal coordinate of Uncoded(theory) and Uncoded(simulation)

% Eb_N0_log=0:0.5:11; %dB
% Eb_N0=10.^(Eb_N0_log/10);
% BER =0.5*erfc(sqrt(2*Eb_N0)/sqrt(2)); %Q function
% hold on;
% semilogy(Eb_N0_log,BER);
% xlabel('Eb/N0 [dB]');
% ylabel('BER');
% title('BER versus Eb/N0');
% 
% R_uncoded=1; %R
% BER_uncoded=zeros(1,23);
% for i=1:23 %calculate BER in different Eb/N0
%     SNR_uncoded=2*R_uncoded*Eb_N0(i); %SNR
%     message=round(rand(1,1980000)); %generate original data
%     message_modul=1-message.*2; %0->1,1->-1 modulation
%     px_w=1; %signal's power
%     pn_w=px_w/SNR_uncoded; %noise's power
%     r=message_modul+sqrt(pn_w)*randn(1,1980000); %add white Gaussian noise
%     error_number=0;
%     r=r<0; %demodulation
%     error_number=error_number+sum(sum(r~=message)); %find the number of error bits
%     BER_uncoded(i)=error_number/1980000;
% end
% hold on;
% semilogy(Eb_N0_log,BER_uncoded);
% xlabel('Eb/N0 [dB]');
% ylabel('BER ');
% title('BER versus Eb/N0');
% legend('RM-coded(simulation)','Uncoded(theory)','Uncoded(simulation)');
