% % AFFT Capstone Project - Shashank S Iyengar (M12934513)
% % MS - Mechanical Engineering
% % University of Cincinnati

close all
clear
clc
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 15)
set(0,'defaultlinelinewidth',1)
set(0,'DefaultLineMarkerSize', 6)
set(0,'defaultAxesFontWeight','bold')

load('TCAP.mat')
tt = xdata;

%========================Case 1========================
%======================================================
Ntime = 1024;
Nc = 0;
Navg = length(tt)/Ntime;
ff = ydata(1:3,:)/0.5;
yy = ydata(4:8,:)/0.1;
h = hann(Ntime)';

% 50% Overlap + Windowing + FFT
for yi=1:5
    for fi=1:3
        for i=1:Navg*2-1
            ff_fft{fi}{i} = fft(h.*ff(fi,(1+(i-1)*Ntime/2):((i+1)*Ntime/2)));
            yy_fft{yi}{i} = fft(h.*yy(yi,(1+(i-1)*Ntime/2):((i+1)*Ntime/2)));
        end
    end
end
% Power Spectrum Averaging
for yi=1:5
    for fi=1:3
        for i=1:Navg*2-1
            Gxf1(yi,fi,i,:) = yy_fft{1,yi}{1,i}.*conj(ff_fft{1,fi}{1,i});
        end
    end
end
for fi=1:3
    for yi=1:5
        for i=1:Navg*2-1
            Gfx1(fi,yi,i,:) = ff_fft{1,fi}{1,i}.*conj(yy_fft{1,yi}{1,i});
        end
    end
end
for yi=1:5
    for yii=1:5
        for i=1:Navg*2-1
            Gxx1(yi,yii,i,:) = yy_fft{1,yi}{1,i}.*conj(yy_fft{1,yii}{1,i});
        end
    end
end
for fi=1:3
    for fii=1:3
        for i=1:Navg*2-1
            Gff1(fi,fii,i,:) = ff_fft{1,fi}{1,i}.*conj(ff_fft{1,fii}{1,i});
        end
    end
end

Gxf = mean(Gxf1,3);
Gff = mean(Gff1,3);
Gfx = mean(Gfx1,3);
Gxx = mean(Gxx1,3);

% H1, H2 Algorithm
for i = 1:Ntime
    H1(:,:,1,i) = Gxf(:,:,1,i)/Gff(:,:,1,i);
    H2(:,:,1,i) = Gxx(:,:,1,i)/Gfx(:,:,1,i);
end
% Hv Algorithm
for yi=1:5
    for fi=1:3
        for cnt_freq = 1:Ntime/2
            [V D] = eig([Gff(:,:,1,cnt_freq).' Gxf(yi,:,1,cnt_freq).';Gfx(:,yi,1,cnt_freq).' Gxx(yi,yi,1,cnt_freq)]);
            eval = diag(D);
            seig = eval(1);
            vec = V(:,1);
            nv=1;
            for k=1:length(eval)
                if(eval(k)<seig)
                    seig=eval(k);
                    vec=V(:,k);
                    nv=k;
                end
            end
            vec = vec*(-1/vec(end,nv));
            Hv1(yi,:,1,cnt_freq) = vec(1:3,nv);
            Hv = permute(Hv1,[2 1 3 4]);
        end
    end
end

% FRFs 
n=1;
for yi=1:5
    for fi=1:3
            figure(n)
            subplot(2,1,1)
            semilogy(0:511,abs(squeeze(H1(yi,fi,1,1:Ntime/2))),'-r')
            hold on
            semilogy(0:511,abs(squeeze(H2(yi,fi,1,1:Ntime/2))),'-b')
            hold on
            semilogy(0:511,abs(squeeze(Hv(fi,yi,1,1:Ntime/2))),'--y')
            grid on
            title(sprintf('Case 1: FRF for H_q_%d_,_p_%d',yi,fi))
            xlabel(['$ Frequency\;\mathrm{[Hz]} $'],'interpreter','latex')
            ylabel(['$ Amplitude\;\mathrm{[g/lbf]} $'],'interpreter','latex')
            legend('H1','H2','Hv')
            subplot(2,1,2)
            plot(0:511,angle(squeeze(H1(yi,fi,1,1:Ntime/2)))*180/pi,'-r')
            hold on
            plot(0:511,angle(squeeze(H2(yi,fi,1,1:Ntime/2)))*180/pi,'-b')
            hold on
            plot(0:511,angle(squeeze(Hv(fi,yi,1,1:Ntime/2)))*180/pi,'--y')
            grid on
            xlabel(['$ Frequency\;\mathrm{[Hz]} $'],'interpreter','latex')
            ylabel(['$ Phase\;\mathrm{[degrees]} $'],'interpreter','latex')
            legend('H1','H2','Hv')
            n=n+1;
    end
end

% Multiple Coherence
M = zeros(5,1,1024);
for p=1:5
    for s=1:3
        for t=1:3
                M(p,1,:) = M(p,1,:) + (H1(p,s,:).*Gff(s,t,:).*conj(H1(p,t,:)))./Gxx(p,p,:);
        end
    end
end
for p=1:5
    figure(n)
    plot(0:511,squeeze(M(p,1,1:Ntime/2)))
    grid on
    title(sprintf('Case 1: Multiple Coherence for Output %d',p))
    xlabel(['$ Frequency\;\mathrm{[Hz]} $'],'interpreter','latex')
    ylabel(['$ \gamma^2\;\mathrm{} $'],'interpreter','latex')
    n = n+1;
end

% Virtual Force
for fi=1:3
    for count=1:Ntime/2
        [V1 D1] = eig([Gff(:,:,count)]);
        vf(fi,count) = D1(fi,fi);
    end
end
for fi=1:3
    figure(n)
    semilogy(0:511,vf(fi,:))
    grid on
    title(sprintf('Case 1: Virtual Force'))
    xlabel(['$ Frequency\;\mathrm{[Hz]} $'],'interpreter','latex')
    ylabel(['$ Force\;\mathrm{[lbf]} $'],'interpreter','latex')
    hold on
end
legend('Input 1','Input 2','Input 3')


% %========================Case 2========================
% %======================================================
Ntime = 1024;
Nc = 4;
Navg = round(length(tt)/(Ntime*Nc))-1;
ff=[];
yy=[];
Gxf1=[]; Gxx1=[]; Gff1=[]; Gfx1=[];
Gxf=[]; Gxx=[]; Gff=[]; Gfx=[];
H1=[]; H2=[]; Hv=[]; V=[]; D=[]; V1=[]; D1=[];
ff = ydata(1:3,:)/0.5;
yy = ydata(4:8,:)/0.1;
h = hann(Ntime*Nc)';

% 50% Overlap + Windowing
for yi=1:5
    for fi=1:3
        for i=1:Navg*2-1
            ff_w{fi}{i} = h.*ff(fi,(1+(i-1)*(Ntime*Nc)/2):((i+1)*(Ntime*Nc)/2));
            yy_w{yi}{i} = h.*yy(yi,(1+(i-1)*(Ntime*Nc)/2):((i+1)*(Ntime*Nc)/2));
        end
    end
end

% Cyclic Averaging (Chop and Sum)
for fi=1:3
    for i=1:Navg*2-1
        ff_w{fi}{i} = (ff_w{fi}{i}(1:(Ntime*Nc)/4) + ff_w{fi}{i}((Ntime*Nc)/4+1:(Ntime*Nc)/4*2) + ff_w{fi}{i}(1+(Ntime*Nc)/2:(Ntime*Nc)/4*3) + ff_w{fi}{i}((Ntime*Nc)/4*3+1:(Ntime*Nc)))/4;
    end
    end
for yi=1:5
    for i=1:Navg*2-1
        yy_w{yi}{i} = (yy_w{yi}{i}(1:(Ntime*Nc)/4) + yy_w{yi}{i}((Ntime*Nc)/4+1:(Ntime*Nc)/4*2) + yy_w{yi}{i}(1+(Ntime*Nc)/2:(Ntime*Nc)/4*3) + yy_w{yi}{i}((Ntime*Nc)/4*3+1:(Ntime*Nc)))/4;
    end
end

% FFT & Power Spectrum Averaging
for fi=1:3
    for i=1:Navg*2-1
        ff_fft{fi}{i} = fft(ff_w{fi}{i});
    end
end
for yi=1:5
    for i=1:Navg*2-1
        yy_fft{yi}{i} = fft(yy_w{yi}{i});
    end
end
for yi=1:5
    for fi=1:3
        for i=1:Navg*2-1
            Gxf1(yi,fi,i,:) = yy_fft{1,yi}{1,i}.*conj(ff_fft{1,fi}{1,i});
        end
    end
end
for fi=1:3
    for yi=1:5
        for i=1:Navg*2-1
            Gfx1(fi,yi,i,:) = ff_fft{1,fi}{1,i}.*conj(yy_fft{1,yi}{1,i});
        end
    end
end
for yi=1:5
    for yii=1:5
        for i=1:Navg*2-1
            Gxx1(yi,yii,i,:) = yy_fft{1,yi}{1,i}.*conj(yy_fft{1,yii}{1,i});
        end
    end
end
for fi=1:3
    for fii=1:3
        for i=1:Navg*2-1
            Gff1(fi,fii,i,:) = ff_fft{1,fi}{1,i}.*conj(ff_fft{1,fii}{1,i});
        end
    end
end

Gxf = mean(Gxf1,3);
Gff = mean(Gff1,3);
Gfx = mean(Gfx1,3);
Gxx = mean(Gxx1,3);

% H1, H2 Algorithm
for i = 1:Ntime
    H1(:,:,1,i) = Gxf(:,:,1,i)/Gff(:,:,1,i);
    H2(:,:,1,i) = Gxx(:,:,1,i)/Gfx(:,:,1,i);    
end
% Hv Algorithm
for yi=1:5
    for fi=1:3
        for cnt_freq = 1:Ntime/2
            [V D] = eig([Gff(:,:,1,cnt_freq).' Gxf(yi,:,1,cnt_freq).';Gfx(:,yi,1,cnt_freq).' Gxx(yi,yi,1,cnt_freq)]);
            eval = diag(D);
            seig = eval(1);
            vec = V(:,1);
            nv=1;
            for k=1:length(eval)
                if(eval(k)<seig)
                    seig=eval(k);
                    vec=V(:,k);
                    nv=k;
                end
            end
            vec = vec*(-1/vec(end,nv));
            Hv1(yi,:,1,cnt_freq) = vec(1:3,nv);
            Hv = permute(Hv1,[2 1 3 4]);
        end
    end
end

% FRFs 
m=22;
for yi=1:5
    for fi=1:3
            figure(m)
            subplot(2,1,1)
            semilogy(0:511,abs(squeeze(H1(yi,fi,1,1:Ntime/2))),'-r')
            hold on
            semilogy(0:511,abs(squeeze(H2(yi,fi,1,1:Ntime/2))),'-b')
            hold on
            semilogy(0:511,abs(squeeze(Hv(fi,yi,1,1:Ntime/2))),'--y')
            grid on
            title(sprintf('Case 2: FRF for H_q_%d_,_p_%d',yi,fi))
            xlabel(['$ Frequency\;\mathrm{[Hz]} $'],'interpreter','latex')
            ylabel(['$ Amplitude\;\mathrm{[g/lbf]} $'],'interpreter','latex')
            legend('H1','H2','Hv')
            subplot(2,1,2)
            plot(0:511,angle(squeeze(H1(yi,fi,1,1:Ntime/2)))*180/pi,'-r')
            hold on
            plot(0:511,angle(squeeze(H2(yi,fi,1,1:Ntime/2)))*180/pi,'-b')
            hold on
            plot(0:511,angle(squeeze(Hv(fi,yi,1,1:Ntime/2)))*180/pi,'--y')
            grid on
            xlabel(['$ Frequency\;\mathrm{[Hz]} $'],'interpreter','latex')
            ylabel(['$ Phase\;\mathrm{[degrees]} $'],'interpreter','latex')
            legend('H1','H2','Hv')
            m=m+1;
    end
end

% Multiple Coherence
M = zeros(5,1,1024);
for p=1:5
    for s=1:3
        for t=1:3
                M(p,1,:) = M(p,1,:) + (H1(p,s,:).*Gff(s,t,:).*conj(H1(p,t,:)))./Gxx(p,p,:);
        end
    end
end
for p=1:5
    figure(m)
    plot(0:511,squeeze(M(p,1,1:Ntime/2)))
    grid on
    title(sprintf('Case 2: Multiple Coherence for Output %d',p))
    xlabel(['$ Frequency\;\mathrm{[Hz]} $'],'interpreter','latex')
    ylabel(['$ \gamma^2\;\mathrm{} $'],'interpreter','latex')
    m = m+1;
end

% Virtual Force
for fi=1:3
    for count=1:Ntime/2
        [V1 D1] = eig([Gff(:,:,count)]);
        vf(fi,count) = D1(fi,fi);
    end
end
for fi=1:3
    figure(m)
    semilogy(0:511,vf(fi,:))
    grid on
    title(sprintf('Case 2: Virtual Force'))
    xlabel(['$ Frequency\;\mathrm{[Hz]} $'],'interpreter','latex')
    ylabel(['$ Force\;\mathrm{[lbf]} $'],'interpreter','latex')
    hold on
end
legend('Input 1','Input 2','Input 3')

%========================Case 3========================
%======================================================
Ntime = 4096;
Nc = 4;
Navg = round(length(tt)/(Ntime*Nc))-1;
ff=[];
yy=[];
ff_w={};
yy_w={};
Gxf1=[]; Gxx1=[]; Gff1=[]; Gfx1=[];
Gxf=[]; Gxx=[]; Gff=[]; Gfx=[];
H1=[]; H2=[]; Hv=[]; V=[]; D=[]; V1=[]; D1=[];
ff = ydata(1:3,:)/0.5;
yy = ydata(4:8,:)/0.1;
h = hann(Ntime*Nc)';

% 50% Overlap + Windowing
for yi=1:5
    for fi=1:3
        for i=1:Navg*2-1
            ff_w{fi}{i} = h.*ff(fi,(1+(i-1)*(Ntime*Nc)/2):((i+1)*(Ntime*Nc)/2));
            yy_w{yi}{i} = h.*yy(yi,(1+(i-1)*(Ntime*Nc)/2):((i+1)*(Ntime*Nc)/2));
        end
    end
end

% Cyclic Averaging (Chop and Sum)
for fi=1:3
    for i=1:Navg*2-1
        ff_w{fi}{i} = (ff_w{fi}{i}(1:(Ntime*Nc)/4) + ff_w{fi}{i}((Ntime*Nc)/4+1:(Ntime*Nc)/4*2) + ff_w{fi}{i}(1+(Ntime*Nc)/2:(Ntime*Nc)/4*3) + ff_w{fi}{i}((Ntime*Nc)/4*3+1:(Ntime*Nc)))/4;
    end
end
for yi=1:5
    for i=1:Navg*2-1
        yy_w{yi}{i} = (yy_w{yi}{i}(1:(Ntime*Nc)/4) + yy_w{yi}{i}((Ntime*Nc)/4+1:(Ntime*Nc)/4*2) + yy_w{yi}{i}(1+(Ntime*Nc)/2:(Ntime*Nc)/4*3) + yy_w{yi}{i}((Ntime*Nc)/4*3+1:(Ntime*Nc)))/4;
    end
end

% FFT & Power Spectrum Averaging
for fi=1:3
    for i=1:Navg*2-1
        ff_fft{fi}{i} = fft(ff_w{fi}{i});
    end
end
for yi=1:5
    for i=1:Navg*2-1
        yy_fft{yi}{i} = fft(yy_w{yi}{i});
    end
end
for yi=1:5
    for fi=1:3
        for i=1:Navg*2-1
            Gxf1(yi,fi,i,:) = yy_fft{1,yi}{1,i}.*conj(ff_fft{1,fi}{1,i});
        end
    end
end
for fi=1:3
    for yi=1:5
        for i=1:Navg*2-1
            Gfx1(fi,yi,i,:) = ff_fft{1,fi}{1,i}.*conj(yy_fft{1,yi}{1,i});
        end
    end
end
for yi=1:5
    for yii=1:5
        for i=1:Navg*2-1
            Gxx1(yi,yii,i,:) = yy_fft{1,yi}{1,i}.*conj(yy_fft{1,yii}{1,i});
        end
    end
end
for fi=1:3
    for fii=1:3
        for i=1:Navg*2-1
            Gff1(fi,fii,i,:) = ff_fft{1,fi}{1,i}.*conj(ff_fft{1,fii}{1,i});
        end
    end
end

Gxf = mean(Gxf1,3);
Gff = mean(Gff1,3);
Gfx = mean(Gfx1,3);
Gxx = mean(Gxx1,3);

% H1, H2 Algorithm
for i = 1:Ntime
    H1(:,:,1,i) = Gxf(:,:,1,i)/Gff(:,:,1,i);
    H2(:,:,1,i) = Gxx(:,:,1,i)/Gfx(:,:,1,i);    
end
% HV Algorithm
for yi=1:5
    for fi=1:3
        for cnt_freq = 1:Ntime/2
            [V D] = eig([Gff(:,:,1,cnt_freq).' Gxf(yi,:,1,cnt_freq).';Gfx(:,yi,1,cnt_freq).' Gxx(yi,yi,1,cnt_freq)]);
            eval = diag(D);
            seig = eval(1);
            vec = V(:,1);
            nv=1;
            for k=1:length(eval)
                if(eval(k)<seig)
                    seig=eval(k);
                    vec=V(:,k);
                    nv=k;
                end
            end
            vec = vec*(-1/vec(end,nv));
            Hv1(yi,:,1,cnt_freq) = vec(1:3,nv);
            Hv = permute(Hv1,[2 1 3 4]);
        end
    end
end

% FRFs 
b=43;
for yi=1:5
    for fi=1:3
            figure(b)
            subplot(2,1,1)
            semilogy(0:2047,abs(squeeze(H1(yi,fi,1,1:Ntime/2))),'-r')
            hold on
            semilogy(0:2047,abs(squeeze(H2(yi,fi,1,1:Ntime/2))),'-b')
            hold on
            semilogy(0:2047,abs(squeeze(Hv(fi,yi,1,1:Ntime/2))),'--y')
            grid on
            title(sprintf('Case 3: FRF for H_q_%d_,_p_%d',yi,fi))
            xlabel(['$ Frequency\;\mathrm{[Hz]} $'],'interpreter','latex')
            ylabel(['$ Amplitude\;\mathrm{[g/lbf]} $'],'interpreter','latex')
            legend('H1','H2','Hv')
            subplot(2,1,2)
            plot(0:2047,angle(squeeze(H1(yi,fi,1,1:Ntime/2)))*180/pi,'-r')
            hold on
            plot(0:2047,angle(squeeze(H2(yi,fi,1,1:Ntime/2)))*180/pi,'-b')
            hold on
            plot(0:2047,angle(squeeze(Hv(fi,yi,1,1:Ntime/2)))*180/pi,'--y')
            grid on
            xlabel(['$ Frequency\;\mathrm{[Hz]} $'],'interpreter','latex')
            ylabel(['$ Phase\;\mathrm{[degrees]} $'],'interpreter','latex')
            legend('H1','H2','Hv')
            b=b+1;
    end
end

% Multiple Coherence
M = zeros(5,1,4096);
for p=1:5
    for s=1:3
        for t=1:3
                M(p,1,:) = M(p,1,:) + (H1(p,s,:).*Gff(s,t,:).*conj(H1(p,t,:)))./Gxx(p,p,:);
        end
    end
end
for p=1:5
    figure(b)
    plot(0:2047,squeeze(M(p,1,1:Ntime/2)))
    grid on
    title(sprintf('Case 3: Multiple Coherence for Output %d',p))
    xlabel(['$ Frequency\;\mathrm{[Hz]} $'],'interpreter','latex')
    ylabel(['$ \gamma^2\;\mathrm{} $'],'interpreter','latex')
    b = b+1;
end

% Virtual Force
for fi=1:3
    for count=1:Ntime/2
        [V1 D1] = eig([Gff(:,:,count)]);        % Eigenvalues and eigenvector matrix
        vf(fi,count) = D1(fi,fi);
    end
end
for fi=1:3
    figure(b)
    semilogy(0:2047,vf(fi,:))
    grid on
    title(sprintf('Case 3: Virtual Force'))
    xlabel(['$ Frequency\;\mathrm{[Hz]} $'],'interpreter','latex')
    ylabel(['$ Force\;\mathrm{[lbf]} $'],'interpreter','latex')
    hold on
end
legend('Input 1','Input 2','Input 3')

%====================================================================
%====================================================================
%============================END OF PROGRAM==========================
%====================================================================
%====================================================================


