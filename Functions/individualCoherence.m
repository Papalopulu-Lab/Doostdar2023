%individualCoherence.m Instead of calculating the coherence value of the
%average Fourier transform for a population of cells, this will calculate
%the coherence value of each individual cell's time trace and then return
%the average individual cell coherence. This is an important distinction
%because the population coherence will be affected by spread in oscillation
%frequency.

function [Ydetrend, t, Ysmooth,f,P1,averageIndividualCoherence,f_C1,P_C1, ind_per, FisherPer, Occurrence, coherenceIndividualCells]=individualCoherence(T,Y_raw,polyOrder,frac,frameLength,Nt,dt)
    
    [Ydetrend, t, Ysmooth]=detrendSgolay(Y_raw, polyOrder, 0, T, frameLength);

    if frac==0
        startIDX=1;
    else
        startIDX=frac*Nt;
    end
    
    Ydetrend=Ydetrend(:,startIDX:end);
  
    Ysmooth=Ysmooth(:,startIDX:end);
    t=t(startIDX:end);
    
    %Fast Fourier transform (FT) of the detrended signal
    Y=fft(Ydetrend,[],2);
    L=length(Y);
    Fs=1/dt;
    f = Fs*(0:(L/2))/L;
    P2 = abs(Y/L).^2;
    P1 = P2(:,1:L/2+1);
    P1(:,2:end-1) = 2*P1(:,2:end-1);

    [FisherPer,Occurrence] = Fisher(P1',f);

    FisherPer = dt*FisherPer/60;


    [~,i]=max(P1,[],2);
    ind_per=1./(60*f(i)); %Individual cell fourier dominant periods

    numCells = size(Y_raw,1);

    for cell = 1:numCells

        cellFourier = P1(cell, :);
    
        [~,I]=max(cellFourier); %Peak value index of the individual cell's FT (largest frequency contribution to signal)
        
        df=f(2)-f(1);
        tenP=0.1*f(I);
        tenP_idx=0.1*f(I)/df; %Ten percent in index value of peak frequency value in power spectrum, used in coherence calculation
    
        if I==1
            C1=0;
            f_C1=0;
            P_C1=cellFourier(I);
            
        elseif round(tenP_idx)==tenP_idx
            C1=trapz( f(I-tenP_idx:I+tenP_idx),cellFourier(I-tenP_idx:I+tenP_idx) );
            f_C1=f(I-tenP_idx:I+tenP_idx);
            P_C1=cellFourier(I-tenP_idx:I+tenP_idx);
            
        else
            [f_C1, P_C1]=accurateC1(f,cellFourier,I,tenP,tenP_idx); %Interpolates between fourier values to get more accurate coherence area  
            C1=trapz(f_C1,P_C1);
        end
        
        C2=trapz(f,cellFourier);
        coherence(cell)=C1/C2;

    end
coherenceIndividualCells = coherence;
averageIndividualCoherence = mean(coherence);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f_C1, P_C1]=accurateC1(f,P,I,tenP,tenP_idx)
    

    idx_im=I-floor(tenP_idx);
    idx_om=I-ceil(tenP_idx);
    idx_ip=I+floor(tenP_idx);
    idx_op=I+ceil(tenP_idx);
    
    if idx_im < 1 || idx_om < 1
        idx_im = 1;
        idx_om = 1;
    end
    
    if idx_ip > length(f) || idx_op > length(f)
        idx_ip = length(f);
        idx_op = length(f);
    end
    
    f_im=f(idx_im); %Frequency value at inner, minor (Subscripts: i=inner, m=minor)
    f_ip=f(idx_ip); %Frequency value at inner, plus (Subscripts: i=inner, p=plus)
    f_om=f(idx_om);
    f_op=f(idx_op);
    
    P_im=P(idx_im);
    P_ip=P(idx_ip); 
    P_om=P(idx_om);
    P_op=P(idx_op);
    
    f_mm=f(I)-tenP;
    f_mp=f(I)+tenP; 
    
    smidge_m=abs((f_mm-f_im))/abs((f_om-f_im));
    smidge_p=abs((f_mp-f_ip))/abs((f_op-f_ip));
    P_mm=P_om*smidge_m + P_im*(1-smidge_m);
    P_mp=P_op*smidge_p + P_ip*(1-smidge_p);
    
    f_C1=[f_mm, f(idx_im:idx_ip), f_mp];
    P_C1=[P_mm,  P(idx_im:idx_ip), P_mp ];
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [FisherPer,Occurrence]=Fisher(P1,f)

II=size(P1,2);

    PXX2=P1(2:end,:);
    F2=f(2:end);
    %Preallocation for storing stats
    fisherG=zeros(II,1);
    pval=zeros(II,1);
    peakLoc=zeros(II,1);
    peakHeight=zeros(II,1);
    FisherPer=zeros(II,1);

    for ii=1:II

        pxx=PXX2(:,ii); % periodogram
        
        [fisherG(ii),pval(ii),idx]=GetFisherG(pxx); % Find the peak and collect Fisher G-statistic
        peakLoc(ii)=idx;
        peakHeight(ii)=pxx(idx);
        if pval(ii)<0.05
            FisherPer(ii)=1/F2(idx);
        else
            FisherPer(ii)=nan;
        end
    end
    
    Occurrence=1-sum(isnan(FisherPer))/II;
    
end

%__________________________________________________________________________

%GetFisherG is based on the matlab page and refs found here: 
%https://www.mathworks.com/help/signal/ug/significance-testing-for-periodic-component.html

function [fisher_g,pval,idx]=GetFisherG(Pxx)
    idx=find(Pxx==max(Pxx),1,'first');
    fisher_g=Pxx(idx)/sum(Pxx);
    N = length(Pxx);
    nn = 1:floor(1/fisher_g);

    I = (-1).^(nn-1).*exp(gammaln(N+1)-gammaln(nn+1)-gammaln(N-nn+1)).*(1-nn*fisher_g).^(N-1);
    
    pval = sum(I);
end