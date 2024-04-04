function h=showfigure1(x,m,raw,y1,xd,md,M1,S1,S2,LLR,par2,i)
% adjust to fl intensity
m=m*S1;%+M1;
raw=raw*S1;%+M1;
y1=y1*S2*S1;
md=md*S2*S1;
%SHOWFIGURE Summary of this function goes here
%   Detailed explanation goes here
    % plots raw data with trend and detrended data
    h=figure();
    subplot(2,1,1)
    hold on
    plot(x,m,'r','LineWidth',2)
    plot(x,raw,'r');%,'+')
    xlabel('Time (hours)')
    ylabel('Fluorescence (a.u.)')
    str = sprintf('Raw fluorescent intensity cell %.0f',i);
    title(str,'fontweight','normal');
    %set(gca,'XLim',[0,12]);% was72[15,50]
    %set(gca,'YLim',[1000,3000]); % was 5000
    %text(x(2),4300,['Amplitude ', num2str(round(2*S1))],'FontSize',12);
    subplot(2,1,2)    
    plot(x,y1)
    hold on
    plot(xd,md,'b','LineWidth',2)
    xlabel('Time (hours)')
    %set(gca,'XLim',[0,12]);
    %text(x(2),750,['Amplitude ', num2str(round(2*S2*S1))],'FontSize',12);
    %set(gca,'YLim',[-1000,1000]);% was 1000
    ylabel('Detrended fluorescence (a.u.)');
    if LLR<0
    str = sprintf('LLR score %f',LLR);
    title(str);
    else
        q=par2(2)/(par2(1));
        %snr=par2(3)/par2(4);
        %q=q*snr;
        if q>10^8;
            q=Inf;
        end
    str = sprintf('LLR score %f, Period = %f, Quality = %f',LLR,(2*pi()/par2(2)),q);
    title(str,'fontweight','normal');     
   
end

