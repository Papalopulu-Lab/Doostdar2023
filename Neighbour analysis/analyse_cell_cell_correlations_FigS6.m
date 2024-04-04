clear all, close all, clc
% HV examples over time
num=xlsread('190822 HV ALL spots based on Keima_Statistics\Position_INTENSITY.csv');
X=num(:,1);
Y=num(:,2);
Z=num(:,3);
T=num(:,7);
Venus=num(:,9); % 9 is venus 
%
HILL=[];
LIST1=[];
LIST2=[];
n=0;
vect=[1,20,40,60,80];
for i=1:numel(vect)
    HILL=[];
    idx=find(T==vect(i));
    x=X(idx,:); y=Y(idx,:); z=Z(idx,:);
    venus=Venus(idx,:);
    [~,dist]=getCellCellDist(x,y,z);
    % compute cell pairs by taking the closest neighbour
    rec_hill(1:numel(x),1)=venus;
    for j=1:numel(x)
        vals=dist(j,:); 
        nidx=find(vals==min(vals),1,'first');
        rec_hill(j,2)=venus(nidx);
        dn(i,j)=vals(nidx);
    end
    HILL=[HILL; rec_hill];
    figure(1),subplot(2,5,i)
    plot(HILL(:,1),HILL(:,2),'o')
    hold on
    xlabel('Cell 1 (Fl Intensity)');
    ylabel('Cell 2 (Fl Intensity)');
    if i==1
        title(['HV: Time = 20h']); 
    else
        title(['HV: Time = ', num2str(vect(i)/10+20) 'h']); 
    end
    % use this for venus
    set(gca,'XLim',[0,5000],'YLim',[0,5000]);
    axis square
    % compute correlation coefficient
    rc=corr(HILL(:,1),HILL(:,2));
    text(3500,4500,['r= ', num2str(round(rc*100)/100)]);
    % compute intensity ratio
    M=mean(HILL(:,1)./HILL(:,2));
    text(500,4500,['ratio= ' num2str(round(M*100)/100)]);
    clear rec_hill HILL dist dn
end

%% HVP examples over time
num=xlsread('190822 HVP ALL spots based on Keima_Statistics\Position_INTENSITY.csv');
X=num(:,1);
Y=num(:,2);
Z=num(:,3);
T=num(:,7);
Venus=num(:,9); % 9 is venus 
%
HILL=[];
LIST1=[];
LIST2=[];
n=0;
vect=[1,20,40,60,80];
for i=1:numel(vect)
    HILL=[];
    idx=find(T==vect(i));
    x=X(idx,:); y=Y(idx,:); z=Z(idx,:);
    venus=Venus(idx,:);
    [~,dist]=getCellCellDist(x,y,z);
    % compute cell pairs by taking the closest neighbour
    rec_hill(1:numel(x),1)=venus;
    for j=1:numel(x)
        vals=dist(j,:); 
        nidx=find(vals==min(vals),1,'first');
        rec_hill(j,2)=venus(nidx);
        dn(i,j)=vals(nidx);
    end
    HILL=[HILL; rec_hill];
    figure(1),subplot(2,5,i+5)
    plot(HILL(:,1),HILL(:,2),'o')
    hold on
    xlabel('Cell 1 (Fl Intensity)');
    ylabel('Cell 2 (Fl Intensity)');
    if i==1
        title(['HVP: Time = 20h']); 
    else
        title(['HVP: Time = ', num2str(vect(i)/10+20) 'h']); 
    end
    % use this for venus
    set(gca,'XLim',[0,5000],'YLim',[0,5000]);
    axis square
    % compute correlation coefficient
    rc=corr(HILL(:,1),HILL(:,2));
    text(3500,4500,['r= ', num2str(round(rc*100)/100)]);
    % compute intensity ratio
    M=mean(HILL(:,1)./HILL(:,2));
    text(500,4500,['ratio= ' num2str(round(M*100)/100)]);
    clear rec_hill HILL dist dn
end
%% HV mKeima examples over time
num=xlsread('190822 HV ALL spots based on Keima_Statistics\Position_INTENSITY.csv');
X=num(:,1);
Y=num(:,2);
Z=num(:,3);
T=num(:,7);
mKeima=num(:,10); % 10 is keima-nuclear marker 
%
HILL=[];
LIST1=[];
LIST2=[];
n=0;
vect=[1,20,40,60,80];
for i=1:numel(vect)
    HILL=[];
    idx=find(T==vect(i));
    x=X(idx,:); y=Y(idx,:); z=Z(idx,:);
    mkeima=mKeima(idx,:);
    [~,dist]=getCellCellDist(x,y,z);
    % compute cell pairs by taking the closest neighbour
    rec_hill(1:numel(x),1)=mkeima;
    for j=1:numel(x)
        vals=dist(j,:); 
        nidx=find(vals==min(vals),1,'first');
        rec_hill(j,2)=mkeima(nidx);
        dn(i,j)=vals(nidx);
    end
    HILL=[HILL; rec_hill];
    figure(2),subplot(2,5,i)
    plot(HILL(:,1),HILL(:,2),'o')
    hold on
    xlabel('Cell 1 (Fl Intensity)');
    ylabel('Cell 2 (Fl Intensity)');
    if i==1
        title(['HV-mKeima: Time = 20h']); 
    else
        title(['HV-mKeima: Time = ', num2str(vect(i)/10+20) 'h']); 
    end
    % use this for venus
    set(gca,'XLim',[0,6000],'YLim',[0,6000]);
    axis square
    % compute correlation coefficient
    r=corrcoef(HILL(:,1),HILL(:,2))
    r=HILL(:,1)./HILL(:,2);
    M(i)=mean(r);
    % compute correlation coefficient
    rc=corr(HILL(:,1),HILL(:,2));
    text(4500,5500,['r= ', num2str(round(rc*100)/100)]);
    % compute intensity ratio
    M=mean(HILL(:,1)./HILL(:,2));
    text(1000,5500,['ratio= ' num2str(round(M*100)/100)]);
    clear rec_hill HILL dist dn
end 
%% HVP mKeima examples over time
num=xlsread('190822 HVP ALL spots based on Keima_Statistics\Position_INTENSITY.csv');
X=num(:,1);
Y=num(:,2);
Z=num(:,3);
T=num(:,7);
mKeima=num(:,10); % 10 is keima-nuclear marker 
%
HILL=[];
LIST1=[];
LIST2=[];
n=0;
vect=[1,20,40,60,80];
for i=1:numel(vect)
    HILL=[];
    idx=find(T==vect(i));
    x=X(idx,:); y=Y(idx,:); z=Z(idx,:);
    mkeima=mKeima(idx,:);
    [~,dist]=getCellCellDist(x,y,z);
    % compute cell pairs by taking the closest neighbour
    rec_hill(1:numel(x),1)=mkeima;
    for j=1:numel(x)
        vals=dist(j,:); 
        nidx=find(vals==min(vals),1,'first');
        rec_hill(j,2)=mkeima(nidx);
        dn(i,j)=vals(nidx);
    end
    HILL=[HILL; rec_hill];
    figure(2),subplot(2,5,i+5)
    plot(HILL(:,1),HILL(:,2),'o')
    hold on
    xlabel('Cell 1 (Fl Intensity)');
    ylabel('Cell 2 (Fl Intensity)');
    if i==1
        title(['HVP-mKeima: Time = 20h']); 
    else
        title(['HVP-mKeima: Time = ', num2str(vect(i)/10+20) 'h']); 
    end
    % use this for venus
    set(gca,'XLim',[0,3000],'YLim',[0,3000]);
    axis square
    % compute correlation coefficient
    rc=corr(HILL(:,1),HILL(:,2));
    text(2300,2500,['r= ', num2str(round(rc*100)/100)]);
    % compute intensity ratio
    M=mean(HILL(:,1)./HILL(:,2));
    text(300,2500,['ratio= ' num2str(round(M*100)/100)]);
    clear rec_hill HILL dist dn
end 
