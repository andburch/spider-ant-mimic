clear all; close all
cd 'C:/Users/Andrew Burchill/Documents/Side Projects/spider-ant-mimic/sync' % Set directory
dinfo = dir;    %opens the directory
A = {dinfo.name};   %gets names of files in the directory

%A=sort_nat(A)   %sorts them
j=1;
for m=3:length(A)   %goes through every file... except three?
%figure
filed=[char(A(m))];
 activityt4=csvread(filed);     %reads the file
 
 
 time=(1:length(activityt4))/240;  %converts index to hours
                                            %1 sec = 0.0002778 hrs
                                            %each time pt 30 sec apart
 total_time = time(length(time));

 % uncomment next line to normalize time series 

 % activityt4=rescale(activityt4);
 
% plot(time,activityt4,'k','LineWidth',5)
% xlabel('Seconds')
% set(gca,'fontsize',15)
% set(gca,'linewidth',2)
% set(gca,'TickDir','out')
% set(gca, 'FontName', 'Times')


% figure;
% cwt(activityt4,seconds(1/240));   %plots the wavelet thing
% xlabel(A(m))



%gets the wavelet transform, the sampling freq in Hz, & cone of influence
[dpoaeCWT,f,coi] = cwt(activityt4,240);
%breaks apart 0-9hrs by the number of the wavelet transforms
%t(1)=0, unlike the `time` variable
t=linspace(0,total_time,length(dpoaeCWT))';
%get the lefthand unique cone of infl. values 
[coil, index] = unique(coi); 
%coil=x, t(index)=y, f= sampling points
YIsx1 = interp1(coil, t(index), f);

midway = floor(length(dpoaeCWT)/2)
%repeat for the righthand side. (WHY 534?)
[coir, index] = unique(coi(midway:end)); 
%repeat, but add 4.5 hrs to time
YIsx2 = interp1(coir, t(index)+total_time/2, f);

%for each wavelet transform (freq)
for i=1:length(f)
    cfsOAE = dpoaeCWT(i,:); %get its fourier series?
    q=abs(cfsOAE);      
    h=vertcat(cfsOAE,t');   %two row matrix, time is second row
    
    %only get time values within/outside cone of influence
    q2 = h(2,:) < YIsx1(i) | h(2,:) > YIsx2(i)
    h2=h;   %not sure why this is called....
    q(q2) = [];     %only keep the "good" values of q?
    power(i)=max(q);    %find the maximum value for the wt
end
[M,Y] = max(power)
period(j)=1/f(Y)


rhythmicity(j)=M
sync(j)=[(std(activityt4)).^2/mean(activityt4)]
id(j)=A(m)

power=[];
j=j+1;
end

ant_time_series=table(id',sync',rhythmicity',period')
ant_time_series.Properties.VariableNames ={'Colony' 'Synchrony' 'Rhythmicity' 'Period'}
writetable(ant_time_series)