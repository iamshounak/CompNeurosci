close all;
load('data_cn_project_iii_a17.mat');

%% DATA SEPARATION
train_stm=Stimulus(1:15000);
test_stm=Stimulus(15001:20000);

for i=1:4
    for j=1:50
        train_spikes{i,j}=[];
        test_spikes{i,j}=[];
        curr_spike_times=All_Spike_Times{i,j};
        for k=1:length(curr_spike_times)
            if curr_spike_times(1,k)<15
                train_spikes{i,j}=[train_spikes{i,j},curr_spike_times(1,k)];
            else
                test_spikes{i,j}=[test_spikes{i,j},curr_spike_times(1,k)-15]; 
                %Spike timings considered from time after 15 seconds for test data
            end
        end
    end
end

%% PART1
max_delay=50;
[acf,lags]=xcorr(Stimulus,max_delay,'normalized');
figure
stem(lags,acf)
xlabel("Time-delay");
ylabel("ACF Value");
title("Autocorrelation function of stimulus");
%PLOT LOOKS LIKE DIRAC-DELTA FUNCTION

%% PART 2
binsize=0.001;
countspikes_train=zeros(4,20/binsize);
for i=1:4
    for j=1:50
        curr_spike_times_train=train_spikes{i,j};
        for k=1:length(curr_spike_times_train)
            countspikes_train(i,1+floor(curr_spike_times_train(1,k)*1000))=countspikes_train(i,1+floor(curr_spike_times_train(1,k)*1000))+1;
        end            
    end
end

countspikes_test=zeros(4,20/binsize);
for i=1:4
    for j=1:50
        curr_spike_times_test=test_spikes{i,j};
        for k=1:length(curr_spike_times_test)
            countspikes_test(i,1+floor(curr_spike_times_test(1,k)*1000))=countspikes_test(i,1+floor(curr_spike_times_test(1,k)*1000))+1;
        end            
    end
end

figure
subplot(2,2,1);
plot(countspikes_train(1,1:20000))
ylabel("Firing rate");
xlabel("Time");
title("PSTH for training data - neuron 1");
subplot(2,2,2);
plot(countspikes_train(2,1:20000))
ylabel("Firing rate");
xlabel("Time");
title("PSTH for training data - neuron 2");
subplot(2,2,3);
plot(countspikes_train(3,1:20000))
ylabel("Firing rate");
xlabel("Time");
title("PSTH for training data - neuron 3");
subplot(2,2,4);
plot(countspikes_train(4,1:20000))
ylabel("Firing rate");
xlabel("Time");
title("PSTH for training data - neuron 4");


figure
subplot(2,2,1);
plot(countspikes_test(1,1:20000))
ylabel("Firing rate");
xlabel("Time");
title("PSTH for test data - neuron 1");
subplot(2,2,2);
plot(countspikes_test(2,1:20000))
ylabel("Firing rate");
xlabel("Time");
title("PSTH for test data - neuron 2");
subplot(2,2,3);
plot(countspikes_test(3,1:20000))
ylabel("Firing rate");
xlabel("Time");
title("PSTH for test data - neuron 3");
subplot(2,2,4);
plot(countspikes_test(4,1:20000))
ylabel("Firing rate");
xlabel("Time");
title("PSTH for test data - neuron 4");


%% PART 3
binsizes=[0.01,0.02,0.05,0.1,0.2,0.5];
for d=1:length(binsizes)
    countspikes=zeros(4,50,20/binsizes(d));
    figure
    sgtitle(["Mean-variance scatter plot for bin size of",num2str(binsizes(d))]);
    for i=1:4    
        for j=1:50
            curr_spike_times=All_Spike_Times{i,j};
            for k=1:length(curr_spike_times)
                countspikes(i,j,1+floor(curr_spike_times(1,k)*1./binsizes(d)))=countspikes(i,j,1+floor(curr_spike_times(1,k)*1./binsizes(d)))+1;
            end            
        end
        meanval=mean(countspikes(i,:,:),2);
        varval=var(countspikes(i,:,:),1,2);
        subplot(2,2,i);
        scatter(meanval,varval)
        xlabel("mean");
        ylabel("variance");
        title(["Neuron",num2str(i)]);
    end    
end
%Poisson process if mean vs variance curve looks like y=x line because mean
%and variance are equal for Poisson
%On increasing bin size, Poisson property is lost

%% PART 4
windowsize=100;
sta_array=zeros(4,100);
train_stm_part4=train_stm;
train_stm_part4(1:windowsize)=0;
numspike=0;
figure
sgtitle("Spike Triggered Averages for four neurons");
for i=1:4
    avg=zeros(1,100);
    for j=1:50
        curr_spike_times_train=train_spikes{i,j};
        for k=1:length(curr_spike_times_train)
            window=1+floor(curr_spike_times_train(k)*1000)-(windowsize-1):1+floor(curr_spike_times_train(k)*1000);
            if(window(1)>=1)
                avg=avg+train_stm_part4(flip(window)); %Flip so that delays from left to right increasing
            end 
        end
        numspike=numspike+length(curr_spike_times_train);     
    end
    avg=avg./numspike;
    sta_array(i,:)=avg;
    subplot(2,2,i)
    plot(avg)
    xlabel("Time");
    ylabel("STA");
    ylim([-0.06 0.06]);
    title(["Neuron",num2str(i)]);
end

%% CORRECTION FOR NON-GAUSSIANITY
corr_new=zeros(101);
for i = 0:100
    for l=1:20000-i
        corr_new(i+1)=corr_new(i+1)+(Stimulus(l)*Stimulus(l+i));
    end
    corr_new(i+1)=corr_new(i+1)/(20000-i);
end
corr_mat=zeros(100,100);

for row=1:100
    for col=1:100
        corr_mat(row,col)=corr_new(abs(row-col)+1);
    end
end
cssinv=inv(corr_mat);

figure
sgtitle("Spike Triggered Averages for four neurons after non-Gaussianity correction with Css inverse");
new_ar=zeros(100,4);
for i=1:4
    new_ar(:,i)=cssinv*sta_array(i,:)';
    subplot(2,2,i)
    plot(new_ar(:,i))
    xlabel("Time");
    ylabel("STA");
    %ylim([-0.06 0.06]);
    title(["Neuron",num2str(i)]);
end

%NO GOOD RESULT WITH CSS AND CORRECTION AS ALREADY GAUSSIAN NATURE 
%AND ONLY SOME ERROR TERMS COME INTO CSS-INVERSE 
%% PART 5
y_est=zeros(4,15099);
for i=1:4
    y_est(i,:)=conv(sta_array(i,:),train_stm);
end
lambda_est=zeros(4,15000);
lambda_est(:,1:15000)=countspikes_train(:,1:15000);
x=y_est(:,1:15000);
y=lambda_est(:,1:15000);

bin_size = 100; %estimated y and estimated lambda grouped into bins 
%for better visibility of scatter plots and better curve fitting
    for i = 1:ceil(15000/bin_size)
        e = i*bin_size;
        if e>15000
            e=15000;
        end
        x1(i)= mean(x(1,(1+(i-1)*bin_size):e));
        x2(i)= mean(x(2,(1+(i-1)*bin_size):e));
        x3(i)= mean(x(3,(1+(i-1)*bin_size):e));
        x4(i)= mean(x(4,(1+(i-1)*bin_size):e));
    
        y1(i)= mean(y(1,(1+(i-1)*bin_size):e));
        y2(i)= mean(y(2,(1+(i-1)*bin_size):e));
        y3(i)= mean(y(3,(1+(i-1)*bin_size):e));
        y4(i)= mean(y(4,(1+(i-1)*bin_size):e));
    end
    
    figure()
    subplot(2,2,1)
    scatter(x1, y1)
    title('Measured \lambda(t) vs y(t) for neuron 1')
    xlabel('y(t) = convolution of s(t) and h(t)')
    ylabel('\lambda(t)')
    subplot(2,2,2)
    scatter(x2, y2)
    title('Measured \lambda(t) vs y(t) for neuron 2')
    xlabel('y(t) = convolution of s(t) and h(t)')
    ylabel('\lambda(t)')
    subplot(2,2,3)
    scatter(x3, y3)
    title('Measured \lambda(t) vs y(t) for neuron 3')
    xlabel('y(t) = convolution of s(t) and h(t)')
    ylabel('\lambda(t)')
    subplot(2,2,4)
    scatter(x4, y4)
    title('Measured \lambda(t) vs y(t) for neuron 4')
    xlabel('y(t) = convolution of s(t) and h(t)')
    ylabel('\lambda(t)')
    
    poly_order=1;
    %%Sigmoid fitting function approximation is given below:
    fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)));
    figure
    p1=polyfit(x1,y1,poly_order);
    %p1=sigm_fit(x1,y1,[],[],0);
    x_axis = linspace(min(x1),max(x1));
    y_axis = polyval(p1,x_axis);
    %y_axis = fsigm(p1,x_axis);
    subplot(2,2,1)
    plot(x1,y1,'o')
    hold on
    plot(x_axis,y_axis)
    hold off
    title('f(t) estimated from scatter plot between y(t) and \lambda(t)')
    xlabel('y(t) = convolution of s(t) and h(t)')
    ylabel('\lambda(t)')
    
    %p2=sigm_fit(x2,y2,[],[],0);
    p2=polyfit(x2,y2,poly_order);
    x_axis = linspace(min(x2),max(x2));
    y_axis = polyval(p2,x_axis);
    %y_axis = fsigm(p2,x_axis);
    subplot(2,2,2)
    plot(x2,y2,'o')
    hold on
    plot(x_axis,y_axis)
    hold off
    title('f(t) estimated from scatter plot between y(t) and \lambda(t)')
    xlabel('y(t) = convolution of s(t) and h(t)')
    ylabel('\lambda(t)')
    
    p3=polyfit(x3,y3,poly_order);
    %p3=sigm_fit(x3,y3,[],[],0);
    x_axis = linspace(min(x3),max(x3));
    %y_axis = fsigm(p3,x_axis);
    y_axis = polyval(p3,x_axis);
    subplot(2,2,3)
    plot(x3,y3,'o')
    hold on
    plot(x_axis,y_axis)
    hold off
    title('f(t) estimated from scatter plot between y(t) and \lambda(t)')
    xlabel('y(t) = convolution of s(t) and h(t)')
    ylabel('\lambda(t)')
    
    p4=polyfit(x4,y4,poly_order);
    %p4=sigm_fit(x4,y4,[],[],0);
    x_axis = linspace(min(x4),max(x4));
    %y_axis = fsigm(p4,x_axis);
    y_axis = polyval(p4,x_axis);
    subplot(2,2,4)
    plot(x4,y4,'o')
    hold on
    plot(x_axis,y_axis)
    hold off
    title('f(t) estimated from scatter plot between y(t) and \lambda(t)')
    xlabel('y(t) = convolution of s(t) and h(t)')
    ylabel('\lambda(t)')
    
%% PART 6
%PREDICTION OF RESPONSE FOR LAST 5 SECONDS
windowsize=100;
sta_array_test=zeros(4,100);
test_stm_part4=zeros(1,5100);
test_stm_part4(windowsize+1:5100)=test_stm;
test_stm_part4(1:windowsize)=0;
numspike=0;
for i=1:4
    avg=zeros(1,100);
    for j=1:50
        curr_spike_times_train=test_spikes{i,j};
        for k=1:length(curr_spike_times_train)
            window=1+floor(curr_spike_times_train(k)*1000)-(windowsize-1):1+floor(curr_spike_times_train(k)*1000);
            if(window(1)>=1)
                avg=avg+test_stm_part4(flip(window)); %Flip so that delays from left to right increasing
            end 
        end
        numspike=numspike+length(curr_spike_times_train);     
    end
    avg=avg./numspike;
    sta_array_test(i,:)=avg;
end

y_pred=zeros(4,5099);
for i=1:4
   y_pred(i,:)=conv(sta_array_test(i,:),test_stm);
end

%lambda_est=zeros(4,5000);
fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)));    
x=countspikes_test(:,1:5000);
y=y_pred(:,1:5000);
bin_size = 100; %estimated y and estimated lambda grouped into bins 
%for better visibility of scatter plots and better curve fitting
    for i = 1:ceil(5000/bin_size)
        e = i*bin_size;
        if e>5000
            e=5000;
        end
        x1(i)= mean(x(1,(1+(i-1)*bin_size):e));
        x2(i)= mean(x(2,(1+(i-1)*bin_size):e));
        x3(i)= mean(x(3,(1+(i-1)*bin_size):e));
        x4(i)= mean(x(4,(1+(i-1)*bin_size):e));
    
        y1(i)= mean(y(1,(1+(i-1)*bin_size):e));
        y2(i)= mean(y(2,(1+(i-1)*bin_size):e));
        y3(i)= mean(y(3,(1+(i-1)*bin_size):e));
        y4(i)= mean(y(4,(1+(i-1)*bin_size):e));
    end
    %y1=fsigm(p1,y1); %p1,p2,p3,p4 are parameters learnt from training data
    %y2=fsigm(p2,y2);
    %y3=fsigm(p3,y3);
    %y4=fsigm(p4,y4);
    y1=polyval(p1,y1); %p1,p2,p3,p4 are parameters learnt from training data
    y2=polyval(p2,y2);
    y3=polyval(p3,y3);
    y4=polyval(p4,y4);
    figure()
    subplot(2,2,1)
    scatter(x1, y1)
    title('Measured vs predicted PSTH for neuron 1')
    xlabel('Measured \lambda(t)')
    ylabel('Predicted \lambda(t)')
    subplot(2,2,2)
    scatter(x2, y2)
    title('Measured vs predicted PSTH for neuron 2')
    xlabel('Measured \lambda(t)')
    ylabel('Predicted \lambda(t)')
    subplot(2,2,3)
    scatter(x3, y3)
    title('Measured vs predicted PSTH for neuron 3')
    xlabel('Measured \lambda(t)')
    ylabel('Predicted \lambda(t)')
    subplot(2,2,4)
    scatter(x4, y4)
    title('Measured vs predicted PSTH for neuron 4')
    xlabel('Measured \lambda(t)')
    ylabel('Predicted \lambda(t)')

    
%FROM THE PLOTS FOR BOTH SIGMOID AND LINEAR CURVE FITTINGS, WE FIND THAT A
%PART OF THE DATA HAS MEASURED AND PREDICTED VALUES CORRESPONDING VERY WELL
%WHILE ANOTHER PORTION OF THE DATA HAS MEASURED AND PREDICTED VALUES
%SIGNIFICANTLY DIFFERENT AND THESE CORRESPOND TO THE NOISE POINTS IN THE
%DATA. THE TWO REGIONS ARE VERY WELL SEPARATED IN THESE CASES

%CORRELATION COEFF
r_square=zeros(1,4);
r=corrcoef(x1,y1);
r_square(1)=r(1,2).^2;
r=corrcoef(x2,y2);
r_square(2)=r(1,2).^2;
r=corrcoef(x3,y3);
r_square(3)=r(1,2).^2;
r=corrcoef(x4,y4);
r_square(4)=r(1,2).^2;
%IDEALLY CORRELLATION COEFFICIENTS SHOULD BE +1

%REST OF PART 6 involving modification of h(t) vector to reduce
%overfitting
sta_array_new=sta_array;

%%NEED TO REMOVE MIN ELEMENT, NOT SORT THE ORIGINAL ARRAY
r_square_ar=zeros(4,100);    
for j=1:100
    for p=1:4
       ar=abs(sta_array_new(p,:));  
       for d=1:100
           if(ar(d)==0)
               ar(d)=inf;
           end
       end
       [minvalue,pos]=min(ar); 
       sta_array_new(p,pos)=0; %min abs value element in h(t)(STA) is set to 0
    end
    y_est=zeros(4,15099);
    for i=1:4
        y_est(i,:)=conv(sta_array_new(i,:),train_stm);
    end
    lambda_est=zeros(4,15000);
    lambda_est(:,1:15000)=countspikes_train(:,1:15000);
    x=y_est(:,1:15000);
    y=lambda_est(:,1:15000);
    
    bin_size = 100; %estimated y and estimated lambda grouped into bins 
    %for better visibility of scatter plots and better curve fitting
        for i = 1:ceil(15000/bin_size)
            e = i*bin_size;
            if e>15000
                e=15000;
            end
            x1(i)= mean(x(1,(1+(i-1)*bin_size):e));
            x2(i)= mean(x(2,(1+(i-1)*bin_size):e));
            x3(i)= mean(x(3,(1+(i-1)*bin_size):e));
            x4(i)= mean(x(4,(1+(i-1)*bin_size):e));
        
            y1(i)= mean(y(1,(1+(i-1)*bin_size):e));
            y2(i)= mean(y(2,(1+(i-1)*bin_size):e));
            y3(i)= mean(y(3,(1+(i-1)*bin_size):e));
            y4(i)= mean(y(4,(1+(i-1)*bin_size):e));
        end
        
        poly_order=1;
        %%Sigmoid fitting function approximation is given below:
        fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)));
        p1=polyfit(x1,y1,poly_order);
        %p1=sigm_fit(x1,y1,[],[],0);
        
        %p2=sigm_fit(x2,y2,[],[],0);
        p2=polyfit(x2,y2,poly_order);
        
        p3=polyfit(x3,y3,poly_order);
        %p3=sigm_fit(x3,y3,[],[],0);
        
        p4=polyfit(x4,y4,poly_order);
        %p4=sigm_fit(x4,y4,[],[],0);
        
    windowsize=100;
    sta_array_test=zeros(4,100);
    test_stm_part4=zeros(1,5100);
    test_stm_part4(windowsize+1:5100)=test_stm;
    test_stm_part4(1:windowsize)=0;
    numspike=0;
    for i=1:4
        avg=zeros(1,100);
        for f=1:50
            curr_spike_times_train=test_spikes{i,f};
            for k=1:length(curr_spike_times_train)
                window=1+floor(curr_spike_times_train(k)*1000)-(windowsize-1):1+floor(curr_spike_times_train(k)*1000);
                if(window(1)>=1)
                    avg=avg+test_stm_part4(flip(window)); %Flip so that delays from left to right increasing
                end 
            end
            numspike=numspike+length(curr_spike_times_train);     
        end
        avg=avg./numspike;
        sta_array_test(i,:)=avg;
    end

    y_pred=zeros(4,5099);
    for i=1:4
       y_pred(i,:)=conv(sta_array_test(i,:),test_stm);
    end

    %lambda_est=zeros(4,5000);
    fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)));    
    x=countspikes_test(:,1:5000);
    y=y_pred(:,1:5000);
    bin_size = 100; %estimated y and estimated lambda grouped into bins 
    %for better visibility of scatter plots and better curve fitting
        for i = 1:ceil(5000/bin_size)
            e = i*bin_size;
            if e>5000
                e=5000;
            end
            x1(i)= mean(x(1,(1+(i-1)*bin_size):e));
            x2(i)= mean(x(2,(1+(i-1)*bin_size):e));
            x3(i)= mean(x(3,(1+(i-1)*bin_size):e));
            x4(i)= mean(x(4,(1+(i-1)*bin_size):e));

            y1(i)= mean(y(1,(1+(i-1)*bin_size):e));
            y2(i)= mean(y(2,(1+(i-1)*bin_size):e));
            y3(i)= mean(y(3,(1+(i-1)*bin_size):e));
            y4(i)= mean(y(4,(1+(i-1)*bin_size):e));
        end
        %y1=fsigm(p1,y1); %p1,p2,p3,p4 are parameters learnt from training data
        %y2=fsigm(p2,y2);
        %y3=fsigm(p3,y3);
        %y4=fsigm(p4,y4);
        y1=polyval(p1,y1); %p1,p2,p3,p4 are parameters learnt from training data
        y2=polyval(p2,y2);
        y3=polyval(p3,y3);
        y4=polyval(p4,y4);
        
        
        r=corrcoef(x1,y1);
        r_square_ar(1,j)=r(1,2).^2;
        r=corrcoef(x2,y2);
        r_square_ar(2,j)=r(1,2).^2;
        r=corrcoef(x3,y3);
        r_square_ar(3,j)=r(1,2).^2;
        r=corrcoef(x4,y4);
        r_square_ar(4,j)=r(1,2).^2;
        
end
figure
for k=1:4
    subplot(2,2,k);
    plot(r_square_ar(k,:));
    title("r-squared vs number of parameters removed");
    xlabel("Number of parameters removed");
    ylabel("square of correlation coeff.");
end

%In the cases of all 4 neurons, a peak in test accuracy from correlation coeff.
%is obtained for certain number of parameters removed


%% SIGMOID FUNCTION FOR CURVE FITTING
%%THIS FUNCTION IS DOWNLOADED FROM MATHWORKS WEBSITE FOR SIGMOID CURVE
%%FITTING
function [param,stat]=sigm_fit(x,y,fixed_params,initial_params,plot_flag)
% Optimization of parameters of the sigmoid function
%
% Syntax:
%       [param]=sigm_fit(x,y)       
%
%       that is the same that
%       [param]=sigm_fit(x,y,[],[],[])     % no fixed_params, automatic initial_params
%
%       [param]=sigm_fit(x,y,fixed_params)        % automatic initial_params
%       [param]=sigm_fit(x,y,[],initial_params)   % use it when the estimation is poor
%       [param]=sigm_fit(x,y,fixed_params,initial_params,plot_flag)
%
% param = [min, max, x50, slope]
%
% if fixed_params=[NaN, NaN , NaN , NaN]        % or fixed_params=[]
% optimization of "min", "max", "x50" and "slope" (default)
%
% if fixed_params=[0, 1 , NaN , NaN]
% optimization of x50 and slope of a sigmoid of ranging from 0 to 1
%
%
% Additional information in the second output, STAT
% [param,stat]=sigm_fit(x,y,fixed_params,initial_params,plot_flag)
%
%
% Example:
% %% generate data vectors (x and y)
% fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)))
% param=[0 1 5 1];  % "min", "max", "x50", "slope"
% x=0:0.1:10;
% y=fsigm(param,x) + 0.1*randn(size(x));
%
% %% standard parameter estimation
% [estimated_params]=sigm_fit(x,y)
%
% %% parameter estimation with forced 0.5 fixed min
% [estimated_params]=sigm_fit(x,y,[0.5 NaN NaN NaN])
%
% %% parameter estimation without plotting
% [estimated_params]=sigm_fit(x,y,[],[],0)
%
%
% Doubts, bugs: rpavao@gmail.com
% Downloaded from http://www.mathworks.com/matlabcentral/fileexchange/42641-sigmoid-logistic-curve-fit
% warning off
x=x(:);
y=y(:);
if nargin<=1 %fail
    fprintf('');
    help sigm_fit
    return
end
automatic_initial_params=[quantile(y,0.05) quantile(y,0.95) NaN 1];
if sum(y==quantile(y,0.5))==0
    temp=x(y==quantile(y(2:end),0.5));    
else
    temp=x(y==quantile(y,0.5));
end
automatic_initial_params(3)=temp(1);
if nargin==2 %simplest valid input
    fixed_params=NaN(1,4);
    initial_params=automatic_initial_params;
    plot_flag=1;    
end
if nargin==3
    initial_params=automatic_initial_params;
    plot_flag=1;    
end
if nargin==4
    plot_flag=1;    
end
if exist('fixed_params','var')
    if isempty(fixed_params)
        fixed_params=NaN(1,4);
    end
end
if exist('initial_params','var')
    if isempty(initial_params)
        initial_params=automatic_initial_params;
    end
end
if exist('plot_flag','var')
    if isempty(plot_flag)
        plot_flag=1;
    end
end
%p(1)=min; p(2)=max-min; p(3)=x50; p(4)=slope como em Y=Bottom + (Top-Bottom)/(1+10^((LogEC50-X)*HillSlope))
%f = @(p,x) p(1) + (p(2)-p(1)) ./ (1 + 10.^((p(3)-x)*p(4)));
f_str='f = @(param,xval)';
free_param_count=0;
bool_vec=NaN(1,4);
for i=1:4;
    if isnan(fixed_params(i))
        free_param_count=free_param_count+1;
        f_str=[f_str ' param(' num2str(free_param_count) ')'];
        bool_vec(i)=1;
    else
        f_str=[f_str ' ' num2str(fixed_params(i))];
        bool_vec(i)=0;
    end
    if i==1; f_str=[f_str ' + (']; end
    if i==2;
        if isnan(fixed_params(1))            
            f_str=[f_str '-param(1) )./ (   1 + 10.^( (']; 
        else
            f_str=[f_str '-' num2str(fixed_params(1)) ')./ (1 + 10.^((']; 
        end
    end    
    if i==3; f_str=[f_str ' - xval ) *']; end
    if i==4; f_str=[f_str ' )   );']; end
end
eval(f_str)
[BETA,RESID,J,COVB,MSE] = nlinfit(x,y,f,initial_params(bool_vec==1));
stat.param=BETA';
% confidence interval of the parameters
stat.paramCI = nlparci(BETA,RESID,'Jacobian',J);
% confidence interval of the estimation
[stat.ypred,delta] = nlpredci(f,x,BETA,RESID,'Covar',COVB);
stat.ypredlowerCI = stat.ypred - delta;
stat.ypredupperCI = stat.ypred + delta;
% plot(x,y,'ko') % observed data
% hold on
% plot(x,ypred,'k','LineWidth',2)
% plot(x,[lower,upper],'r--','LineWidth',1.5)
free_param_count=0;
for i=1:4;
    if isnan(fixed_params(i))
        free_param_count=free_param_count+1;
        param(i)=BETA(free_param_count);
    else
        param(i)=fixed_params(i);
    end    
end
    
if plot_flag==1 
    x_vector=min(x):(max(x)-min(x))/100:max(x);
    plot(x,y,'k.',x_vector,f(param(isnan(fixed_params)),x_vector),'r-')
    xlim([min(x) max(x)])
end


%% PART 7
figure()
for neuron_cnt=1:4
    num_trial=10;
    MI=zeros(num_trial,7);
    for vp_trial=1:num_trial
        q=[0, 0.001, 0.01, 0.1, 1, 10, 100];
        rng(vp_trial);
        start_pt=rand(1,8);
        start_pt=start_pt*15;
        end_pt=start_pt+0.1;


        for trial_no = 1:50
            for stpt=1:8
                temp1=[];
                temp1=All_Spike_Times{neuron_cnt,trial_no}>start_pt(stpt);
                spktemp1=All_Spike_Times{neuron_cnt,trial_no}(temp1);
                temp2=[];
                temp2=All_Spike_Times{neuron_cnt,trial_no}<=start_pt(stpt)+0.1;
                spktemp2=All_Spike_Times{neuron_cnt,trial_no}(temp2);
                vpmat{stpt,trial_no}=intersect(spktemp1,spktemp2);
            end
        end

        confusion=zeros(7,8,8);

        for qno=1:7
            vpmindist=zeros(8,50);
            vpmindist=inf+vpmindist;
            for m=1:50
                for n=1:8
                    minj=[];
                    for i = 1:50
                        for j = 1:8
                            if ~(m==i && n==j)
                                if(spkd(vpmat{j,i},vpmat{n,m},q(qno))<=vpmindist(n,m))
                                    vpmindist(n,m)=spkd(vpmat{j,i},vpmat{n,m},q(qno));
                                    if(spkd(vpmat{j,i},vpmat{n,m},q(qno))==vpmindist(n,m))
                                        minj=[minj,j];
                                    else
                                        minj=[j];
                                    end

                                end
            %                     vpdist(n,m)=min(a,vpdist(m,n));
                            end
                        end
                    end
                    for it=1:length(minj)
                        confusion(qno,n,minj(it))=confusion(qno,n,minj(it))+1/length(minj);
                    end
                end
            end



        end
        confusion=confusion/400;
        for qno=1:7
            for x = 1:8
                    for y=  1:8
                        MI(vp_trial,qno)=MI(vp_trial,qno)+confusion(qno,x,y)*log2(confusion(qno,x,y)/(sum(confusion(qno,x,:))*sum(confusion(qno,:,y))));
                    end
            end
        end
    end
    MIavg=sum(MI(:,:));
    MIavg=MIavg/num_trial;
    subplot(2,2,neuron_cnt)
    plot(linspace(-3,2,6),MIavg(2:7));
    title(strcat('neuron ',int2str(neuron_cnt),' MI vs q'))
    xlabel('q')
    ylabel('MI')

end

        
function d=spkd(tli,tlj,cost)
%
% d=spkd(tli,tlj,cost) calculates the "spike time" distance
% (Victor & Purpura 1996) for a single cost
%
% tli: vector of spike times for first spike train
% tlj: vector of spike times for second spike train
% cost: cost per unit time to move a spike
%
%  Copyright (c) 1999 by Daniel Reich and Jonathan Victor.
%  Translated to Matlab by Daniel Reich from FORTRAN code by Jonathan Victor.
%
nspi=length(tli);
nspj=length(tlj);

if cost==0
   d=abs(nspi-nspj);
   return
elseif cost==Inf
   d=nspi+nspj;
   return
end

scr=zeros(nspi+1,nspj+1);
%
%     INITIALIZE MARGINS WITH COST OF ADDING A SPIKE
%
scr(:,1)=(0:nspi)';
scr(1,:)=(0:nspj);
if nspi & nspj
   for i=2:nspi+1
      for j=2:nspj+1
         scr(i,j)=min([scr(i-1,j)+1 scr(i,j-1)+1 scr(i-1,j-1)+cost*abs(tli(i-1)-tlj(j-1))]);
      end
   end
end
d=scr(nspi+1,nspj+1);
end