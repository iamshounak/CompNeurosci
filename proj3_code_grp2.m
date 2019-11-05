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
    
    figure
    p1=polyfit(x1,y1,poly_order);
    x_axis = linspace(min(x1),max(x1));
    y_axis = polyval(p1,x_axis);
    subplot(2,2,1)
    plot(x1,y1,'o')
    hold on
    plot(x_axis,y_axis)
    hold off
    title('f(t) estimated from scatter plot between y(t) and \lambda(t)')
    xlabel('y(t) = convolution of s(t) and h(t)')
    ylabel('\lambda(t)')
    
    p2=polyfit(x2,y2,poly_order);
    x_axis = linspace(min(x2),max(x2));
    y_axis = polyval(p2,x_axis);
    subplot(2,2,2)
    plot(x2,y2,'o')
    hold on
    plot(x_axis,y_axis)
    hold off
    title('f(t) estimated from scatter plot between y(t) and \lambda(t)')
    xlabel('y(t) = convolution of s(t) and h(t)')
    ylabel('\lambda(t)')
    
    p3=polyfit(x3,y3,poly_order);
    x_axis = linspace(min(x3),max(x3));
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
    x_axis = linspace(min(x4),max(x4));
    y_axis = polyval(p4,x_axis);
    subplot(2,2,4)
    plot(x4,y4,'o')
    hold on
    plot(x_axis,y_axis)
    hold off
    title('f(t) estimated from scatter plot between y(t) and \lambda(t)')
    xlabel('y(t) = convolution of s(t) and h(t)')
    ylabel('\lambda(t)')