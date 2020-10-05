% one dimensional user distribution scenario
clear all;
close all;
clc;


d = 1; % User space dimension
no_samples = 2000;
r = 2; % Pathloss component
P =10^(75/10); % Transmitter power in Watts
delta = 0.9; % Attenuation factor
h = 1000; % UAV Altitude
radius = 1000; % Area radius
% Rate parameters
b = 0.43;
b1 = (1-delta)*b;
c = 4.88;
c1 = c*exp(-b*(pi/2-c));
c2 = .5*(b^2*(1-delta)+r*(c1-1));
n = 1:20;
x = radius*rand(d,no_samples); % user samples uniformly distributed in area of interest.
deploy_opt_PSO = zeros(length(n),length(n)); % Initializing the optimal deployment of UAVs Derived with PSO method
deploy_opt_asymp = zeros(length(n),length(n)); % Initializing the optimal deployment of UAVs Derived with our proposed asymptotic approach
K=1; % Number of iterations


% Solving with PSO method
for k=1:K
    for j=1:length(h)
        for i=1:length(n)
            [cost_opt, deploy_opt_PSO(i,1:n(i),j)] = PSO(r,x,h(j),n(i),P,delta)
            R_opt_PSO(i,j,k) = cost_opt
        end
    end
end

% assymptotic case
for j=1:length(h)
    for i=1:length(n)
        [idx,u,sumd,D] = kmeans(x,n(i));
        distance = sqrt(min(D,[],2));
        deploy_opt_asymp(i,1:n(i),j) = u'
	% Asymptotic 
        cost_opt_asymp(i,j) = (1-delta)*(b*c1/(h(j)))*costcalculator(deploy_opt_asymp(i,1:n(i),j),n(i),x,idx); 
        COST = costcalculator(deploy_opt_asymp(i,1:n(i),j),n(i),x,idx)
        R_opt_largeH(i,j) =(P./(log(2)*(1+c1)^2*(h(j)^r)))*((1+delta*c1)*(1+c1) +  cost_opt_asymp(i,j))
        R_opt_largeN(i,j) =(log(1+P/h(j)^r)*(1/(c1+1))+log(1+P*delta/h(j)^r)*(c1/(c1+1))+( log((h(j)^r+P*delta)/(h(j)^r+P))*(b*c1/(h(j)*(1+c1)^2))*costcalculator(deploy_opt_asymp(i,1:n(i),j),n(i),x,idx)))/log(2);
        R_opt_asymp(i,j) = -ObjectiveFunction ( deploy_opt_asymp(i,1:n(i),j),n(i),r,x,h(j),delta,P,b,c );
        R_opt_quant(i,j) = (log(1+P/h(j)^r)*(1/(c1+1))+log(1+P*delta/h(j)^r)*(c1/(c1+1))+( log((h(j)^r+P*delta)/(h(j)^r+P))*(b*c1/(h(j)*(1+c1)^2))*250/n(i)))/log(2);
        Plos = 1./(1+c*exp(-b*(atan(h(j)./distance)-c)));
        R_los = log2(1+(P./(distance.^2+h(j).^2).^(r/2)));
        R_nlos = log2(1+(P*delta./(distance.^2+h(j).^2).^(r/2)));
        R_new(i,j) = sum(R_los.*Plos+R_nlos.*(1-Plos))/no_samples;
    end
end


plot(R_opt_largeH,'red')
hold on;
plot(R_opt_largeN,'black')
plot(R_opt_PSO)
plot(R_new,'green')

Logreal = log2(1+(P./(distance.^2+h(j).^2).^(r/2)));
Logasymp = (P/h(j)^r).*(1-r/2*distance.^2/h(j)^2)/log(2);
%%
% two dim twitter 
clear all;
close all;
clc;

d = 2;
r = 2;
no_samples = 1000;
P =10^(50/10);
delta = 0.9;
h =[100,300,500];
radius = 10000;
sigma = 100;
b = 0.43;
b1 = (1-delta)*b;
c = 4.88;
c1 = c*exp(-b*(pi/2-c));
c2 = .5*(b^2*(1-delta)+r*(c1-1));
n = 1:20;
%x = sigma*randn(d,no_samples);
%x = radius*rand(d,no_samples)
data = load('X_train.mat');
x_total = data.x;
no_samples = size(x_total,1);
x_radian = x_total(:,:,1)';
%https://www.baeldung.com/java-convert-latitude-longitude, Lattitude phi,
%longitude lambda
earthRadius = 6378137;
lat = x_radian(1,:)*pi/180;
long = x_radian(2,:)*pi/180;
x(1,:) = earthRadius*long;
x(2,:) = earthRadius*log(tan(pi/4+lat/2));
deploy_opt_PSO = zeros(2,length(n),length(n),length(h));
deploy_opt_asymp = zeros(2,length(n),length(n),length(h));
kd = 0.3772;

for i=1:no_samples
    L(i) = norm(x(:,i));
end

falpha = mean(sum((exp(-L./(2*pi*(sigma)))).^(2/3)))^(3/2);

% Solving with PSO method
for j=1:length(h)
    for i=1:length(n)
        [cost_opt, deploy_opt_PSO(:,1:n(i),i,j)] = PSOxy(r,x,h(j),n(i),P,delta,radius,sigma)
        R_opt_PSO(i,j) = cost_opt
    end
end

% assymptotic case
for j=1:length(h)
    for i=1:length(n)
        [idx,u,sumd,D] = kmeans(x',n(i));
        %[distance,u] = MyKmeans(x,n(i))
        distance = sqrt(min(D,[],2));
        %distance = sumd;
        deploy_opt_asymp(:,1:n(i),i,j) = u'
        %cost_opt_asymp(i,j) = (1-delta)*(b*c1/(h(j)))*costcalculatorXY(deploy_opt_asymp(i,1:n(i),j),n(i),x,idx); %- sum(sumd)*(1+delta*c1)*(1+c1)/h(j)^2%/no_samples;
        R_opt_largeN(i,j) =(log(1+P/h(j)^r)*(1/(c1+1))+log(1+P*delta/h(j)^r)*(c1/(c1+1))+( log((h(j)^r+P*delta)/(h(j)^r+P))*(b*c1/(h(j)*(1+c1)^2))*costcalculatorXY(deploy_opt_asymp(:,1:n(i),i,j),n(i),x,idx)))/log(2);
        R_opt_quant(i,j) =(log(1+P/h(j)^r)*(1/(c1+1))+log(1+P*delta/h(j)^r)*(c1/(c1+1))+( log((h(j)^r+P*delta)/(h(j)^r+P))*(b*c1/(h(j)*(1+c1)^2))*kd*falpha*(n(i)^(-0.5))))/log(2);
        %R_opt_asymp(i,j) = -ObjectiveFunction ( deploy_opt_asymp(:,1:n(i),i,j),n(i),r,x,h(j),delta,P,b,c );
    end
end

for j=1:length(h)
    hold on
    plot(n,R_opt_PSO(:,j),'blue')
    plot(n,R_opt_largeN(:,j),'black')
    plot(n,R_opt_quant(:,j),'green')
end
hold off;

legend('Iterative Approach, h = 100m','Quantization Theory Approach, h = 100m','PSO, h = 100m','Iterative Approach, h = 300m','Quantization Theory Approach, h = 300m','PSO, h = 300m','Iterative Approach, h = 500m','Quantization Theory Approach, h = 500m','PSO, h = 500m' )
%%
% one dim TIME variate - trajectory optimization
clear all;
close all;
clc;

d = 1;
no_samples = 5000;
r = 2;
P = 10^(50/10);
delta = 0.5;
h = 100%[300,100,50];
radius = 1000;
b = 0.43;
b1 = (1-delta)*b;
c = 4.88;
c1 = c*exp(b*c);
c2 = .5*(b^2*(1-delta)+r*(c1-1));
n = 5;
T = -1:.1:1;
deploy_opt_PSO = zeros(n,length(T));
deploy_opt_asymp = zeros(n,length(T));
deploy_opt_quant = zeros(n,length(T));
R_opt_PSO = zeros(length(T));
K = 100;

for t = 1:length(T)
    q = rand(1,no_samples);
    x = 2-2*abs(T(t))+q.^(1/(1+2*abs(T(t))));
    %PSO result
    for k=1:K
        [cost_opt, temp] = PSO(r,x,h,n,P,delta);
        deploy_opt_PSO_temp(:,t,k) = sort(temp)
    end
    %asymp result 
    [idx,u] = kmeans(x,n);
     deploy_opt_asymp(:,t) = sort(u);
    %Quantization result
    deploy_opt_quant(:,t) =  2-2*abs(T(t))+((2*(1:n)-1)/(2*n)).^(1/(1+abs(T(t))));
end

R_opt_PSO(t) = cost_opt;
deploy_opt_PSO = mean(deploy_opt_PSO_temp,3);
%deploy_opt_PSO = min(deploy_opt_PSO_temp,[],3);
for t=1:length(T)
    if mod(t,2) == 1
        temp1(:,(t+1)/2) = deploy_opt_PSO(:,t);
        temp2(:,(t+1)/2) = deploy_opt_asymp(:,t);
        temp3(:,(t+1)/2) = deploy_opt_quant(:,t);
        T_new((t+1)/2) = T(t); 
    end
end

figure(2)
hold on
for i = 1:n
    plot(T,deploy_opt_PSO(i,:),'black')
    plot(T,deploy_opt_asymp(i,:),'blue')   
    plot(T,deploy_opt_quant(i,:),'red')    
end

figure(3)
hold on
for i = 1:n
    plot(T_new,temp1(i,:),'black','LineWidth',2)
    plot(T_new,temp2(i,:),'blue','LineWidth',2)   
    plot(T_new,temp3(i,:),'red','LineWidth',2)    
end
xlabel('Trajectories of 5 UAVs in a one-dimensional network')
ylabel('UAV location')
legend('PSO','Thorem 2','Quantization theory')
%%

% two dim trentino
clear all;
close all;
clc;

d = 2;
r = 2;
no_samples = 10000;
P =10^(50/10);
delta = 0.9;
h =10%[50,100,300,500];
radius = 10000;
sigma = 100;
b = 0.43;
b1 = (1-delta)*b;
c = 4.88;
c1 = c*exp(-b*(pi/2-c));
c2 = .5*(b^2*(1-delta)+r*(c1-1));
n = 1:1:10;
grid_side = 235;
% %x = sigma*randn(d,no_samples);
% %x = radius*rand(d,no_samples)
% %lat_long_data = load('trentino_centroid_lat_long.mat');
% %lat_long_data = lat_long_data.x;
% Internet_data = load('dataset_internet_t0_final.mat');
% Internet_data = Internet_data.x;
% Internet_data_scaled = Internet_data;
% Internet_data_scaled(:,4)= abs(Internet_data_scaled(:,4));
% Internet_data_scaled(:,5) = (Internet_data_scaled(:,4))/(sum(Internet_data_scaled(:,4)));
% Internet_data_scaled(:,5)= abs(Internet_data_scaled(:,5));
% earthRadius = 6378137;
% 
% lat = Internet_data_scaled(:,2)*pi/180;
% long = Internet_data_scaled(:,3) *pi/180;
% Internet_data_scaled(:,2) = earthRadius*long; 
% Internet_data_scaled(:,3) = earthRadius*log(tan(pi/4+lat/2));
% %Internet_data_scaled(:,2) = (Internet_data_scaled(:,2)-min(Internet_data_scaled(:,2)))/(max(Internet_data_scaled(:,2))-min(Internet_data_scaled(:,2)));
% %fq = zeros(size(Internet_data_scaled,1),1,2);
% 
% %Creating CDF
% Internet_data_scaled(:,6) = zeros(length(Internet_data_scaled(:,1)),1); 
% temp(1) = Internet_data_scaled(1,5);
% for i=2:length(Internet_data_scaled(:,1))
%     temp(i) = temp(i-1)+Internet_data_scaled(i,5)
% end
% Internet_data_scaled(:,6) = temp'

data = load('Internet_data_scaled_finalNew.mat');
Internet_data_scaled = data.Internet_data_scaled;
Internet_data_scaled(:,7) = zeros(length(Internet_data_scaled(:,1)),1);
temp = Internet_data_scaled(:,6);
rnd_num = rand(no_samples,1);
for i = 1:no_samples
    temp_index = max(find(temp<=rnd_num(i)));
    Internet_data_scaled( temp_index ,7) = Internet_data_scaled( temp_index ,7)+1;
    
end


x=[];
for i=1:length(Internet_data_scaled(:,1))
    num_users = Internet_data_scaled(i,7);
    x_sub = grid_side/2*(2*rand(2,num_users)-1)+repmat([Internet_data_scaled(i,2);Internet_data_scaled(i,3)],1,num_users);
    x = cat(2,x,x_sub);
end

x_new=x;
x_new(1,:)=(x(1,:)-min(x(1,:)))/(max(x(1,:))-min(x(1,:)));
x_new(2,:)=(x(2,:)-min(x(2,:)))/(max(x(2,:))-min(x(2,:)));

% Solving with PSO method
for j=1:length(h)
    for i=1:length(n)
        [cost_opt, deploy_opt_PSO(:,1:n(i),i,j)] = PSOxy(r,x_new,h(j),n(i),P,delta,radius,sigma)
        R_opt_PSO(i,j) = cost_opt
        
    end
end
deploy_opt_PSO_final=deploy_opt_PSO;
deploy_opt_PSO_final(1,:,:) = deploy_opt_PSO(1,:,:)*(max(x(1,:))-min(x(1,:)))+min(x(1,:));
deploy_opt_PSO_final(2,:,:) = deploy_opt_PSO(2,:,:)*(max(x(2,:))-min(x(2,:)))+min(x(2,:));
for i=1:length(n)
cost_PSO_opt(i) = ObjectiveFunctionXY ( deploy_opt_PSO_final(:,:,i),n(i),r,x,h,delta,P,b,c );
end
fq = Internet_data_scaled(:,5);
falpha = sum((fq.^(2/3))/length(Internet_data_scaled(:,5)))^(3/2);
kd = 0.3772;
for j=1:length(h)
    for i=1:length(n)
        [idx,u,sumd,D] = kmeans(x',n(i));
        %[distance,u] = MyKmeans(x,n(i))
        distance = sqrt(min(D,[],2));
        %distance = sumd;
        deploy_opt_asymp(:,1:n(i),i,j) = u'
        %cost_opt_asymp(i,j) = (1-delta)*(b*c1/(h(j)))*costcalculatorXY(deploy_opt_asymp(i,1:n(i),j),n(i),x,idx); %- sum(sumd)*(1+delta*c1)*(1+c1)/h(j)^2%/no_samples;
        R_opt_largeN(i,j) =(log(1+P/h(j)^r)*(1/(c1+1))+log(1+P*delta/h(j)^r)*(c1/(c1+1))+( log((h(j)^r+P*delta)/(h(j)^r+P))*(b*c1/(h(j)*(1+c1)^2))*costcalculatorXY(deploy_opt_asymp(:,1:n(i),i,j),n(i),x,idx)))/log(2);
        R_opt_quant(i,j) =(log(1+P/h(j)^r)*(1/(c1+1))+log(1+P*delta/h(j)^r)*(c1/(c1+1))+( log((h(j)^r+P*delta)/(h(j)^r+P))*(b*c1/(h(j)*(1+c1)^2))*kd*falpha*(n(i)^(-0.5))))/log(2);
        %R_opt_asymp(i,j) = -ObjectiveFunction ( deploy_opt_asymp(:,1:n(i),i,j),n(i),r,x,h(j),delta,P,b,c );
    end
end

for j=1:length(h)
    hold on
    
    plot(n,R_opt_quant(:,j),'green')
    plot(n,R_opt_largeN(:,j),'black')
    plot(n,R_opt_PSO(:,j),'blue')
end
hold off;

legend('Quantization Theory Approach, h = 100m','Iterative Approach, h = 100m','PSO, h = 100m')
xlabel('n')
ylabel('cost')