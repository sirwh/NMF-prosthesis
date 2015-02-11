%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: The function that uses NMF algorithm to decompose EMG matrix into
% synergy matrix W and H, where H for weight, W for activation level.
%--------------------------------------------------------------------------------------------------------------------------------
% Used in: NMF_Stop.m
%--------------------------------------------------------------------------------------------------------------------------------
% Update : 2015/2/12 by Jiaxin MA @Kyoto Univ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [W H] = Custom_NMF(Feature_Matrix, Provided_Cues, headcut,tailcut,trail_len,regularize_data,training_method,remove_baseline)
% Feature_Matrix is t*24, Provided_Cues is t*1

%% trim
step=15;
spr=600;
len = round(trail_len*spr/step);

head=round(headcut*spr/step);  % 1000 is the sample rate
tail=round(tailcut*spr/step);
rows_to_keep = zeros(1,length(Feature_Matrix));

for i=1:len:length(Feature_Matrix)
    rows_to_keep(i+head:i+len-1-tail)=1;
end
% only trim the data for training

%% remove baseline 
BLremover = zeros(1,size(Feature_Matrix,2));
if remove_baseline==1
    Rest_Feature_Matrix = Feature_Matrix(Provided_Cues==1,:);
    BLremover = mean(Rest_Feature_Matrix); % or mean?
    for i=1:size(Rest_Feature_Matrix,2)
        Feature_Matrix(:,i) = Feature_Matrix(:,i) - BLremover(i);
    end
    Feature_Matrix(Feature_Matrix<0)=0;
end
save('BLremover.mat','BLremover');

%% feature normalization
maxf=ones(size(Feature_Matrix,2)/8,1);
if regularize_data == 1
    j=1;
    for i=1:8:size(Feature_Matrix,2)
        maxf(j) = max(max(Feature_Matrix(:,i:i+7)));
        Feature_Matrix(:,i:i+7)=Feature_Matrix(:,i:i+7)./maxf(j);
        j=j+1;
    end
end
save('maxf.mat','maxf');

%% Sort the label, if necessary

% % Class 1 = Rest
% % Class 2 = Open
% % Class 3 = Close
% % Class 4 = Open and Pronate
% % Class 5 = Close and Pronate
% % Class 6 = Open and Supinate
% % Class 7 = Close and Supinate
% % Class 8 = Pronate
% % Class 9 = Supinate
%
% movement_list=[1:9;1 2 3 4 5 6 7 12 13]';
% for i=1:length(Provided_Cues)
%     Provided_Cues(i)=movement_list(Provided_Cues(i)==movement_list(:,2),1);
% end
%
figure;

%% calculate NMF for Open, Close, Pronate, and Supinate
trimmed_Feature_Matrix = Feature_Matrix(rows_to_keep==1,:);
trimmed_Provided_Cues = Provided_Cues(rows_to_keep==1,:);

Open_Feature_Matrix = trimmed_Feature_Matrix(trimmed_Provided_Cues==2,:);
Close_Feature_Matrix = trimmed_Feature_Matrix(trimmed_Provided_Cues==3,:);
Pronate_Feature_Matrix = trimmed_Feature_Matrix(trimmed_Provided_Cues==4,:);
Supinate_Feature_Matrix = trimmed_Feature_Matrix(trimmed_Provided_Cues==5,:);

if strcmp(training_method,'single')
    [W_O H_O] = nnmf(Open_Feature_Matrix,1,'w0',ones(length(Open_Feature_Matrix),1),'alg','mult');
    [W_C H_C] = nnmf(Close_Feature_Matrix,1,'w0',ones(length(Close_Feature_Matrix),1),'alg','mult');
    [W_P H_P] = nnmf(Pronate_Feature_Matrix,1,'w0',ones(length(Pronate_Feature_Matrix),1),'alg','mult');
    [W_S H_S] = nnmf(Supinate_Feature_Matrix,1,'w0',ones(length(Supinate_Feature_Matrix),1),'alg','mult');
    H = [H_O*mean(W_O);H_C*mean(W_C);H_P*mean(W_P);H_S*mean(W_S)];
elseif strcmp(training_method,'mult')
    [W_O H_O] = custom_nnmf(Open_Feature_Matrix,1,'alg','mult');
    [W_C H_C] = custom_nnmf(Close_Feature_Matrix,1,'alg','mult');
    [W_P H_P] = custom_nnmf(Pronate_Feature_Matrix,1,'alg','mult');
    [W_S H_S] = custom_nnmf(Supinate_Feature_Matrix,1,'alg','mult');
    H = [H_O*0.5*max(W_O);H_C*0.5*max(W_C);H_P*0.5*max(W_P);H_S*0.5*max(W_S)];
end

subplot(5,4,8);
barweb(mean(W_O),std(W_O),.5);
set(gca, 'XTickLabel', 'Open');
ylabel({'Synergy 1 ' ;'Activation'});
set(gcf,'color','w');
subplot(5,4,12);
barweb(mean(W_C),std(W_C),.5);
set(gca, 'XTickLabel', 'Close');
ylabel({'Synergy 2 ' ;'Activation'});
set(gcf,'color','w');
subplot(5,4,16);
barweb(mean(W_P),std(W_P),.5);
set(gca, 'XTickLabel', 'Pronate');
ylabel({'Synergy 3 ' ;'Activation'});
set(gcf,'color','w');
subplot(5,4,20);
barweb(mean(W_S),std(W_S),.5);
set(gca, 'XTickLabel', 'Supinate');
ylabel({'Synergy 4 ' ;'Activation'});
set(gcf,'color','w');

% H = [H_O;H_C;H_P;H_S];
W = Feature_Matrix/H;
%% Plot the synergy decomposition using the individually calculated synergies
%from the four primary movement components.
subplot(5,4,1:4);

% Now sort the various movements and plot them based on their synergy
% involvements.
W_mean = zeros(5,4);
W_std = zeros(5,4);
for i = 1:5
    W_mean(i,:) = mean(W(Provided_Cues==i,:));
    W_std(i,:) = std(W(Provided_Cues==i,:));
end
% barweb(W_mean,W_std,.5);
bar(W_mean,'grouped');
title('Non-negative Matrix Factorization for Synergy Detection in Hand and Wrist Movements: Individually Calculated Synergies Held Constant');
movement_names = char('Rest','Open','Close','Pronate','Supinate');
set(gca, 'XTickLabel', movement_names);
ylabel({'Synergy ' ;'Activation'});
set(gcf,'color','w');
ylim([-1 2]);
legend('Synergy 1','Synergy 2','Synergy 3','Synergy 4','location','northeast');

%% calulate base and Th

base = mean(W(41:240,:));  %take 1~6s data as baseline
Th = 0.1*max(W(41:end,:))+base;

%baseline process
for i=1:4
    W(:,i)=(W(:,i)-base(i))*(1/(1-base(i)));
end
Th = [0.1 0.1 0.1 0.1]; % test

save('para.mat','Th','base');

%% plot the time-variance synergy activations
% W=[W(Provided_Cues==1,:);W(Provided_Cues==2,:);W(Provided_Cues==3,:);W(Provided_Cues==8,:); W(Provided_Cues==9,:)];
% Provided_Cues = [Provided_Cues(Provided_Cues==1,:);Provided_Cues(Provided_Cues==2,:);Provided_Cues(Provided_Cues==3,:);Provided_Cues(Provided_Cues==8,:); Provided_Cues(Provided_Cues==9,:)];

j=1;
for i=1:length(Provided_Cues)-1
    if Provided_Cues(i)~=Provided_Cues(i+1)
        change(j)=i+1;
        j=j+1;
    end
end

for i=1:4
    subplot(5,4,1+i*4:3+i*4);
    plot(W(:,i));
    hold on
    % plot(diff(W(:,i)),'g');
    for j=1:length(change)
        plot([change(j),change(j)],[-2,2],'g');
    end
    ylim([-2 2]);
    set(gca,'xtick',[change]);
    set(gca,'xticklabel',int2str(Provided_Cues([change])));%style 4
end