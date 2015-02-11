%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: A function triggered when the model Simultaneous_Wrist_and_Hand.mdl  is stopped.
% If the synergy matrix does not existed (which means the model has not been trained), 
% this function will use NMF alogrithm to calculate the synergy matrix.
%--------------------------------------------------------------------------------------------------------------------------------
% Used in: Simultaneous_Wrist_and_Hand.mdl
%--------------------------------------------------------------------------------------------------------------------------------
% Update : 2015/2/12 by Jiaxin MA @Kyoto Univ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('Feature_Matrix.mat');
load('Provided_Cues_for_NMF.mat');
if NMF_threshold < 0
    
    % Use the recorded feature matrix to calculate a new H matrix.
    
    Feature_Matrix = Feature_Matrix(2:end,2:end)';
    Provided_Cues = Provided_Cues(2:end,2:end)';
    
    %%
    [W H] = Custom_NMF(Feature_Matrix, Provided_Cues, 0, 0, 9, 0, 'single',0);
    % headcut, tailcut, trail_length, feature_nomalize,
    % train_method(single/mult),remove_baseline
    
    % Save H.mat in this folder.
    save('H.mat','H');
    
%     base = mean(W(41:240,:));  %take 1~6s data as baseline
%     Th = 0.1*max(W(41:end,:))+base;
%     save('para.mat','Th');
    
    % Draw the smoothed output
    ControlScheme;
    %    state=log10(state*9+1); % take log conversion if necessary
    
%     Fpass=8;
%     fs=600;
%     wn=Fpass/(fs/2);
%     N=ceil(6.6/wn);
%     b=fir1(N,wn);
%     for i=1:4
%     state(:,i)=filtfilt(b,1,W(:,i));
%     end
%     state(state < 0)=0;
%     state(state > 1)=1;
    
    for i=1:4
        subplot(5,4,1+i*4:3+i*4);
        hold on
        plot(state(:,i),'r');
    end
    %%
    % Save this session.
    time = datestr(now);
    time = time(time~='-');
    time(time==':')='_';
    time(time==' ')='_';
    cd('./data');
    save(strcat('NMFtest_',time,'.mat'),'Feature_Matrix','Provided_Cues','H','W','state');
    cd('..');
    disp('NMF matrix genarated successful.');
else
    %%
    load('W.mat');
    load('state.mat');
    W=W(2:end,2:end)';
    state=state(2:end,2:end)';
    
%     Th=[0.1 0.1 0.1 0.1];
%     % if want to modify the control scheme, uncomment this to test
%     ControlScheme;
    
    figure;
    for i=1:4
        subplot(4,1,i);
        plot(W(:,i));
        hold on
        plot(state(:,i),'r');
        ylim([-0.5,1.5]);
    end
    
end
