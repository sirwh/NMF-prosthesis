%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: An s-function to present the visual instruction of training process
%--------------------------------------------------------------------------------------------------------------------------------
% Used in: Simultaneous_Wrist_and_Hand.mdl
%--------------------------------------------------------------------------------------------------------------------------------
% Update : 2015/2/12 by Jiaxin MA @Kyoto Univ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sys,x0,str,ts] = Cue_Presentation(t,x,u,flag,sequence)

if flag==0
    
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 0;
    sizes.NumOutputs     = 0;
    sizes.NumInputs      = -1; % dynamically sized
    sizes.DirFeedthrough = 0;  % has no direct feedthrough
    sizes.NumSampleTimes = 1;
    
    sys = simsizes(sizes);
    str = [];
    x0  = [];
    ts  = [1 0];   % sample time
    
    
elseif flag==2
    if u==0,return,end
    DelayBegin = 0;
    if t < DelayBegin,return,end
    if round(t)==t;
        time = formatTimeFcn(t);
        session_timer = findobj('Tag','timerStatus');
        set(session_timer,'String',time);
    end
    
elseif flag==9
    if u==0,return,end
    close Cue_Presentation_GUI;
    
end

    function output = formatTimeFcn(float_time)
        cd('Hand Images');
        % Format the Time String
        float_time = floor(float_time);
        Cue_Presentation_GUI_handle = findobj('type', 'figure', 'Name', 'Cue_Presentation_GUI');
        % Cue_Presentation_GUI_handle.text1=text(0,0,'ready','Fontsize',30,'HorizontalAlignment','center','Visible','Off','EraseMode','Background','Color','r');
        set(0,'CurrentFigure',Cue_Presentation_GUI_handle);
        hint='';
        total_time = length(sequence);
        tcn=float_time-DelayBegin+1;
        if tcn >= 1 && tcn <= total_time;
            switch(sequence(tcn,3))
                case 1
                    disp(strcat('hand_',num2str(sequence(tcn,1)),'.jpg'));
                    hint='ready';
                    imshow(strcat('hand_',num2str(sequence(tcn,1)),'.jpg'));
                case 2
                    disp(strcat('hand_',num2str(sequence(tcn,1)),'.jpg'));
                    hint='begin';
%                     imshow(strcat('hand_',num2str(sequence(tcn,1)),'.jpg'));
                case 3
                    disp(strcat('red_hand_',num2str(sequence(tcn,1)),'.jpg'));
                    hint='hold';
                    imshow(strcat('red_hand_',num2str(sequence(tcn,1)),'.jpg'));
                case 4
                    disp(strcat('hand_',num2str(sequence(tcn,1)),'.jpg'));
                    hint='return';
                    imshow(strcat('hand_',num2str(sequence(tcn,1)),'.jpg'));
                otherwise
                    error('check Cue_Presentation.m');                    
            end
        end
        cd('..');
        
%         mins = floor(float_time/60);
%         secs = float_time - 60*(mins);
%         time_left = total_time - (mins*60+secs) + DelayBegin;
%         
%         mins = floor(time_left/60);
%         secs = mod(time_left,60);
%         m = sprintf('%1.0f:',mins);
%         s = sprintf('%1.0f',secs);
%         
%         if mins < 10
%             m = sprintf('0%1.0f:',mins);
%         end
%         if secs < 9.9995
%             s = sprintf('0%1.0f',secs);
%         end
%         output = [m s];
        output = sprintf('%s',hint);
        
        %         if mins == 0 && secs == 0
        %             STOPPED = 1;
        %         end
    end
end


