%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: An s-function of an original control scheme which is used to smooth the output signal.
%--------------------------------------------------------------------------------------------------------------------------------
% Used in: Simultaneous_Wrist_and_Hand.mdl
%--------------------------------------------------------------------------------------------------------------------------------
% Update : 2015/2/12 by Jiaxin MA @Kyoto Univ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [sys,x0,str,ts] = sfunction_magic(t,x,u,flag,Th,NMF_threshold)

if flag==0
    
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 4;
    sizes.NumOutputs     = 4;
    sizes.NumInputs      = -1;  % dynamically sized
    sizes.DirFeedthrough = 1;   % has direct feedthrough
    sizes.NumSampleTimes = 1;
    
    sys = simsizes(sizes);
    str = [];
    x0  = [0 0 0 0];   % [state time_skip time_uait2decide]
    ts  = [0.025 0];   % inherited sample time
    
elseif flag==2
    if NMF_threshold<=0,return,end
    u2=u;
    state=zeros(size(u));
    upbuff=0.4;
    downbuff=0.4;
    incr=0.05;
    
    for j=1:length(state)
        if u2(j)<max(Th(j),0)
            % set a base line
            u2(j)=0;
        end
        
        state(j)=x(j);
        
        if u2(j)-state(j)>upbuff
            % buffer of rising
            state(j)=state(j)+incr;
            %         elseif state(j)-u2(j)>downbuff
            %             % buffer of declining
            %             state(j)=state(j)-incr;
        elseif u2(j)==0 || state(j)-u2(j)>downbuff
            % to make sure it uill go doun to zero
            %             state(j)=0;
            state(j)=state(j)-incr;
        end
        
        % upper threshold
        if state(j)>1
            state(j)=1;
        elseif state(j)<0
            state(j)=0;
        end
    end
    
    % deal with the ambigious state (when both exclusive synergy > exth)
    exth=0.7;
    if state(1) ~= state(2) && state(1) > exth && state(2) > exth
        state(1)=x(1);
        state(2)=x(2);
    end
    if state(3) ~= state(4) && state(3) > exth && state(4) > exth
        state(3)=x(3);
        state(4)=x(4);
    end
    
    %deal with the situation when two exclusive synergy level is both max
%     if state(1)==1 && state(2)==1
%         if x(1)>x(2)
%             state(2)=1-incr;
%         else
%             state(1)=1-incr;
%         end
%     end
%     if state(3)==1 && state(4)==1
%         if x(3)>x(4)
%             state(4)=1-incr;
%         else
%             state(3)=1-incr;
%         end
%     end
    
    sys=state;
    
elseif flag==3
    sys = x;
end


