%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: An initialization function for the block "Cue_Presentation"
%--------------------------------------------------------------------------------------------------------------------------------
% Used in: Simultaneous_Wrist_and_Hand.mdl
%--------------------------------------------------------------------------------------------------------------------------------
% Update : 2015/2/12 by Jiaxin MA @Kyoto Univ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('Cue_Presentation Initiating');
Cue_Presentation_GUI;
load('Provided_Cues.mat'); % the file containing the sequence of cues
movement_list=[1 1; 2 2; 3 3; 4 12; 5 13];
%time for ready/begin/hold/return
t_step=[2 2 3 2];

spr=600;

total_time = length(Provided_Cues)/spr;
% Find out how long the cues are...
different_cue = Provided_Cues(2,1);
counter = 1;
while different_cue == Provided_Cues(2,counter)
    counter = counter+1;
end
cue_length = (counter-1) / spr;
if cue_length ~= sum(t_step)
    error('Please check the length of cues or change t_step.');
end

%Determine the image sequence that will be shown corresponding to the desired cues.
image_sequence = zeros(total_time,1);
transition_image_sequence = zeros(total_time,1);
for i = 1:total_time-1
    image_sequence(i) = movement_list(Provided_Cues(2,(i-1)*spr+1),2);
    if Provided_Cues(2,(i-1)*spr+spr-1)~=Provided_Cues(2,(i*spr)+1)
        transition_image_sequence(i) = 1;
    end
end
transition_image_sequence(total_time) = 1;
image_sequence(end) = image_sequence(end-1);
hint_sequence = repmat([ones(t_step(1),1);2*ones(t_step(2),1);...
    3*ones(t_step(3),1);4*ones(t_step(4),1)],total_time/cue_length,1);
sequence = [image_sequence transition_image_sequence hint_sequence];

save('sequence.mat','sequence');
