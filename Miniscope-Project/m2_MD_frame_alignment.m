clearvars
miniscope = readtable("m2_MD_timeStamps.csv");

beh = readtable("m2_MD_timeStamps_beh.csv");

miniscope = miniscope{:,:}; % turn table into matrix
beh = beh{:,:};

% Kyle start frame
start_frame = [480, 3641, 5588, 8333, 10428, 13484, 16655, 19080, 21621,...
    23873, 26002, 28298, 30396, 32711];

% Kyle choice frame
choice_frame = [508, 3671, 5612, 8358, 10449, 13509, 16683, 19101, 21645,...
    23905, 26039, 28325, 30418, 32734];

% Kyle end frame
end_frame = [608, 3740, 5679, 8426, 10518, 13758, 16748, 19185, 21719,...
    23985, 26115, 28425, 30484, 32804];


% behavior timestamp
beh_timestamp_start = [];
beh_timestamp_choice = [];
beh_timestamp_end = [];

%% checkpoint 1
for (i = 1:length(start_frame))
    %index = find(start_frame(i) == beh(:, 1));
    beh_timestamp_start(i) = beh(find(start_frame(i) == beh(:, 1)), 2);
    beh_timestamp_choice(i) = beh(find(choice_frame(i) == beh(:, 1)), 2);
    beh_timestamp_end(i) = beh(find(end_frame(i) == beh(:, 1)), 2);
end


% Miniscope index
miniscope_start_index = [];
miniscope_choice_index = [];
miniscope_end_index = [];

for i = 1:length(beh_timestamp_start)
    % start
    [temp, start_index] = min(abs(miniscope(:,2) - beh_timestamp_start(i))); % temp has minimum number, index has index of min number
    miniscope_start_frame(i) = miniscope(start_index, 1);
    miniscope_start_index(length(miniscope_start_index)+1) = start_index;
    
    % choice
    [temp, choice_index] = min(abs(miniscope(:,2) - beh_timestamp_choice(i))); % temp has minimum number, index has index of min number
    miniscope_choice_frame(i) = miniscope(choice_index, 1);
    miniscope_choice_index(length(miniscope_choice_index)+1) = choice_index;
    
    % end
    [temp, end_index] = min(abs(miniscope(:,2) - beh_timestamp_end(i))); % temp has minimum number, index has index of min number
    miniscope_end_frame(i) = miniscope(end_index, 1);
    miniscope_end_index(length(miniscope_end_index)+1) = end_index;
end

m2_MD_miniscope_frames = [miniscope_start_frame', miniscope_choice_frame', miniscope_end_frame'];


