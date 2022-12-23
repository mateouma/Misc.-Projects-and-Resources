%% download x_data_processed mat file
load('x_data_processed.mat')
clearvars -except m2_MD_miniscope_frames dff

frames_per_sec = 15; % don't forget temporal downsampling

%% Frame load
% convert behavior timestamps to miniscope timestamps
start_frame = m2_MD_miniscope_frames(:, 1);
choice_frame = m2_MD_miniscope_frames(:, 2);
end_frame = m2_MD_miniscope_frames(:, 3);

%temporal downsampling
start_frame = round(start_frame / 2);
choice_frame = round(choice_frame / 2);
end_frame = round(end_frame / 2);


%% prompt user to enter data type

prompt = "How long do you want the interval to be? (x seconds before and after event)\n";
seconds = input(prompt);

while true
    prompt = "Which event do you want? Type 1, 2, or 3 \n 1: Start \n 2: Choice \n 3: End \n";
    data_type = input(prompt);
    
    % enter which timestamps you want to look at
    if data_type == 1
        event = start_frame;
        event_name = "Maze Start";
        break
    elseif data_type == 2
        event = choice_frame;
        event_name = "Choice";
        break
    elseif data_type == 3
        event = end_frame;
        event_name = "Food Cup";
        break
    end
end






%% Calculate average activity x seconds before and after an event
% beginning = (frames x seconds) before event
% ending = (frames x seconds) after event
%% OLD
beginning = event - (seconds*frames_per_sec)+1;
ending = event + (seconds*frames_per_sec);

%% NEW KYLE
beginning = event - (seconds*frames_per_sec);
ending = event + (seconds*frames_per_sec);

% number neurons x event length (seconds*frames_per_sec*2)
% ex. 349 neurons (rows) x 60 frames = 2 seconds (columns)
neuron_activity = zeros([size(dff,1), seconds*frames_per_sec*2 + 1]);




%% loop through each neuron to gather data
for neuron = 1:size(dff,1)
    % for one neuron
    % number trials x event length
    one_neuron = zeros([length(start_frame), seconds*frames_per_sec*2] + 1);
    
    % loop through each trial, normalize neuron signal
    for trial = 1:length(beginning)
        one_neuron(trial,:) = normalize(dff(neuron, beginning(trial):ending(trial)));
    end

    % take average of single neuron across all trials
    neuron_activity(neuron,:) = mean(one_neuron);
end



%% sorting
% sort neurons by peak activity
Arowmax = max(neuron_activity, [], 2);
[~,idx] = sort(Arowmax, 'descend');
neuron_activity_sorted = neuron_activity(idx,:);


% sort neurons by peak activity attempt 2
[~, max_pos] = max(neuron_activity, [], 2);
[~, max_order] = sort(max_pos);
nas2 = neuron_activity(max_order, :);

% sort new_data by peak activity
% Arowmax = max(new_data, [], 2);
% [~,idx] = sort(Arowmax, 'descend');
% new_data_sorted = new_data(idx,:);

%% plotting
%plot unsorted neurons
figure;
imagesc(neuron_activity);
colormap(jet);
hold on;
ylabel("Neuron Index");
ticksarray=[1 seconds*frames_per_sec seconds*frames_per_sec*2];
tickslabels={"-"+seconds+" seconds", event_name, "+"+seconds+" seconds"};
set(gca,'XTick',ticksarray)
set(gca,'XTickLabel',tickslabels)
colorbar;
caxis([-1.8, 3.2]);
title("Trial-averaged neural activity around " + event_name);

% plot sorted neurons
figure;
imagesc(nas2); %imagesc(neuron_activity_sorted);
colormap(jet);
hold on;
ylabel("Sorted Neurons");
ticksarray=[1 seconds*frames_per_sec seconds*frames_per_sec*2];
tickslabels={"-"+seconds+" seconds", event_name, "+"+seconds+" seconds"};
set(gca,'XTick',ticksarray)
set(gca,'XTickLabel',tickslabels)
colorbar();
caxis([-1.8, 3.2]);
title("Trial-averaged neural activity around " + event_name + " (sorted by frame of max activity)");


% %% Individual Neurons
% row = 10;
% col = row;
% n = 1; %useful couter
% full_num = 100; % 100 unless last iteration
% 
% while true % loop until reach end (full_num != 0)
%     
%     extra = (n-1) * 100;    
%     figure;    
%     for i = 1:100 % 1:49 for 301-349
%         
%         if ((i + extra) > size(neuron_activity, 1))
%             full_num = i - 1;
%             break
%         end
%         
%         subplot(row, col, i);
% 
%         plot(neuron_activity(i+extra, :));
%         ylim([0, 0.4]);
% 
%         hold on;
% 
%         % add midline
%         neuron_max = max(neuron_activity(i+extra, :));
%         t = "Neuron " + (i+extra);
% 
%         title(t);
%     end
%         
%     sgtitle("Neurons " + (extra + 1) + "-" + (extra + full_num)); %+ " (min of each neuron is 0)");
%     
%     if (full_num ~= 100)
%         break
%     end
%     n = n + 1;
% end














%% Significance test without normalization
num_shuffles = 10000;

shuffled_data = zeros(num_shuffles, (2 * seconds*frames_per_sec + 1));

significant_neuron_index = [];

for neuron = 1:size(dff,1)
    for i = 1:num_shuffles
        % random center frame 
        rand_center_frame = (seconds*frames_per_sec + 1) + ceil((size(dff, 2) - 2 * seconds*frames_per_sec - 1) * rand(1));
        
        % ensure minimum start_frame is 1
        rand_start_frame = rand_center_frame - seconds*frames_per_sec;

        rand_end_frame = rand_center_frame + seconds*frames_per_sec;
    
        shuffled_data(i,:) = dff(neuron, rand_start_frame:rand_end_frame);
    end

    
    one_neuron = zeros([length(start_frame), seconds*frames_per_sec*2 + 1]);
    
    % extract correct neural activity
    for trial = 1:length(beginning)
        one_neuron(trial,:) = dff(neuron, beginning(trial):ending(trial));
    end

    % take average of single neuron across all trials
    neuron_data = mean(one_neuron);

    %shuffled_data
    %neuron_data

    sorted_average_shuffled_data = sort(mean(shuffled_data));
    %sorted_average_shuffled_data(ceil(0.95 * length(sorted_average_shuffled_data)))
    
    if (mean(neuron_data) > sorted_average_shuffled_data(ceil(0.95 * length(sorted_average_shuffled_data))))
        significant_neuron_index(end+1) = neuron;
    end    
end

%% Significance test with normalization
% num_shuffles = 10000;
% 
% shuffled_data = zeros(num_shuffles, (2 * seconds*frames_per_sec + 1));
% 
% significant_neuron_index = [];
% 
% for neuron = 1:size(dff,1)
%     for i = 1:num_shuffles
%         % random center frame 
%         rand_center_frame = (seconds*frames_per_sec + 1) + ceil((size(dff, 2) - 2 * seconds*frames_per_sec - 1) * rand(1));
%         
%         % ensure minimum start_frame is 1
%         rand_start_frame = rand_center_frame - seconds*frames_per_sec;
% 
%         rand_end_frame = rand_center_frame + seconds*frames_per_sec;
%     
%         shuffled_data(i,:) = normalize(dff(neuron, rand_start_frame:rand_end_frame)); % this is wrong
%     end
% 
%     
%     one_neuron = zeros([length(start_frame), seconds*frames_per_sec*2 + 1]);
%     
%     % extract correct neural activity
%     for trial = 1:length(beginning)
%         one_neuron(trial,:) = normalize(dff(neuron, beginning(trial):ending(trial)));
%     end
% 
%     % take average of single neuron across all trials
%     neuron_data = mean(one_neuron);
% 
%     %shuffled_data
%     %neuron_data
% 
%     sorted_average_shuffled_data = sort(mean(shuffled_data));
%     %sorted_average_shuffled_data(ceil(0.95 * length(sorted_average_shuffled_data)))
%     
%     if (mean(neuron_data) > sorted_average_shuffled_data(ceil(0.95 * length(sorted_average_shuffled_data))))
%         significant_neuron_index(end+1) = neuron;
%     end    
% end


%% Neuron activity ignoring non-signficant neurons

neuron_activity_significant_separated = zeros([size(dff,1), seconds*frames_per_sec*2 + 1]);

%% loop through significant neurons to gather data
for neuron = 1:size(dff,1)
    if (any(significant_neuron_index(:) == neuron))
        % for one neuron
        % number trials x event length
        one_neuron = zeros([length(start_frame), seconds*frames_per_sec*2] + 1);
        
        % loop through each trial, normalize neuron signal
        for trial = 1:length(beginning)
            one_neuron(trial,:) = normalize(dff(neuron, beginning(trial):ending(trial)));
        end
    
        % take average of single neuron across all trials
        neuron_activity_significant_separated(neuron,:) = mean(one_neuron);
    else
        % Turn 0 to -1.5
        neuron_activity_significant_separated(neuron,:) = -1.5;
    end
end




%% sorting according to first
% sort neurons by peak activity
Arowmax = max(neuron_activity, [], 2);
[~,idx] = sort(Arowmax, 'descend');
neuron_activity_sorted = neuron_activity(idx,:);


% sort neurons by peak activity attempt 2
[~, max_pos] = max(neuron_activity, [], 2);
[~, max_order] = sort(max_pos);
nas_new = neuron_activity_significant_separated(max_order, :);





%% plotting
%plot unsorted neurons
figure;
imagesc(neuron_activity_significant_separated);
colormap(jet);
hold on;
ylabel("Neuron Index");
ticksarray=[1 seconds*frames_per_sec seconds*frames_per_sec*2];
tickslabels={"-"+seconds+" seconds", event_name, "+"+seconds+" seconds"};
set(gca,'XTick',ticksarray)
set(gca,'XTickLabel',tickslabels)
colorbar;
caxis([-1.8, 3.2]);
title("SIGNIFICANT NEURONS: Trial-averaged neural activity around " + event_name);



% plot sorted neurons
figure;
imagesc(nas_new); %imagesc(neuron_activity_sorted);
colormap(jet);
hold on;
ylabel("Sorted Neurons");
ticksarray=[1 seconds*frames_per_sec seconds*frames_per_sec*2];
tickslabels={"-"+seconds+" seconds", event_name, "+"+seconds+" seconds"};
set(gca,'XTick',ticksarray)
set(gca,'XTickLabel',tickslabels)
colorbar();
caxis([-1.8, 3.2]);
title("SIGNIFICANT NEURONS: Trial-averaged neural activity around " + event_name + " (sorted by frame of max activity)");



