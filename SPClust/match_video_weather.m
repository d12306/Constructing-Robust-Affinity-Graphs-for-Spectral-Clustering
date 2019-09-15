%% Video and weather data are not necessary to be synchronised, necessatating matching between them
%
% @Author: Eddy Zhu
% @Date: 20 Feb. 2013

clc;

addpath('../../');
global_vars;
global GAP_VIDEO_WEATHER;

% The pre-define ideal time gap between video and weather data
ideal_time_gap = 4; % unit: miniutes

training_set = '_2';
video_dir = ['/import/geb-experiments/Eddy/dataset/webcam/NY_TS/training' training_set '/'];
%video_dir = ['Z:/Eddy/dataset/webcam/NY_TS/training' training_set '/'];
%video_dir = ['/export/beware/thumper/Eddy/dataset/webcam/NY_TS/training' training_set '/']; % thumper
%video_dir = ['/homes/xz303/training' training_set '/']; % home dir

%weather_dir = '/export/beware/thumper/Eddy/dataset/webcam/NY_TS/weather/'; % thumper
%weather_dir = 'Z:/Eddy/dataset/webcam/NY_TS/weather/'; % experiments server
weather_dir = '/import/geb-experiments/Eddy/dataset/webcam/NY_TS/weather/'; % experiments server
tic;

video_folders = dir(video_dir);
video_folders(1) = [];
video_folders(1) = [];
display('-- video folder #: %d\n', length(video_folders));

weather_files = dir(weather_dir);
weather_files(1) = [];
weather_files(1) = [];
display('-- weather files #: %d\n', length(weather_files));

video_weather_pair = cell(length(video_folders), 4); % Format: video, weather, time gap (normal and abnormal)

% to match weather data for all video folders
for video_idx = 1 : length(video_folders)
	v_folder = video_folders(video_idx);
    if (v_folder.isdir == 1 && strcmp(v_folder.name(1:3), '201') == 1)
        
        target_file = [video_dir, v_folder.name '/weather_matched.xml'];
%         if(exist(target_file, 'file') > 0)
%             continue;
%         end

        dir_date_time_str = regexp(v_folder.name, '_', 'split');
        video_date_time = [regexp(dir_date_time_str{1}, '-', 'split') regexp(dir_date_time_str{2}, '-', 'split')];
        video_date_time = [str2double(video_date_time{1}) str2double(video_date_time{2}) str2double(video_date_time{3}) str2double(video_date_time{4}) str2double(video_date_time{5}) str2double(video_date_time{6})];

        min_time_gap = 100000 * 60; % for those unable to find weather
        closest_weather = 0;
        isfound = 0;
        
        for weather_idx = 1 : length(weather_files)
            w_file = weather_files(weather_idx);

            if (w_file.isdir == 0 && strcmp(w_file.name(end-10:end), 'weather.xml') == 1)
                stamp = regexp(w_file.name, '_weather', 'split');
                date_time_str = regexp(stamp{1}, '_', 'split');
                weather_date_time = [regexp(date_time_str{1}, '-', 'split') regexp(date_time_str{2}, '-', 'split')];
                weather_date_time = [str2double(weather_date_time{1}) str2double(weather_date_time{2}) str2double(weather_date_time{3}) str2double(weather_date_time{4}) str2double(weather_date_time{5}) str2double(weather_date_time{6})];
                
                % the time gap between two time points
                time_elps = etime(video_date_time, weather_date_time);
                time_elps = time_elps + (GAP_VIDEO_WEATHER * 60);

                if(abs(time_elps) <= ideal_time_gap * 60) % if the time gap between the video and weather is small enough, stop searching
                    % copy weather data
                    try
                        copyfile([weather_dir w_file.name], target_file, 'f');
                    catch me
                        %disp(me.message);
                    end
                    
                    video_weather_pair{video_idx, 1} = v_folder.name;
                    video_weather_pair{video_idx, 2} = w_file.name;
                    video_weather_pair{video_idx, 3} = time_elps;
                    
                    isfound = 1;
                    disp(['got it for video: ' v_folder.name, ': ' w_file.name]);
                    break;
                else
                    % record the weather with minimum time distance
                    if(abs(time_elps) < abs(min_time_gap))
                        min_time_gap = time_elps;
                        closest_weather = w_file.name;
                    end
                end
            end
        end
        
        if(isfound == 0) % no ideal weather data;
            disp(['Note **: fail to find a good weather data: ' v_folder.name]);
            % copy weather data
            try
                copyfile([weather_dir closest_weather], target_file, 'f');
            catch me
                %disp(me.message);
            end
            video_weather_pair{video_idx, 1} = v_folder.name;
            video_weather_pair{video_idx, 2} = closest_weather; 
            video_weather_pair{video_idx, 4} = min_time_gap/60; % unit: minute
        end      

    end

end

save('video_weather_pair.mat', 'video_weather_pair');
    
toc;
