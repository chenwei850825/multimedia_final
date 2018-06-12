%all 2 channels
clc;clear all;close all;
[y_input,fs]=audioread('input.mp3'); 

start_offset = 176216; % ms
end_offset = 178968; % ms
BPM = 109;
beat_snap = 8;

cut_input = y_input( floor(start_offset * (fs/1000) + 1) : floor(end_offset * fs/1000) + 1, : );
beat_snap_ms = 60/BPM/beat_snap*1000; % ms range you want to judge

player = audioplayer(cut_input,fs);
play(player);

%figure(1);
subplot(1,1,1);
plot(cut_input);

diff_per_ms = [];
sum_diff_per_ms = [];
prev = [0 0];
diff = [0 0];
now = [0 0];
for i = 1 : length(cut_input(:,1))

    if mod( i , fs/1000) < 0.999999
            diff_per_ms = [diff_per_ms ; diff(1)/(fs/1000) diff(2)/(fs/1000)];
            sum_diff_per_ms = [sum_diff_per_ms ; diff_per_ms(end,1)+diff_per_ms(end,2) ];
            diff = [0 0];
    end

    for j = 1:2
        now( j ) = cut_input( i , j);
        diff( j ) = diff( j ) + abs(prev( j ) - now( j ));
        prev( j ) = now( j );
    end
    
end

%figure(2);
subplot(1,1,1);
plot(diff_per_ms);

%figure(3);
subplot(1,1,1);
plot(sum_diff_per_ms);

fileID = fopen( 'test2.txt' , 'w' );
for i = 1 : length(diff_per_ms(:,1))
    fprintf(fileID, '%d,%f\r\n', i -1, sum_diff_per_ms(i) ); % i - 1  for ms
end
fclose(fileID);

% choose points above  a threshold every  ( judge_time )ms from sum_diff_per_ms

fileID = fopen( 'test3.txt' , 'w' );

for i = 1:length(sum_diff_per_ms)
    if mod(i, beat_snap_ms) < 0.999999

        max = 0;
        max_ms = 0;
        for j = i - floor(0.5*beat_snap_ms) : i + floor(0.5*beat_snap_ms)
        
            if j < length(sum_diff_per_ms) && max < sum_diff_per_ms( j )
                max_ms = j;
                max = sum_diff_per_ms( j );
            end
            
        end

        fprintf(fileID, '256,256,%d,1,2,0:0:0:0:\r\n', start_offset + max_ms );
        
    end
end

fclose(fileID);