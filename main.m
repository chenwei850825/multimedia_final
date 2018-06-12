clc;clear all;close all;
addpath(genpath('C:\Users\ycw\Desktop\multimedia_final\mp3_toolbox'))
addpath(genpath('C:\Users\ycw\Desktop\multimedia_final\MATLAB_TSM-Toolbox_2.01'))

begin = 162864; % Define where you want to start playing
now = 163684 ;
ending = 135172;
%start = begin;
%startSecond = start / 1000;
beginBeat = 0;
SPB = 60/(208*4);


[origin, sr] = audioread('test.mp3'); % Get sampling rate



fid = fopen('beatTime2.txt');
while ~feof(fid)
%{
    Str = fgetl(fid);
    Key   = '256,';
    Index = strfind(Str, Key);
    now = sscanf(Str(Index(2) + length(Key):end), '%g', 1)
    [begin, origin, sr] = PhaseVocoderMATLABExample(beging, now, origin, sr);
%}
    Str = fgetl(fid);
    Key   = '256,';
    Index = strfind(Str, Key);
    now = sscanf(Str(Index(2) + length(Key):end), '%g', 1)
    Key   = ': ';
    Index = strfind(Str, Key);
    nowBeat = sscanf(Str(Index(1) + length(Key):end), '%g', 1);
    
    beginSecond = begin / 1000; % Define where you want to start playing
    nowSecond = now / 1000;
    
    %beginBeat = floor((beginSecond - startSecond) / SPB);
    %nowBeat = floor((nowSecond - startSecond) / SPB);
    beatDiff =  nowBeat - beginBeat
    

    beginIndex = floor(beginSecond*sr); % Find beginning index of where to sample
    nowIndex = floor(nowSecond*sr);
    y = origin(beginIndex:nowIndex, :); % Extract out the signal from the starting point to the end of the file
    
    if (nowSecond - beginSecond) == 0
        continue;
    else
        [tt,sideinfo] = wsolaTSM(y, (SPB * beatDiff) /  (nowSecond - beginSecond) );
    end

    
    if 1 / ((nowSecond - beginSecond) / (SPB * beatDiff)) <= 1
        origin(beginIndex : beginIndex+length(tt) - 1, :) = tt;
        origin(beginIndex+length(tt): nowIndex, :) = [];
    else
        diff = length(tt) - (nowIndex - beginIndex + 1);
        tmp = zeros(length(origin) + diff, 2);
        tmp(1:beginIndex - 1, :) = origin(1:beginIndex - 1, :);
        tmp(beginIndex:beginIndex+length(tt)-1, :) = tt;
        tmp(beginIndex+length(tt):end, :) = origin(nowIndex + 1:end, :);
        origin = tmp;
    end
    
    beginBeat = nowBeat;
    begin = now;
end

%mp3write(origin,sr,32,'audio.mp3');
mp3write(origin,sr,32,'test2.mp3');

fclose(fid);
