clc;clear all;close all;
% Load a song
 [d,sr] = audioread('test.wav');
 % Calculate the beat times
 b = beatdyn(d,sr);
 % Resynthesize a blip-track of the same length
 db = mkblips(b,sr,length(d));
 % Listen to them mixed together
 soundsc(d+db,sr);