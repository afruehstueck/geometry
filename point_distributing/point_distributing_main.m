% @file     point_distributing_main.m
% @author   afruehstueck
% @date     23/03/2017
%
% start the implementations of the three point distributing algorithms with
% custom specific parameters from this main file
% -> you can also start the implementations from the respective function 
% files, which will set the parameters to default values

close all
clc;
clear;

lloyd(128, 500, 1e-8, false)
dart_throwing(0.08, 200, true, 0.0)
ccvt(2048, 64*2048, 50, true, true, '../data/img/sloth_drawing.jpg')

