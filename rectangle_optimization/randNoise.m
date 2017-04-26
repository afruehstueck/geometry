% @file     randNoise.m
% @author   afruehstueck
% @date     03/04/2017
%
% takes a float range value and returns a random value between (-noise/2, noise/2)

function randVal = randNoise(noise)
    randVal = (rand * noise) - (noise / 2);
end