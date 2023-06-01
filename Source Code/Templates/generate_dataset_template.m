%% generate_dataset_template.m
% Denne templaten er for å vise hvordan generere datasett ved å bruke
% fieldii_generate_dataset-funksjonen implementert under Funksjoner.
% Denne metoden er bare en generell visning, og skal kun kopieres inn for å
% kunne implementeres. Navn på datasett må endres.

clear;
close all;

path = "Source Code\Field II Dataset";
filename = "Template";

degs_separated = -3;

phantom_positions(1,:)  = [20/1000*sind( degs_separated ), 0, 20/1000*cosd( degs_separated )];
phantom_positions(2,:)  = [-20/1000*sind( degs_separated ), 0, 20/1000*cosd( degs_separated )];
phantom_positions(3,:)  = [40/1000*sind( degs_separated ), 0, 40/1000*cosd( degs_separated )];
phantom_positions(4,:)  = [-40/1000*sind( degs_separated ), 0, 40/1000*cosd( degs_separated )];

h.degs_separated = degs_separated;

% Denne må settes dersom man ønsker å direkte lagre dataen. Dersom man
% ønsker et valg underveis settes ikke denne variabelen.
h.save = 1;


channel_data = fieldii_generate_dataset(path, filename, phantom_positions,h);
