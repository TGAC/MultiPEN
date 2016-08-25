function [chemicalID, mSamples, mLevels, mIdentif] = loadMetabolomics()

%% Import data from text file.
% Script for importing data from the following text file:
%
%    /Users/troncosp/Documents/Public-Datasets/2016-8-24_Fatty-Liver-Disease_Gene-Expression_Metabolomics/Processed-Data/MetaboLights/m_hna_fld_metabolite_profiling_NMR_spectroscopy_v2_maf.txt
%

%% Initialize variables.
filename = '/Users/troncosp/Documents/Public-Datasets/2016-8-24_Fatty-Liver-Disease_Gene-Expression_Metabolomics/Processed-Data/MetaboLights/m_hna_fld_metabolite_profiling_NMR_spectroscopy_v2_maf.txt';
delimiter = '\t';

%% chemicalID
startRow = 2;

% Format string for each line of text:
%   column1: text (%q)
formatSpec = '%q%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

% Close the text file.
fclose(fileID);

% Create output variable
chemicalID = [dataArray{1:end-1}];


%% mSamples
% Initialize variables.
endRow = 1;

% Format string for each line of text:
%   column19: text (%q)
%	. . .
%	column36: text (%q)
formatSpec = '%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, endRow, 'Delimiter', delimiter, 'ReturnOnError', false);

% Close the text file.
fclose(fileID);

% Create output variable
mSamples = [dataArray{1:end-1}];
mSamples = mSamples'; % for a 18-by-1 vector


%%  mLevels

% Initialize variables.
startRow = 2;

% Format string for each line of text:
%   column19: double (%f)
%   . . .
%	column36: double (%f)
formatSpec = '%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

% Close the text file.
fclose(fileID);

% Create output variable
mLevels = [dataArray{1:end-1}];
% Transpose mLevels as it originally has metabolites in rows and samples in
% columns
mLevels = mLevels';  %n-by-p



%% mIdentif

% Initialize variables.
startRow = 2;

% Format string for each line of text:
%   column5: text (%q)
formatSpec = '%*q%*q%*q%*q%q%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

% Close the text file.
fclose(fileID);

% Create output variable
mIdentif = [dataArray{1:end-1}];