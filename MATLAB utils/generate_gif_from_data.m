%% Init: Run before any other cell

% sys = 'Windows';
sys = 'Mac';

switch sys
    case 'Mac'
        prefix = '/Users/dmitrykazakov/';
    case 'Windows'
        prefix = 'C:/Users/dmk637/';
end

addpath([prefix 'Dropbox (Harvard University)/authored/TuringPatterns/cbrewer'])
addpath([prefix 'Dropbox (Harvard University)/courses/QUALS/linspecer'])
colors = cbrewer('div', 'RdBu', 512);
C=linspecer(30);

%%

fp = '/Users/dmitrykazakov/Dropbox (Harvard University)/courses/AM230/project/AM230-project';
cd(fp)

files = dir('*.mat');

pos_x_evolution = zeros(40,length(files));
pos_y_evolution = zeros(40,length(files));


% ii = 1;

for ii = 1:length(files)
    filename = files(ii).name;
    step = split(filename,'_');
    step = split(step{2},'.');
    step = str2double(step{1});
    load(filename)

    pos_x_evolution(:,step+1) = x';
    pos_y_evolution(:,step+1) = y';


end


figure
hold on
particles = plot(pos_x_evolution(:,1),pos_y_evolution(:,1),'.','MarkerSize',20,'Color',[0 0 0]);
 
for ii = 1:length(files)
    delete(particles)
    particles = plot(pos_x_evolution(:,ii),pos_y_evolution(:,ii),'.','MarkerSize',20,'Color',[0 0 0]);
    pause(0.1)
end