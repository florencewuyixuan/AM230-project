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
colors = flipud(cbrewer('div', 'RdBu', 512));
C=linspecer(30);

% colors = cbrewer('seq', 'BuGn', 512);
% C=linspecer(30);

colors(colors<0) = 0;
colors(colors>1) = 1;
%% particle position
% fp = '/Users/dmitrykazakov/Dropbox (Harvard University)/courses/AM230/project/AM230-project/from Florence/J 2.5';
fp = '/Users/dmitrykazakov/Dropbox (Harvard University)/courses/AM230/project/AM230-project/from Florence/content';
% fp = '/Users/dmitrykazakov/Downloads/replicate';
% fp = '/Users/dmitrykazakov/Dropbox (Harvard University)/courses/AM230/project/AM230-project';
cd(fp)

files = dir('*.mat');
% files = dir('testing_*.mat');

pos_x_evolution = zeros(250,length(files));
pos_y_evolution = zeros(250,length(files));


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
xlim([-10 10])
ylim([-10 10])
 
for ii = 1:length(files)
    delete(particles)
    particles = plot(pos_x_evolution(:,ii),pos_y_evolution(:,ii),'.','MarkerSize',20,'Color',[0 0 0]);
    pause(0.01)
end


%% Compute MSD
close all

clear SD

% for ii = 1:length(files)/2
%     SD(:,ii) = (pos_x_evolution(:,length(files)/2+ii-1) - pos_x_evolution(:,length(files)/2)).^2 + (pos_y_evolution(:,length(files)/2+ii-1) - pos_y_evolution(:,length(files)/2)).^2;
%     
% end

for ii = 1:length(files)/2
    SD(:,ii) = (pos_x_evolution(:,ii) - pos_x_evolution(:,1)).^2 + (pos_y_evolution(:,ii) - pos_y_evolution(:,1)).^2;
    
end

MSD = mean(SD,1);

fz = 50;
ar = 2;


figure
set(gcf,'Position',[400   600   fz*ar   fz])
pbaspect([ar 1 1])
plot(MSD)
xlabel('Time step')
ylabel('MSD')
% set(gca,'Yscale','log')
set(gcf,'Color',[1 1 1])
set(gca,'Position',[0 0 1 1])
set(gca,'XColor','None')
set(gca,'YColor','None')



figure
set(gcf,'Position',[600   600   fz*ar   fz])
pbaspect([ar 1 1])
histogram(MSD,'Normalization','probability')
xlabel('Displacement')
ylabel('Probability')
set(gcf,'Color',[1 1 1])
% set(gca,'Position',[0 0 1 1])
% set(gca,'XColor','None')
% set(gca,'YColor','None')


%% intensity
close all

fp = '/Users/dmitrykazakov/Dropbox (Harvard University)/courses/AM230/project/AM230-project';
cd(fp)

files = dir('intensity*');

field_size = 2^9;

intensity_vs_time = zeros(field_size,field_size,length(files));


for ii = 1:length(files)
    filename = files(ii).name;
    step = split(filename,'_');
    step = split(step{2},'.');
    step = str2double(step{1});
    load(filename)

    intensity_vs_time(:,:,step) = intensity;


end
 

figure
hold on
int_plot = imagesc(intensity_vs_time(:,:,1));
% xlim([1 field_size])
% ylim([1 field_size])
colormap(colors)
xlim([150 350])
ylim([150 350])
axis square


 
for ii = 1:length(files)
    delete(int_plot)
    int_plot = imagesc(intensity_vs_time(:,:,ii));
    pause(0.2)
end


%% intensity and particles together

filename = 'motility_speckle.gif';

fig = figure;
hold on
set(gcf,'Color',[1 1 1])

subplot(1,2,1)
hold on
particles = plot(pos_x_evolution(:,1),pos_y_evolution(:,1),'.','MarkerSize',20,'Color',[0 0 0]);
xlim([-10 10])
ylim([-10 10])
axis square

subplot(1,2,2)
hold on
int_plot = imagesc(intensity_vs_time(:,:,1));
colormap(colors)

% xlim([1 field_size])
% ylim([1 field_size])
xlim([150 350])
ylim([150 350])
axis square


for ii = 1:length(files)
    frame = getframe(fig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
%     if ii == 1 
%         imwrite(imind,cm,filename,'gif','DelayTime',0.2,'Loopcount',inf); 
%     else 
%         imwrite(imind,cm,filename,'gif','DelayTime',0.2,'WriteMode','append');
%     end
    subplot(1,2,1)
    delete(particles)
    particles = plot(pos_x_evolution(:,ii),pos_y_evolution(:,ii),'.','MarkerSize',20,'Color',[0 0 0]);

    subplot(1,2,2)
    delete(int_plot)
    int_plot = imagesc(intensity_vs_time(:,:,ii));
    pause(0.2)

    
end

%% intensity only

filename = 'motility_speckle_1.gif';

fig = figure;
hold on
set(gcf,'Color',[1 1 1])

hold on
particles = plot(pos_x_evolution(:,1),pos_y_evolution(:,1),'.','MarkerSize',20,'Color',[0 0 0]);
xlim([-10 10])
ylim([-10 10])
axis square

% subplot(1,2,2)
% hold on
% int_plot = imagesc(intensity_vs_time(:,:,1));
% colormap(colors)

% xlim([1 field_size])
% ylim([1 field_size])
% xlim([150 350])
% ylim([150 350])
% axis square


for ii = 1:length(files)/2
    frame = getframe(fig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
%     Write to the GIF File 
    if ii == 1 
        imwrite(imind,cm,filename,'gif','DelayTime',0.2,'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','DelayTime',0.2,'WriteMode','append');
    end
    delete(particles)
    particles = plot(pos_x_evolution(:,length(files)/2+ii-1),pos_y_evolution(:,length(files)/2+ii-1),'.','MarkerSize',20,'Color',[0 0 0]);

%     subplot(1,2,2)
%     delete(int_plot)
%     int_plot = imagesc(intensity_vs_time(:,:,ii));
%     pause(0.2)

    
end




