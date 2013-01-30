
% Test whether it is necessary to rotate VB PCA
function nc2010_publish_datasets

n = 200;
m = 50;
d = 10;

ncomps = 30;

randn('state', 1);
rand('state', 1);

% Generate multivariate normal data
disp('Using weak subspace')
eigs = (1 + [d:-1:1 zeros(1,m-d)]) .^ 2
datastring = 'weak';
plot_eigenvalues(eigs, datastring);

disp('Using strong subspace')
eigs = ([5*ones(1,d) 1*ones(1,m-d)]) .^ 2
datastring = 'strong';
plot_eigenvalues(eigs, datastring);

disp('Using no subspace')
eigs = (m:-1:1) .^ 2
datastring = 'no';
plot_eigenvalues(eigs, datastring);


function plot_eigenvalues(eigs, datastring)

filename = sprintf(['/home/jluttine/papers/neurocomputing2010/' ...
                    'fig_eigenvalues_dataset=%s'], datastring);

figure
plot(sqrt(eigs), 'k-');
set(gcf, 'units', 'centimeters', 'paperunits', 'centimeters');
pos = get(gcf, 'position');
set(gcf, 'position', [pos(1:2), 6,4]);
pos = get(gcf, 'paperposition');
set(gcf, 'paperposition', [pos(1:2),6,4])

ylim([0, max(sqrt(eigs))+1]);
xlim([1, length(eigs)]);

%ylabel('standard deviation')

print(gcf, '-depsc2', filename);

