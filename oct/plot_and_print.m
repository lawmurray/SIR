% Copyright (C) 2013
% Author: Lawrence Murray <lawrence.murray@csiro.au>
% $Rev$
% $Date$

% -*- texinfo -*-
% @deftypefn {Function File} plot_and_print ()
%
% Produce plots and print for manuscript.
% @end deftypefn
%
function plot_and_print ()
    FIG_DIR = 'figs';
    mkdir(FIG_DIR);

    sz = [ 12 9 ];
    set (figure(1), 'papersize', sz);
    set (figure(1), 'paperposition', [0 0 sz]);
    set (figure(1), 'paperorientation', 'portrait');

    clf;
    subplot(4,4,1:4);
    bi_plot_quantiles('results/posterior.nc', 's');
    ylabel('S');
    axis('tight');
    grid on;
    
    subplot(4,4,5:8);
    bi_plot_quantiles('results/posterior.nc', 'i');
    hold on;
    bi_plot_paths('data/obs.nc', 'y_i');
    hold off;
    ylabel('I');
    axis('tight');
    grid on;
    
    subplot(4,4,9:12);
    bi_plot_quantiles('results/posterior.nc', 'r');
    xlabel('t');
    ylabel('R');
    axis('tight');
    grid on;
    
    subplot(4,4,13);
    bi_hist('results/posterior.nc', 'beta');
    xlabel('\beta');
    axis('tight');
    grid on;

    subplot(4,4,14);
    bi_hist('results/posterior.nc', 'nu');
    xlabel('\nu');
    axis('tight');
    grid on;
 
    subplot(4,4,15);
    bi_hist('results/posterior.nc', 'sigma1');
    xlabel('\sigma_1');
    axis('tight');
    axis('tight');
    grid on;
 
    subplot(4,4,16);
    bi_hist('results/posterior.nc', 'sigma2');
    xlabel('\sigma_2');
    axis('tight');
    grid on;

    file = sprintf('%s/posterior.pdf', FIG_DIR);
    saveas (figure(1), file);
    system(sprintf('pdfcrop %s %s', file, file));
end
