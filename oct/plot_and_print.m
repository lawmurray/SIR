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

    sz = [ 16 12 ];
    set (figure(1), 'papersize', sz);
    set (figure(1), 'paperposition', [0 0 sz]);
    set (figure(1), 'paperorientation', 'portrait');

    clf;
    subplot(4,4,1:2);
    bi_plot_quantiles('results/posterior.nc', 's');
    ylabel('S');
    axis([0 13]);
    grid on;
    
    subplot(4,4,3:4);
    bi_plot_quantiles('results/posterior.nc', 'i');
    hold on;
    bi_plot_paths('data/obs.nc', 'y_i');
    hold off;
    ylabel('I');
    axis([0 13]);
    grid on;
    
    subplot(4,4,5:6);
    bi_plot_quantiles('results/posterior.nc', 'r');
    xlabel('t');
    ylabel('R');
    axis([0 13]);
    grid on;

    subplot(4,4,7:8);
    bi_plot_quantiles('results/posterior.nc', 'ln_beta', 1);
    hold on;
    bi_plot_quantiles('results/posterior.nc', 'ln_beta', 2);
    hold off;
    
    xlabel('t');
    ylabel('log \beta, log \nu');
    axis([0 13]);
    grid on;

    subplot(4,3,7);
    bi_hist('results/posterior.nc', 'theta', [1,1]);
    xlabel('\theta_{\beta,1}');
    %axis('tight');
    grid on;

    subplot(4,3,8);
    bi_hist('results/posterior.nc', 'theta', [2,1]);
    xlabel('\theta_{\beta,2}');
    %axis('tight');
    grid on;

    subplot(4,3,9);
    bi_hist('results/posterior.nc', 'theta', [3,1]);
    xlabel('\theta_{\beta,3}');
    %axis('tight');
    grid on;

    subplot(4,3,10);
    bi_hist('results/posterior.nc', 'theta', [1,2]);
    xlabel('\theta_{\nu,1}');
    %axis('tight');
    grid on;

    subplot(4,3,11);
    bi_hist('results/posterior.nc', 'theta', [2,2]);
    xlabel('\theta_{\nu,2}');
    %axis('tight');
    grid on;

    subplot(4,3,12);
    bi_hist('results/posterior.nc', 'theta', [3,2]);
    xlabel('\theta_{\nu,3}');
    %axis('tight');
    grid on;
    
    file = sprintf('%s/posterior.pdf', FIG_DIR);
    saveas(figure(1), file);
    system(sprintf('pdfcrop %s %s', file, file));
end
