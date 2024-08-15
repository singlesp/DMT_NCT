function [m,sem] = plot_bounded_line(data,color)

%assumes data is subjects by time
%color is either text or [R G B] format

m = nanmean(data);
sem = std(data) / sqrt(size(data,1));

[l,p] = boundedline(1:size(data,2),m, sem,'alpha', 'transparency', 0.28);
set(l, 'linewidth', 2, 'color', color);
set(p, 'facecolor', color);