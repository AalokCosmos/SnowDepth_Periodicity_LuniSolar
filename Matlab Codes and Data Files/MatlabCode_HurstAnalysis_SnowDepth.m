clear; clc; close all;


%==================================================
% Load the Data

snowDepth_data = readtable('Data/SnowDepth_Data_700Days_NsnowN.csv', 'NumHeaderLines', 0);
snowDepth_data = snowDepth_data{:, :};

observation_num = snowDepth_data(:, 1);
snow_depth = snowDepth_data(:, 2);


%==================================================
% ASSIGN WINDOW AS TIME SERIES


X = snow_depth;
X = transpose(X);
X = diff(X);

figure(1);
plot(1:length(X), X);


%==================================================
% IF THE TIME SERIES IS NOISE LIKE CONVERT TO RANDOM WALK LIKE
% IF ITS RANDOM WALK LIKE DO NOT CHANGE ANYTHING

RW_X = cumsum(X - mean(X));

RW_scale = 2 / (max(RW_X) / max(X));

figure(2);
hold on
plot(1:length(X), X);
plot(1:length(X), RW_X * RW_scale, 'LineWidth', 1.5, 'Color', 'red');


%==================================================
% RMS OF THE WHOLE SERIES

mean_X = mean(X);
RMS_X = rms(X);

figure(3);
hold on
plot(1:length(X), X);
plot(1:length(X), ones(length(X), 1) * mean_X, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'red');
plot(1:length(X), mean_X + ones(length(X), 1) * RMS_X, 'LineWidth', 2, 'Color', 'red');
plot(1:length(X), mean_X - ones(length(X), 1) * RMS_X, 'LineWidth', 2, 'Color', 'red');


%==================================================
% LOCAL DETRENDING OF THE TIME SERIES

figure(4);
tiledlayout(3, 1);

m_list = [1, 2, 3];

for i = 1:3
    % Plot local fluctuations for given m
    nexttile
    hold on
    plot(1:length(X), X);
    plot(1:length(X), RW_X * RW_scale, 'LineWidth', 1.5, 'Color', 'red');

    scale = 2048;
    m = m_list(i);
    segments = floor(length(RW_X) / scale);

    for v = 1:segments
        m_Idx_start = ((v-1) * scale) + 1;
        m_Idx_stop = v * scale;
        m_Index{v} = m_Idx_start:m_Idx_stop;
        m_X_Idx = RW_X(m_Index{v});
        m_C = polyfit(m_Index{v}, RW_X(m_Index{v}), m);
        m_fit{v} = polyval(m_C, m_Index{v});
        m_RMS{1}(v) = sqrt(mean((m_X_Idx - m_fit{v}).^2));

        plot(m_Index{v}, m_fit{v} * RW_scale, 'LineStyle', '--', 'Color', 'black');
        plot(m_Index{v}, (m_fit{v} + m_RMS{1}(v)) * RW_scale, 'LineStyle', '-', 'Color', 'black');
        plot(m_Index{v}, (m_fit{v} - m_RMS{1}(v)) * RW_scale, 'LineStyle', '-', 'Color', 'black');

    end
    
    m_F(i) = sqrt(mean(m_RMS{1}.^2));

end


%==================================================
% MONOFRACTAL DETRENDED FLUCTUATION ANALYSIS

scale = [16, 32, 64, 128, 256, 512, 1024, 2048];
m = 1;
scale = flip(scale);

figure(5);
tiledlayout(length(scale), 1);

for ns = 1:length(scale)
    segments(ns) = floor(length(X) / scale(ns));
    
    for v = 1:segments(ns)
        Idx_start = ((v-1) * scale(ns)) + 1;
        Idx_stop = v * scale(ns);
        Index{v, ns} = Idx_start:Idx_stop;
        X_Idx = RW_X(Index{v, ns});
        C = polyfit(Index{v, ns}, RW_X(Index{v, ns}), m);
        fit{v, ns} = polyval(C, Index{v,ns});
        RMS{ns}(v) = sqrt(mean((X_Idx - fit{v, ns}).^2));
        
        RMS_display{ns}(Index{v, ns}) = RMS{ns}(v) * ones(1, length(Index{v, ns}));
    end
    
    F(ns) = sqrt(mean(RMS{ns}.^2));
    
    % Plot RMS of each segment and F for that scale
    nexttile
    hold on
    plot(1:length(RMS_display{ns}), RMS_display{ns}, 'LineWidth', 1.5);
    plot(1:length(RMS_display{ns}), ones(1, length(RMS_display{ns})) * F(ns), 'LineWidth', 1.5);
    ylabel(sprintf('Scale = %d', scale(ns)));
    
end


%==================================================
% PLOT DFA LINE

C = polyfit(log2(scale), log2(F), 1);
H = C(1);
RegLine = polyval(C, log2(scale));

figure(6);
hold on
plot(log2(scale), log2(F), 'LineStyle', 'none', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black');
plot(log2(scale), RegLine, 'LineWidth', 1.5, 'Color', 'blue');

xlim([log2(min(scale)), max(log2(scale))]);
ylim([min([min(log2(F)), min(RegLine)]), max([max(log2(F)), max(RegLine)])]);


disp(H);
% Clear the variables created
clear scale m ns segments v Idx_start Idx_stop Index X_Idx C fit RMS RMS_display  F H RegLine


%==================================================
% MULTIFRACTAL DETRENDED FLUCTUATION ANALYSIS OF TIME SERIES
% FOR PARTICULAR SCALE

scale = 16;
m = 2;
segments = floor(length(RW_X) / scale);

for v = 1:segments
    Idx_start = ((v-1) * scale) + 1;
    Idx_stop = v * scale;
    Index{v} = Idx_start:Idx_stop;
    X_Idx = RW_X(Index{v});
    C = polyfit(Index{v}, RW_X(Index{v}), m);
    fit{v} = polyval(C, Index{v});
    RMS{1}(v) = sqrt(mean((X_Idx - fit{v}).^2));
    
    RMS_display{1}(Index{v}) = RMS{1}(v) * ones(1, length(Index{v}));
end


q = [-3, -1, 1, 3];

figure(7);
tiledlayout(length(q) + 1, 1);

nexttile
plot(1:length(X), X);

for nq = 1:length(q)
    qRMS{1} = RMS{1}.^q(nq);
    q_RMS_display{1} = RMS_display{1}.^q(nq);
    Fq(nq) = mean(qRMS{1}).^(1/q(nq));
    
    nexttile
    plot(1:length(q_RMS_display{1}), q_RMS_display{1}, 'LineWidth', 1.5);

end
Fq(q == 0) = exp(0.5 * mean(log(RMS{1}.^2)));


%==================================================
% MULTIFRACTAL DETRENDED FLUCTUATION ANALYSIS OF TIME SERIES

clear scale m v nq segments q Idx_start Idx_stop Index X_Idx C fit RMS RMS_display qRMS Fq;

% scale = [16, 32, 64, 128, 256, 512, 1024];
scale_min = 32;
scale_max = 4096;
scale_res = 19;
scale_exponents = linspace(log2(scale_min), log2(scale_max), scale_res);
scale = round(2.^scale_exponents);

m = 3;
scale = flip(scale);
q = [-3, -1, 0, 1, 3];

for ns = 1:length(scale)
    segments(ns) = floor(length(X) / scale(ns));
    
    for v = 1:segments(ns)
        Idx_start = ((v-1) * scale(ns)) + 1;
        Idx_stop = v * scale(ns);
        Index{v, ns} = Idx_start:Idx_stop;
        X_Idx = RW_X(Index{v, ns});
        C = polyfit(Index{v, ns}, RW_X(Index{v, ns}), m);
        fit{v, ns} = polyval(C, Index{v,ns});
        RMS{ns}(v) = sqrt(mean((X_Idx - fit{v, ns}).^2));
        
        RMS_display{ns}(Index{v, ns}) = RMS{ns}(v) * ones(1, length(Index{v, ns}));
    end
    
    for nq = 1:length(q)
        qRMS{nq, ns} = RMS{ns}.^q(nq);
        Fq(nq, ns) = mean(qRMS{nq, ns}).^(1/q(nq));
    end
    Fq(q == 0, ns) = exp(0.5 * mean(log(RMS{ns}.^2)));
    
end


%==================================================
% PLOT MFDFA LINE

figure(8);
hold on

for nq = 1:length(q)
    C = polyfit(log2(scale), log2(Fq(nq,:)), 1);
    Hq(nq) = C(1);
    qRegLine{nq} = polyval(C, log2(scale));
    
    % Plot the MFDFA lines
    plot(log2(scale), log2(Fq(nq, :)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black');
    plot(log2(scale), qRegLine{nq}, 'LineWidth', 1.5, 'Color', 'blue');
    
end

return
% xlim([log2(min(scale)), max(log2(scale))]);
% ylim([0, max([max(max(log2(Fq))), max(max(qRegLine))])]);

































