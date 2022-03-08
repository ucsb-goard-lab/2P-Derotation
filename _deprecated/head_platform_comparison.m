% inject heading?
disp('Choose your floating:')
fn = uigetfile('*.mat');
load(fn);
disp('Choose supplementary file:')
load(uigetfile('*.mat'));

% corrected_heading = heading - nanmean(floating.heading);
% corrected_heading = wrapTo180(corrected_heading);
floating.alpha = heading;
floating.heading = heading;
save(fn, 'floating');
%%%%%%%


head = LightDarkAnalyzer();

plat = LightDarkAnalyzer();

save('processed_data.mat', 'head', 'plat')

%%%%%%%%%%%%%%%%%%%%%%%5
pref1 = head.calculatePreferredDirection('fit');
pref2 = plat.calculatePreferredDirection('fit');
% 
pref1 = pref1(:, 2);
pref2 = pref2(:, 2);

head.calculateHeadDirection_new();
plat.calculateHeadDirection_new();

% ishd = (pref1(:, 3) > 0.5 | pref2(:, 3) > 0.5) & (head.is_head_direction | plat.is_head_direction) & is_matched';
ishd =  is_matched;
% ishd = true(1, length(head.is_head_direction));
% ishd = is_matched';

figure
subplot(1, 2, 1)
scatter(pref1(ishd), pref2(ishd))
axis square
ylim([0, 2*pi])
xlim([0, 2*pi])
refline([1, 0])

subplot(1, 2, 2)
histogram(abs(angdiff(pref1(ishd), pref2(ishd))))

angular_distance = abs(angdiff(pref1(ishd), pref2(ishd)));
%%% checks


p = plat.tuning_curves;
h = head.tuning_curves;


x = linspace(-pi, pi, 60);
cmap = lines(2);
for ii = find(ishd)'
    plot(x, rescale(h(ii, :)), 'LineWidth', 2, 'Color', cmap(1, :))
    hold on
    plot(x, rescale(p(ii, :)), 'LineWidth', 2, 'Color', cmap(2, :))
    
    xline(pref1(ii)-pi, ':', 'LineWidth', 2, 'Color', cmap(1, :))
    xline(pref2(ii)-pi, ':', 'LineWidth', 2, 'Color', cmap(2, :))
    
    hold off
    legend({'Head rotation', 'Platform rotation'})
    ylabel('normalized \DeltaF/F')
    xlabel('heading')
    title(angular_distance(ii));
    prettyPlot
    
    
    pause
end

%%%%%%%%%%%


include = true(1, size(pref1, 1));
resp1 = head.tuning_curves(include, :);
resp2 = plat.tuning_curves(include, :);


shift_query = [0:1:59]; 
% phi_opt = zeros(1, size(resp1, 1));
phi_candidates = zeros(size(resp1, 1), length(shift_query));
for c = 1:size(resp1, 1)
    for s = 1:length(shift_query)
        phi_candidates(c, s) = corr(resp1(c, :)', circshift(resp2(c , :), shift_query(s))');
    end
end
[~, phi_opt] = max(phi_candidates, [], 2);


uniform_distribution = makedist('uniform', 'lower', 0, 'upper', 60);
[~, p] = chi2gof(phi_opt, 'cdf', uniform_distribution); 
fprintf('p value: %0.5f\n', p)


% plot
% figure
h = histogram(phi_opt, 10, 'Normalization', 'probability');








p = plat.tuning_curves;
h = head.tuning_curves;

for ii = 1:size(h, 1)
    hpcor(ii) = corr(h(ii, :)', p(ii, :)');
end

ishd = plat.calculateHeadDirection_trial();

for ii = find(ishd)'
    plot(rescale(h(ii, :)), 'LineWidth', 2)
    hold on
    plot(rescale(p(ii, :)), 'LineWidth', 2)
%     title(ii)
    hold off
    legend({'Head rotation', 'Platform rotation'})
    ylabel('normalized \DeltaF/F')
    xticklabels(-180:60:180)
    xlabel('heading')
    prettyPlot
%     export_fig(sprintf('example_cell_%d.png', ii), '-transparent', '-m8');
    pause
end
head_pref = head.calculatePreferredDirection('fit');
plat_pref = plat.calculatePreferredDirection('fit');




















% gauss?
for ii = 1:size(p, 1)
    p_fit{ii} = fit([1:60]', p(ii, :)', 'gauss1', 'Lower', [0, 0, 0], 'Upper', [Inf, 60, Inf]);
    h_fit{ii} = fit([1:60]', h(ii, :)', 'gauss1', 'Lower', [0, 0, 0], 'Upper', [Inf, 60, Inf]);
end
    
for ii = 1:size(p, 1)
    p_pref(ii) = p_fit{ii}.b1;
    h_pref(ii) = h_fit{ii}.b1;
    
    p_std(ii) = p_fit{ii}.c1;
    h_std(ii) = h_fit{ii}.c1;
end

for ii = find(p_std < 20)
    plot(p_fit{ii})
    hold on
    plot(p(ii, :))
    hold off
    pause
end

scatter(p_pref(p_std < 20), h_pref(p_std < 20))