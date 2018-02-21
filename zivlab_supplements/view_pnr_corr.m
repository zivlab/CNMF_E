function view_pnr_corr(obj)
%% view PNR, correations that emerge from parameters used in CNMF-E
gSiz = obj.options.gSiz;
min_corr = obj.options.min_corr;
min_pnr = obj.options.min_pnr;

%%

fprintf('computing the correlation image and the peak-to-noise ratio image....\n');
[cn, pnr] = obj.correlation_pnr_parallel([1, 5000]);
fprintf('Done\n');

tmp_d = max(1,round(gSiz/4));
v_max = ordfilt2(cn.*pnr, tmp_d^2, true(tmp_d));
ind = (v_max==cn.*pnr);

figure('papersize', [12, 3]);
init_fig;
subplot(131);
imagesc(cn, [0, 1]);
title('local corr. image');
axis equal off tight;
subplot(132);
pnr_vmax = max(pnr(:))*0.8;
imagesc(pnr, [3, pnr_vmax]);
axis equal off tight;
title('PNR image');
subplot(133);
imagesc(cn.*pnr, [3, pnr_vmax]);
hold on;

tmp_ind = ind & (cn>=min_corr) & (pnr>=min_pnr);
[r, c] = find(tmp_ind);
ax_seeds = plot(c, r, '.m', 'markersize', 10);
axis equal off tight;
title('candidate seed pixels');
ylabel('PNR');
xlabel('Cn');tmp_d = max(1,round(gSiz/4));
v_max = ordfilt2(cn.*pnr, tmp_d^2, true(tmp_d));
ind = (v_max==cn.*pnr);

end
