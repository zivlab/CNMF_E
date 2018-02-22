function view_pnr_corr(obj)
%% view PNR, correations that emerge from parameters used in CNMF-E
gSig = obj.options.gSig;
gSiz = obj.options.gSiz;
min_corr = obj.options.min_corr;
min_pnr = obj.options.min_pnr;

fprintf(['The current Values: \n (gSiz, gSig) = ', num2str([gSiz, gSig]), ...
    '\n(min_corr, min_pnr) = ', num2str([min_corr, min_pnr]), '\n']);

%% Print the PNR and correlations for the current parametres

% Calculate the PNR and correlation map
fprintf('computing the correlation image and the peak-to-noise ratio image....\n');
[cn, pnr] = obj.correlation_pnr_parallel([1, 5000]);
fprintf('Done\n');

tmp_d = max(1,round(gSiz/4));
v_max = ordfilt2(cn.*pnr, tmp_d^2, true(tmp_d));
ind = (v_max==cn.*pnr);

% Plot the PNR and correlation map
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

% Plot the estimated seeds
tmp_ind = ind & (cn>=min_corr) & (pnr>=min_pnr);
[r, c] = find(tmp_ind);
ax_seeds = plot(c, r, '.m', 'markersize', 10);
axis equal off tight;
title('candidate seed pixels');
ylabel('PNR');
xlabel('Cn');tmp_d = max(1,round(gSiz/4));
v_max = ordfilt2(cn.*pnr, tmp_d^2, true(tmp_d));
ind = (v_max==cn.*pnr);
suptitle(['(gSiz, gSig) = ', num2str([gSiz, gSig]), ...
                 ' (min corr, min pnr) = ', num2str([min_corr, min_pnr])])

%% Insert changes and plot again

temp = input('Do you want to make a change? (y/n)    ', 's');
while strcmpi(temp, 'y')
    % Choose values for (gSiz, gSig)
    temp1 = input('Do you want to change (gSiz, gSig)? (y/n)    ', 's');
    while strcmpi(temp1, 'y')
        temp = input('type new values for (gSiz, gSig):  ', 's');
        try
            new_values = str2num(temp); %#ok<*ST2NM>
            gSiz = new_values(1);
            gSig = new_values(2);
            obj.options.gSig = gSig;
            obj.options.gSiz = gSiz;
            fprintf('Your current selection of parameters (gSiz, gSig) are (%.1f, %.1f).\n', gSiz, gSig);
            % Calculate the PNR and correlation map
            fprintf('computing the correlation image and the peak-to-noise ratio image....\n');
            [cn, pnr] = obj.correlation_pnr_parallel([1, 5000]);
            
            tmp_d = max(1,round(gSiz/4));
            v_max = ordfilt2(cn.*pnr, tmp_d^2, true(tmp_d));
            ind = (v_max==cn.*pnr);
            
            % Plot the PNR and correlation map
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
            
            % Plot the estimated seeds
            tmp_ind = ind & (cn>=min_corr) & (pnr>=min_pnr);
            [r, c] = find(tmp_ind);
            ax_seeds = plot(c, r, '.m', 'markersize', 10);
            axis equal off tight;
            title('candidate seed pixels');
            ylabel('PNR');
            xlabel('Cn');tmp_d = max(1,round(gSiz/4));
            v_max = ordfilt2(cn.*pnr, tmp_d^2, true(tmp_d));
            ind = (v_max==cn.*pnr);
            
            suptitle(['(gSiz, gSig) = ', num2str([gSiz, gSig]), ...
                 ' (min corr, min pnr) = ', num2str([min_corr, min_pnr])])
            
            temp1 = input('Do you want to change (gSiz, gSig)? (y/n)    ', 's');
        catch
            warning('values are bad. try again');
            continue;
        end
    end
    
    % Change values for min_corr, min_pnr (for changing the estimated
    % seeds)
    temp1 = input('Do you want to set new values for (min_corr, min_pnr)? (y/n)', 's');
    if strcmpi(temp1, 'y')
        while true
            temp1 = input('type new values for (min_corr, min_pnr):  ', 's');
            try
                new_values = str2num(temp1); %#ok<*ST2NM>
                min_corr = new_values(1);
                min_pnr = new_values(2);
                obj.options.min_corr = min_corr;
                obj.options.min_pnr = min_pnr;
                
                % Plot the PNR and correlation map
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
                
                % Plot the estimated seeds
                tmp_ind = ind & (cn>=min_corr) & (pnr>=min_pnr);
                [r, c] = find(tmp_ind);
                ax_seeds = plot(c, r, '.m', 'markersize', 10);
                axis equal off tight;
                title('candidate seed pixels');
                ylabel('PNR');
                xlabel('Cn');tmp_d = max(1,round(gSiz/4));
                v_max = ordfilt2(cn.*pnr, tmp_d^2, true(tmp_d));
                ind = (v_max==cn.*pnr);
                suptitle(['(gSiz, gSig) = ', num2str([gSiz, gSig]), ...
                 ' (min corr, min pnr) = ', num2str([min_corr, min_pnr])])
                
                temp1 = input('Do you want to set new values for (min_corr, min_pnr)? (y/n)   ', 's');
                if strcmpi(temp1, 'n')
                    fprintf('Your current selection of parameters (min_corr, min_pnr) are (%2d, %2d).\n', min_corr, min_pnr);
                    break;
                    
                end
            catch
                warning('values are bad. try again');
                continue;
            end
        end
    end
    
    temp = input('Do you want to make more changes? (y/n)    ', 's');
end

end
