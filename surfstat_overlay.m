function surfstat_overlay(U1,U2,mask,avsurf)

    % SurfStat Overlay
    s1(mask==1) = U1;
    s1(abs(s1)<1.96) = 0;
    s2(mask==1) = U2;
    s2(abs(s2)<1.96) = 0;

    X = zeros(size(mask));
    cmap = [0.8 0.8 0.8; cool(3); autumn(3); 0.25 0.25 0.25];

    % Create partial colorbar
    cmapcolours = zeros(8,1);
    for i = 1:size(X,2)
        if (s1(i) > 0 && s2(i) < 0) || (s1(i) < 0 && s2(i) > 0)
            X(i) = 8;
            cmapcolours(8,1) = 1; % dark grey
        elseif s1(i) == 0 && s2(i) == 0
            X(i) = 1;
            cmapcolours(1,1) = 1; % grey
         elseif s2(i) == 0 && s1(i) < 0
            X(i) = 2;
            cmapcolours(2,1) = 1; % light blue
        elseif s1(i) == 0 && s2(i) < 0
            X(i) = 3;
            cmapcolours(3,1) = 1; % purple
        elseif s2(i) == 0 && s1(i) > 0
            X(i) = 5;
            cmapcolours(5,1) = 1; % red
        elseif s1(i) == 0 && s2(i) > 0
            X(i) = 6;
            cmapcolours(6,1) = 1; % orange
        elseif s1(i) < 0
            [~, ind] = max(abs([s1(i) s2(i)]), [], 1, 'linear');
            X(i) = 4;
            cmapcolours(4,1) = 1; % pink
        elseif s1(i) > 0
            [~, ind] = max(abs([s1(i) s2(i)]), [], 1, 'linear');
            X(i) = 7;
            cmapcolours(7,1) = 1; % yellow
        end
    end
    figure; [a, cb] = SurfStatView(X, avsurf);
    SurfStatColormap(cmap)
    SurfStatColLim([1,8])
end