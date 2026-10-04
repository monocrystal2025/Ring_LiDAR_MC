function result = computeBackground(cfg)
%COMPUTEBACKGROUND Pd(wD) under background changes and receiver controls.
arguments
    cfg (1,1) struct
end
bank = discussion.bank(cfg.N,cfg.seed);
backgroundScales = [1,10,100,1e4];
receiverModes = ["matched","envelope","common"];
curves = cell(numel(backgroundScales),numel(receiverModes));
minima = cell(size(curves));
rows = cell(numel(curves),1);
for a = 1:numel(receiverModes)
    for b = 1:numel(backgroundScales)
        caseCfg = cfg;
        caseCfg.receiverMode = receiverModes(a);
        caseCfg.system.skyRadiance_W_m2_sr_nm = ...
            cfg.system.skyRadiance_W_m2_sr_nm*backgroundScales(b);
        if receiverModes(a) == "common"
            % A single fixed circular receiver for the entire width sweep.
            caseCfg.commonFov_rad = max(cfg.widthGrid_mrad)*1e-3/2+3*cfg.system.beamWidth_rad;
        end
        fprintf('tmp_4 receiver %s, sky x%g\n',receiverModes(a),backgroundScales(b));
        curves{b,a} = discussion2.runCurve(caseCfg,bank,"dynamic");
        minima{b,a} = discussion2.minima(curves{b,a},cfg.targetProbability);
        iso = minima{b,a};
        rows{b+(a-1)*numel(backgroundScales)} = struct('Receiver',receiverModes(a), ...
            'SkyScale',backgroundScales(b),'SkyRadiance',caseCfg.system.skyRadiance_W_m2_sr_nm, ...
            'LineWidth_mrad',iso.widths_mrad(1),'RingWidth_mrad',iso.widths_mrad(2), ...
            'Kappa',iso.kappa,'LineStatus',iso.status(1),'RingStatus',iso.status(2), ...
            'LinePeakP',max(curves{b,a}.probability(:,1)), ...
            'RingPeakP',max(curves{b,a}.probability(:,2)),'N',cfg.N);
    end
end
summary = struct2table(vertcat(rows{:}));
result = struct('config',cfg,'backgroundScales',backgroundScales, ...
    'receiverModes',receiverModes,'curves',{curves},'minima',{minima},'summary',summary);
save(fullfile(cfg.outputRoot,"tmp_4.mat"),'result','-v7.3');
writetable(summary,fullfile(cfg.outputRoot,"tmp_4_summary.csv"));
end
