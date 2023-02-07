function saveReconResult(obj, options)
%SAVERECONRESULT saves the reconstruction results for the SartExperiment class
%
% DESCRIPTION:
%     saveReconResult automatically generates a filename and saves the
%     reconstruction results. Alternatively a custom filename can be
%     provided.
%
% USAGE:
%     sart.saveReconResult
%     sart.saveReconResult('..\UST_Imaging\SART\phantom_recon')
%
% INPUTS:
%     obj      - object, instance of the SartExperiment class
%
% OPTIONAL INPUTS:
%     absFname - filename with its absolute path specified. If the suffix
%                '.mat' is not included it will be added.
%
% ABOUT:
%     author   - Morgan Roberts
%     date     - 9th November 2022

arguments
    obj
    options.absFname = 'default'
end

switch options.absFname
    case 'default'
        stamp = char(strrep(strrep(string(datetime("now")), ':', '-'), ' ', '-'));
        fname = ['sart_', stamp, '_Nit-', num2str(obj.Nit), '_d-', num2str(round(1e3*obj.recon_d)), 'mm'];
        fname = [fname, '.mat'];
    otherwise
        fname = options.absFname;
        if strcmp(fname(end-3:end), '.mat')
            fname = fname(1:end-4);
        end
        fname = [fname, '.mat'];
end

save(fname, "obj");

end