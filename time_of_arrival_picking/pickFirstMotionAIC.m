function minAICIndex = pickFirstMotionAIC(data, options)
%PICKFIRSTMOTIONAIC Compute start of signal from minimum of the AIC.
%
% DESCRIPTION:
%     pickFirstMotionAIC identifies the first-motion-sample of a signal
%     using a version of the Akakie information criterion (AIC) described
%     in [1]. The AIC is first computed, and the index of the minimum AIC
%     is returned. The optional input IgnoreRingDown can be used to prevent
%     the AIC from late-picking-errors after the main peak. The optional
%     input NoiseLevel can be used to add gaussian noise to the signal to
%     reduce susceptibility to early-picking errors. Both of these are
%     especially useful for experimental data.
%     
% USAGE:
%     minAICInd = pickFirstMotionAIC(data)
%     minAICInd = pickFirstMotionAIC(data, Plot=true)
%
% INPUTS:
%     data            - [numeric] 1D or 2D matrix of signal data. For 2D
%                       data, the time dimension must be dim 2.
%
% OPTIONAL INPUTS:
%     Plot            - [Boolean] Option to plot the AIC and minimum index.
%                       Default = false.
%     NoiseLevel      - [numeric] Noise standard deviation, relative to the
%                       maximum signal value. Increasing this value reduces
%                       susceptibility to early-picking-errors.
%     IgnoreRingDown  - [Boolean] Whether to ignore the signal after the
%                       main peak.
%
% OUTPUTS:
%     minAICIndex     - [integer] Index of the minimum AIC for each signal.
%
% REFERENCES:
%     [1] Li, C., Huang, L., Duric, N., Zhang, H. and Rowe, C., 2009. An
%     improved automatic time-of-flight picker for medical ultrasound
%     tomography. Ultrasonics, 49(1), pp.61-72.  
%
% ABOUT:
%     original author - Bradley E. Treeby
%     modified by     - Morgan Roberts
%     date            - 25th March 2022
%     last update     - 22nd November 2022

arguments
    data {mustBeNumeric, mustBeReal, mustBeDims(data, 2)}
    options.Plot (1,1) logical = false
    options.NoiseLevel (1,1) {mustBeNumeric, mustBeReal, mustBeDims(options.NoiseLevel, 2)} = 0;
    options.IgnoreRingDown (1,1) logical = false;
end

% check data is in correct orientation if 1D, should be (sig_num, sig_ind)
if size(data, 2) == 1
    data = data.';
end

% add Gaussian noise relative to the maximum value in the signal
rng('default');
rng(1);
noise = options.NoiseLevel * max(data) * randn(size(data));
data  = noise + data;

% calculate the AIC
AIC = inf(size(data));
Nt = size(data, 2);
for k = 2:(Nt - 1)
    AIC(:, k) = k .* log(...
                          var(data(:, 1:k), 0, 2) ...
                        ) + ...
     (Nt - k - 1) .* log(...
                          var(data(:, k + 1:end), 0, 2) ...
                        );
end

% Setting the stop index to ignore the ring down
if options.IgnoreRingDown
    [~, i_stop] = max(abs(data));
else
    i_stop = Nt;
end
% AIC = AIC(1:i_stop);

% get signal arrival from minimum AIC, constrained to be less than i_stop
[~, minAICIndex] = min(AIC(1:i_stop), [], 2);

% plot
if options.Plot
    figure;
    subplot(2, 1, 1);
    plot(data, 'r', 'linewidth', 1.5);
    xline(minAICIndex);
    title('Signal');
    xlim([-Inf, length(data)]);
    xlabel('Sample Index');

    subplot(2, 1, 2);
    hold on;
    plot(AIC, 'b');
    plot(AIC(1:i_stop), 'k-', 'linewidth', 1.5);
    h = xline(minAICIndex, 'k--');
    title('AIC');
    xlim([-Inf, length(data)]);
    xlabel('Sample Index');
    legend(h, {'Minimum AIC'});
end

end

function mustBeDims(input,numDims)
    % Test for number of dimensions    
    if ~isequal(length(size(input)),numDims)
        eid = 'Size:wrongDimensions';
        msg = ['Input must have ',num2str(numDims),' dimension(s).'];
        throwAsCaller(MException(eid,msg))
    end
end
