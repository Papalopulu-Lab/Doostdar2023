function [WT, dominantPeriodWT, significantWT, dominantPeriodWTSignificant, period] = significantCWT(x, dt, numRepeats)
    %This code calculates the continuous wavelet transform using the
    %in-built Matlab cwt function, and determines an average dominant period
    %by identifying the highest power period at each time point and taking
    %the mean of all the highest power periods. Additionally the code
    %calculates the wavelet transform of a randomly shuffles version on the
    %input signal, which is then used as a significance cutoff. A wavelet
    %period is determined significant if it's power is above 95% of the random signal wavelet transform power. 

    % x is an individual timeseries
    % dt is the time interval
    % numRepeats is the number of times to repeat the cwt on the shuffled
    % signal. I have found that only one repeat is sufficient, as the
    % length of the signal gives sufficient data for each period to
    % calculate a 95% confidence threshold.
    
    Fs = hours(dt/60);
    [WT, period] = cwt(x, Fs);    % Compute the continuous wavelet transform using inbuilt Matlab function -- this uses Morse wavelets with set parameters, and the data does not need detrending beforehand.
    WT = abs(WT);                 % Remove complex part of the output
    [~, maxIdx] = max(WT, [], 1); % Peak value at each time point
    dominantPeriodWT = mean(period(maxIdx)); %Average period without any signficance testing
    
    %Generate continuous wavelet transforms from randomly shuffling the
    %input timeseries x.
    WT_RAND = [];

    for n = 1:numRepeats
        randIdx = randperm(length(x));     %Random index values
        xRand = x(randIdx);                %Shuffle the timeseries using random index
        [WTrand, period] = cwt(xRand, Fs); %Calculate the wavelet transform of the randomly shuffled timeseries
        WT_RAND = [WT_RAND, WTrand];
    end

    %For each peak, determine if its above significance level (2 standard deviation should encompass 95% of the data (assumes a normal distribution)
    meanWTrand = mean(abs(WT_RAND), 2);
    stdWTrand = std(abs(WT_RAND), [], 2);
    upper95_WTrand = meanWTrand + 2*stdWTrand; 
    upper95_WTrand = repmat(upper95_WTrand, 1, size(WT,2));
    significantWT = WT;
    significantWT(WT < upper95_WTrand) = NaN; 
    
    [~, maxIdxSignificant] = nanmax(significantWT, [], 1); % Peak value at each time point
    dominantPeriodWTSignificant = mean(period(maxIdxSignificant)); % Find dominant significant period
        
end