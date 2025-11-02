function xs = lpfilt(x, fc, fs, order)
% LPFILT   Digital zero-phase lowpass Butterworth filtering.
%    XS = LPFILT(X, FC, FS, O) performs zero-phase filtering of the signal
%    stored in the vector X, sampled with frequency FS, using a lowpass
%    Butterworth filter of order O and cutoff frequency FC.
% 
%    XS = LPFILT(X, FC) takes FS = 100 Hz (assuming FC given in Hz) and
%    ORDER = 2.
%    
%    NaN values in X are ignored; XS has NaN values in the same places.
    
    if ~exist('fs', 'var'), fs = 100; end
    if ~exist('order', 'var'), order = 2; end
    
    xs = NaN(size(x));
    valid = ~isnan(x);
    [b, a] = butter(order, 2*fc/fs, 'low');
    xs(valid) = filtfilt(b, a, x(valid));
end