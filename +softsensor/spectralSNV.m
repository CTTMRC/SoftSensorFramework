function transformed = spectralSNV(spectra)
% <header>  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
% SPECTRALSNV Apply standard normal variate preprocessing to spectrum rows.
%--------------------------------------------------------------------------
% Syntax
%   transformed = softsensor.spectralSNV(spectra)
%
% Description
%   Subtract each spectrum's row mean and divide by its sample standard
%   deviation. This transformation estimates no shared cross-row statistics.
%
% Inputs
%   spectra - Real double N-by-P matrix, with observations in rows and wavelengths in columns.
%
% Outputs
%   transformed - N-by-P matrix of row-normalized spectra.
%
% Notes
%   - Complete constant spectra become zero; rows containing NaN/Inf become NaN.
%   - The denominator is std(spectra, 0, 2), using sample normalization.
%   - No file writes or changes to RNG, path or working folder.
%
% Example
%   transformed = softsensor.spectralSNV([1, 2, 3; 4, 4, 4]);
%
% Dependencies
%   Base MATLAB; mean, std, isfinite.
%
% Metadata
%   Project:               SoftSensorFrameworkPP3
%   Documentation updated: 2026-10-07
%   Status:                Custom preprocessing; artificial-data checks passed
%
% See also customPLS, generateSyntheticData.
%--------------------------------------------------------------------------
% </header> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
    arguments
        spectra (:,:) double {mustBeReal}
    end
    transformed = NaN(size(spectra));
    finite = all(isfinite(spectra),2);
    complete = spectra(finite,:);
    scale = std(complete,0,2);
    scale(scale==0) = 1;
    transformed(finite,:) = (complete-mean(complete,2))./scale;
end
