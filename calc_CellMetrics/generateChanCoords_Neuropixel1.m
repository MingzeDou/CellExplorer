function chanCoords = generateChanCoords_Neuropixel1
% Channels coordinates for Neuropixels 1.0
% Matches "real" NP1.0: single shank, 960 sites, 2 columns staggered.
% Fields follow your NP2.0 helper for consistency.

chanCoords = struct();
chanCoords.source            = 'generateChanCoords_Neuropixel1';
chanCoords.layout            = 'poly2';          % simple 2D rendering
chanCoords.shankSpacing      = 0;                % single shank → 0
chanCoords.verticalSpacing   = 20;               % µm
chanCoords.horizontalSpacing = 32;               % µm

nSites = 960;           % physical sites on NP1.0 shank
nPerCol = nSites/2;     % 480 per column

% Pre-allocate
x = nan(1, nSites);
y = nan(1, nSites);
shank = ones(1, nSites);      % single shank
channels = 1:nSites;

% Column A (index 1,3,5,...)
idxA = 1:2:nSites;
x(idxA) = 0;
y(idxA) = 0:chanCoords.verticalSpacing:(nPerCol-1)*chanCoords.verticalSpacing;

% Column B (index 2,4,6,...) — staggered by +10 µm
idxB = 2:2:nSites;
x(idxB) = chanCoords.horizontalSpacing;
y(idxB) = (0:chanCoords.verticalSpacing:(nPerCol-1)*chanCoords.verticalSpacing) + chanCoords.verticalSpacing/2;

chanCoords.x        = x;
chanCoords.y        = y;
chanCoords.shank    = shank;
chanCoords.channels = channels;
end
