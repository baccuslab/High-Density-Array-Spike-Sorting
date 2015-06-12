function Cte = Cce2Cte(Cce, nC)
% Compute the time embedded (te) representation of the noise block covariance
% matrix C with nC channels from the channel embedded (ce) representation.
% The time embedding is the one chosen e.g. in Pouzat2002 or Franke2010.
% The channel embedding is chosen more in the signal processing literature,
% e.g. in Luetkepohl2005 and is more numerically robust.

