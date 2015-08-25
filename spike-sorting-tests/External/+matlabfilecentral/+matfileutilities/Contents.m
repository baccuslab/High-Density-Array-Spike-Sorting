% MAT
%
% Files
%   AddDimension     - AddDimension adds a trailing singleton dimension to a variable in a
%   AppendColumns    - AppendColumns appends a 2D matrix to an existing 2D matrix in a MAT-file
%   AppendMatrix     - AppendMatrix appends the contents of a matrix to a variable in a MAT-file
%   AppendVector     - AppendVector adds a vector to an existing vector in a MAT-file
%   CheckIsLastEntry - CheckIsLastEntry checks that a variable is the last entry in a MAT-file
%   CreateMatrix     - CreateMatrix is presently unused
%   endian           - endian returns the endian format for the specified MAT file
%   GetLastEntry     - GetLastEntry checks that a variable is the last entry in a MAT-file
%   MATOpen          - MATOpen opens a MAT file in appropriate endian mode and returns a handle
%   RestoreDiscClass - RestoreDiscClass changes the class of the data in a MAT file (v6)
%   VarRename        - VarRename overwrites a variable name in a Level 5 Version 6 MAT-File
%   where            - returns byte offsets to the variables in a -v6 MAT-file.

% PRIVATE
%
% Files
%   AppendData             - AppendData adds data to an existing variable in a v6 MAT-file
%   ByteAlign              - ByteAlign aligns the file position to an 8 byte boundary
%   ChangeDimensions       - ChangeDimensions is called to increase the size of an existing variable

%   GetSmallDataElement    - GetSmallDataElement returns a small data element from a MAT-file
%   PadToEightByteBoundary - PadToEightByteBoundary does what its name suggests
%   sizeof                 - returns the size in bytes of the class
%   StandardMiCodes        - return Matlab standard codes for data formats
%   StandardMxCodes        - return Matlab standard codes for data formats
%   argcheck               - argcheck does error checking for the MAT-file utilities
