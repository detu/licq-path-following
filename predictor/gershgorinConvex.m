function [Weights] = gershgorinConvex(M,rho)
%GERSHGORINCONVEX Summary of this function goes here
% 
% Convexification using Gershgorin Circle criteria.
%
% [OUTPUTARGS] = GERSHGORINCONVEX(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/04/18 07:24:06 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

% read dimension of M
numM    = size(M,1);
Weights = zeros(numM,numM);

% read nonzero elements from input matrix
[row,column,value] = find(M);

%nnzData = [row column value];

% iterate in the nonzero indeces and employ Gershgorin circle criteria.
%rho = 1;   % added value

% think how to read row, colum, and value !
% because row may have duplicate entries ! read http://se.mathworks.com/help/matlab/math/accessing-sparse-matrices.html?refresh=true

initRow         = row(1);
searchOccurance = find(row==initRow);
numberOccurance = numel(searchOccurance);

% for i=1:row
%     sumRow = 0;
%     diag   = 0;
%     % iterate in the column
%     for j=1:column
%         if (i==j)
%             diag = M(i,j);
%         else
%             sumRow = sumRow + abs(value);
%         end
%     end
%     % check estimated eigenvalue
% end
sum = 0;
while (numberOccurance > 0)
    for i=1:numberOccurance
        if (row(searchOccurance(i)) == column(searchOccurance(i)))
            diag = value(searchOccurance(i));
        else
            sum = sum + abs(value(searchOccurance(i)));
        end
    end
    
    if ((diag - sum ) < 0)
        %Weights(initRow,initRow) = sum - diag + rho;
        Weights(initRow,initRow) = sum + rho;
    end
    
    % delete already processed data
    row(row == initRow) = [];
    
    if isempty(row)
        break;
    else
        % update new data
        initRow = row(1);
        searchOccurance = find(row==initRow);
        numberOccurance = numel(searchOccurance);
        % reset counter
        sum  = 0;
        diag = 0;
    end
    
    
end


end
