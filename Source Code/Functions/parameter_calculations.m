%% parameter_calculations.m
% Function to calculate the values condition number and eigenvalue
% percentages both before and after diagonal loading. Calculates Tse before
% diagonal loading as well.
% Input:
%       R       -       Covariance matrix for a specific pixel
%       R_DL    -       Covariance matrix after diagonal loading
%       w       -       weights calculated from the getCapon function (can
%                       be other functions)
%       M       -       Length of array in specific pixel. Changes if
%                       there is an expanding aperture.
% Output:
%       CN, CN_DL   -   Condition number of matrix before and after
%                       diagonal loading has been applied.
%       Tse         -   Tse values
%       Reigval     -   Percentage of eigenvalues below a certain tolerance
%       Reigval_DL  -   Percentage of eigenvalues below a certain tolerance
%                       after diagonal loading har been applied
%       R_eigval    -   (Will be removed(?)) Eigenvalues of R

function [CN, CN_DL, Tse, Reigval, Reigval_DL] = parameter_calculations(R, R_DL, w, M)
    arguments
        R
        R_DL
        w
        M
    end

    tol = 1e-10;

    % Condition Number Calculation
    CN = cond(R);   % Generate condition number
    CN_DL = cond(R_DL);  % Generate condition number after diagonal loading

    % T_se Calculation
    Tse = M * sum(abs(w).^2);%sum(abs(w).^2);
    if isempty(Tse); Tse = 0; end


    % Eigenvalue Percentage Calculation
    R_eigval = eig(R);
    R_eigval_DL = eig(R_DL);


    % Calculate percentage of eigenvalues below a limit (tol)
    % Both before and after diagonal loading has been applied
    if isempty(R_eigval)
        Reigval = 0;
    else
        len = length(R_eigval);
        ind = find(abs(R_eigval)<tol);
        count = length(ind);
        Reigval = count/len*100;
    end
    if isempty(R_eigval_DL)
        Reigval_DL = 0;
    else
        len = length(R_eigval_DL);
        ind = find(abs(R_eigval_DL)<tol);
        count = length(ind);
        Reigval_DL = count/len*100;
    end

end