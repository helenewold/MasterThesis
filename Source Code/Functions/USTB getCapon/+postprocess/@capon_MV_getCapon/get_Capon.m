%GETCAPON Implementation of the Minimum Variance / Capon filter for linear
%         arrays. Includes FB averaging, time averaging, and the use of
%         subspaces. 
%
% [imAmplitude imPower] = getCapon(dataCube, indsI, indsJ, regCoef, L, nTimeAverage, V, doForwardBackwardAveraging, verbose)
%
% dataCube     : A range x beams x array-elements image cube
% indsI, indsJ : Indices of pixels that are to be processed; indsI==0 / indsJ==0 means all ranges / beams
% regCoef      : Regularization coefficient for diagonal loading
% L            : Length of subarrays
% nTimeAverage : Includes +- this number of time lags to produce 'R' matrix
% V            : If V is an integer, we apply a V-dim beamspace, if V is a matrix, it is used as a projection matrix, if V is empty ([]) than the full space is used (no subspacing)
% doForwardBackwardAveraging : Whether to do forward-backward averaging
% usePinv      : {Default: false} Whether or not to use the 'pinv' command when inverting the matrix
%
% Note I: If using a subspace matrix, V, enabling FB averaging gives a
% substantial increase in the computational load.  This is because the
% reduction of the dimensionality must happen at a much later stage in the
% code. 
%
% Note II: This version assumes a=ones(L,1) in w = Ri*a/(a'*Ri*a), i.e., we
% search for the amplitude/power in the direction perpedicular to the
% linear array.
%
% Note III: Matrix inversion is done like this: Ri = inv(R + regCoef/L*trace(R)*I); 
% or using pinv if so chosen.
%
% Note IV: Whether or not it as faster to do beamspace projection before or
% after the calculation of R is disputed -- we should perhaps just settle
% for post-R-creation projection for simplicity.
%
% Last modified:
% 2009.08.25 - A.C. Jensen {Created the function}
% 2009.08.27 - A.C. Jensen {Robustified indsI in light of nTimeAverage use}
% 2009.09.09 - A.C. Jensen {By popular request, added the factor 1/L in the diagonal loading}
% 2011.06.24 - A.C. Jensen {Added the 'usePinv' parameter}
% 2011.06.29 - A.C. Jensen {Reverse index order when creating the 'd' matrix (reshuffling data) to gain speed}
% 2011.12.02 - A.C. Jensen {Increased user friendliness: If V is an integer, a V-dimensional simple beamspace is applied}
% 2012.06.09 - J.P. �sen {Matrix-size allocation bug fixed for computational gain}

% 2022.08.26 - H. Wold {Updated input syntax, updated verbose}

function [imAmplitude, imPower] = get_Capon(h, dataCube, apod_matrix, progress)

% arguments
%     dataCube                            (:,:,:)                 % Data on the shape; Range, beams, array elements
%     options.indsI                       (1,1) = 0               % Pixels to be processed. 0 means all ranges
%     options.indsJ                       (1,1) = 0               % Pixels to be processed. 0 means all beams
%     options.regCoef                     (1,1) = 1/100           % Diagonal loading regularization coefficient
%     options.nTimeAverage                (1,1) = 10              % Number of time lags to produce 'R' matrix
%     options.L                           (1,1) = 2               % Length of subarray
%     options.V                           (:,:) = []              % If V is an integer, we apply a V-dim beamspace, if V is a matrix, it is used as a projection matrix, if V is empty ([]) than the full space is used (no subspacing)
%     options.doForwardBackwardAveraging  (1,1) = 0               % Whether or not to do forward-backward averaging
%     options.verbose                     (1,1) = 1               % Whether or not to show progress bar
%     options.usePinv                           = false           % Whether or not to use the 'pinv' command when inverting the matrix
% %     options.f_number                    (1,1)                   % F_number added
%     options.apod_matrix                 (:,:,:)                  % 
% end

doForwardBackwardAveraging  =   h.doForwardBackward;
indsI       =   h.indsI;      indsJ           =   h.indsJ;
regCoef     =   h.regCoef;    n    =   h.K_samples;
L           =   round(h.L_elements);
V               =   h.V; 
% verbose     =   h.verbose;    
usePinv         =   h.usePinv;
%     apod_matrix    =   h.apod_matrix;


[N, E, M] = size(dataCube);
L_frac = L/M;

if indsI<=0, indsI = 1:N; end
if indsJ<=0, indsJ = 1:E; end


% Skip pixels that cannot be calculated due to time averaging
indsI = indsI( indsI > n );
indsI = indsI( indsI < (N-n+1) );


% If V is an integer, create a V-dimensional beamspace and place it in V
if size(V,1)==1 && size(V,2)==1
  nSubDimensions = V;
  V = zeros(L,nSubDimensions);
  for m=0:nSubDimensions-1
    V(:,m+1) = exp(-complex(0,1)*((0:L-1)-(L-1)/2)*(0-(m-floor(nSubDimensions/2))*2*pi/L));
  end
  V = V*diag(1./diag(V'*V)).^.5;  % Normalizes the columns of V
end


% a = ones(L,1);
% n = n;
imAmplitude = zeros(N,E);
imPower = zeros(N,E);
if h.calc_param
    CN = zeros(N,E);
    CN_DL = zeros(N,E);
    Tse = zeros(N,E);
    Reigval = zeros(N,E);
    Reigval_DL = zeros(N,E);
    CosineSimilarity = zeros(N,E);
end
h.StructTmp = struct('R',cell(N,E));
h.subarrsize = zeros(N,E);

useSubspace = ~isempty(V);
if useSubspace
    nSubspaceDims = length(V(1,:));
    fprintf('Capon algorithm "subspaced" down to %d dims.\n', nSubspaceDims);
end


tools.workbar();
for i=indsI
    tools.workbar(i/length(indsI), sprintf('%s %s (%s)',h.name,h.version,progress),'getCapon');

    for j=indsJ
        % If apodization matrix is given
        rx_apod = squeeze(apod_matrix(i, j, :));
        idx = find( abs(rx_apod) > h.active_element_criterium);
        M_new = length(idx);

        if h.L_elements_set
            % Kan tas ut dersom L kan være større enn Mnew
            if L<M_new
                L_new = L;
            else
                L_new = M_new-1; % Må gjøres (-1) pga dimension error
            end
        else
        % Alt som normalt
            L_new = floor(L_frac*M_new);
        end

        I_new = eye(L_new);
        a = ones(L_new, 1);
        
        if j == 64 && i == 265
            L_new
        end


        % OBS: Må omdefinere ifht. hva a og L_new skal være
        if useSubspace
          I_new = eye( nSubspaceDims );
          a = V'*a;  % The column of ones (what we seek the 'magnitude' of) represented in the subspace
        end
    
        ar = squeeze(dataCube(i-n:i+n, j, idx));  % Array responses plus-minus n time steps

        if n==0, ar = transpose(ar); end
    
        % Place the array responses in one (K-L+1)*(2n+1) x L matrix:
        ar = transpose(ar);
        d = zeros(L_new, (M_new-L_new+1)*(2*n+1));
        for ii=1:M_new-L_new+1
            d(:, ((ii-1)*(2*n+1)+1):((ii)*(2*n+1))) = ar(ii:ii+L_new-1, :);
        end
        d = d';
        
        % If a subspace matrix V is given and we're _not_ using FB averaging, we
        % can use V to reduce the dimensionality of the data now:
        if (useSubspace) && (~doForwardBackwardAveraging)
            d = d*V;
        end
        
        R = d'*d / ((M_new-L_new+1) * (2*n+1));   % R contains an estimate of the covariance matrix
    
        % Store the sum of the current-time outputs in 'g_singleSnapshot':
        g_singleSnapshot = sum(d(n+1:(2*n+1):end, :))' / (M_new-L_new+1);
    
        if doForwardBackwardAveraging
            R = 0.5*(R + rot90(conj(R), 2));
        end
    
        % If a subspace matrix V is given and we _are_ using FB averaging, we
        % "have" to wait until now to go to the reduced space:
        if (useSubspace) && (doForwardBackwardAveraging)
            R = V'*R*V;
            g_singleSnapshot = V'*g_singleSnapshot;
        end
        
        % Check is the matrix is invertible. If not, ignore inverting it.
        if R == 0
            Ri = R + regCoef/L_new*trace(R)*I_new;
        else
            if usePinv
                Ri = pinv(R + regCoef/L_new*trace(R)*I_new);
            else
                Ri = inv(R + regCoef/L_new*trace(R)*I_new);
            end

        end

        w = Ri*a / (a'*Ri*a);


        if h.USTB_normalization
            imAmplitude(i,j) = w'*g_singleSnapshot*length(idx);   % Note: A bit ad-hoc maybe, but uses only the current time-snapshot to calculate the output/'alpha' value
        else
            imAmplitude(i,j) = w'*g_singleSnapshot;
        end
        imPower(i,j) = abs(w'*R*w);


        if h.calc_param
            [CN(i,j), CN_DL(i,j), Tse(i,j), Reigval(i,j), Reigval_DL(i,j)] = parameter_calculations(R, R + regCoef/L_new*trace(R)*I_new, w, M_new);
        end

        h.StructTmp(i, j).R(:,:) = R(:,:);
        h.StructTmp(i, j).R_DL(:,:) = R(:,:) + regCoef/L_new*trace(R(:,:))*I_new;
        h.StructTmp(i,j).w(:) = w;
        h.StructTmp(i,j).M = M_new;
        h.StructTmp(i,j).R_eigvals(:,:) = eig(R);
        h.StructTmp(i,j).R_eigvals_DL(:,:) = eig(R(:,:) + regCoef/L_new*trace(R(:,:))*I_new);
        h.subarrsize(i,j) = L_new;
    end
end
% hold off
tools.workbar(1);
if h.calc_param
    h.CN        =   CN(:);
    h.CN_DL  = CN_DL(:);
    h.Tse = Tse(:);
    h.Reigval = Reigval(:);
    h.Reigval_DL = Reigval_DL(:);
%     h.CosineSimilarity = CosineSimilarity(:);
end
