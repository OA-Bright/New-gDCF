function [DCF_time] = gDCF_new(k,mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Written by Oluyemi Aboyewa
% % Version 2.0:  October, 2024.
% %    - k = complex kspace trajectory(kx +iKy) normalized between [-0.5,0.5]
% %    - mat = recon matrix size 
% %    Ref:ISMRM abstract 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
Ns = size(k, 1);
Nproj = size(k, 2);
DCF_time = zeros(size(k));
%% lower Beruling Density Calculation
% uncomment for symetric readout about the center e.g radial projection acquisition
% lengths = linspace(0, mat, (mat/2) + 1); 
% widths = linspace(0, mat, (mat/2) + 1);

% use for center out readout e.g spiral as in demo .mat file
lengths = linspace(0, mat, mat); 
widths = linspace(0, mat, mat); 
[~, del_k] = LBD(k * mat, lengths, widths);
%%
UnitDistance = del_k / mat;
input_X = real(k(:,:,1));
input_Y = imag(k(:,:,1));

[Ang, ~] = cart2pol(real(k), imag(k));
Ang = mod(Ang, 2 * pi);

DCF = zeros(Ns, Nproj);

parfor Selected_Proj = 1:Nproj
    xsel = input_X(:, Selected_Proj);
    ysel = input_Y(:, Selected_Proj);

    gDCF = zeros(Ns, 1);

    for Selected_Meas = 1:Ns
        Ref_X = xsel(Selected_Meas);
        Ref_Y = ysel(Selected_Meas);

        % Logical conditions to reduce computation range
        cond = input_X >= Ref_X - UnitDistance & input_X <= Ref_X + UnitDistance & ...
               input_Y >= Ref_Y - UnitDistance & input_Y <= Ref_Y + UnitDistance;
        indx_1 = find(cond);

        % Calculate distance only for valid points
        if ~isempty(indx_1)
            calc_dist = (UnitDistance - sqrt((Ref_X - input_X(indx_1)).^2 + (Ref_Y - input_Y(indx_1)).^2)) / UnitDistance;

            % Only keep positive distances
            calc_dist(calc_dist < 0) = 0;
            % Sum the contributions for this point
            gDCF(Selected_Meas) = 1 / sum(calc_dist); 
        else
            gDCF(Selected_Meas) = 0;  % Handle the case where no valid points were found
        end
    end

    % Store gDCF result for this projection
    DCF(:, Selected_Proj) = gDCF;
end

% Apply phase correction for DCF along time dim
dcf_1 = abs(DCF .* exp(1i * Ang(:,:,2:end)));
DCF_time(:,:,1) = DCF;
DCF_time(:,:,2:end) = dcf_1;

elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);

end

function [D_B, delt_k] = LBD(sample_points, lengths, widths)
    % sample_points: k-space trajectory
    % lengths: Vector of length values corresponding to each rectangle
    % widths: Vector of width values corresponding to each rectangle
    
    % Ensure lengths and widths have the same number of elements
    sample_x =real(sample_points(:,:,1));
    sample_y =imag(sample_points(:,:,1));
    assert(length(lengths) == length(widths), 'Lengths and widths must have the same number of elements.');
    
    num_rectangles = length(lengths);
    densities = zeros(num_rectangles-1, 1);
    np= zeros(num_rectangles-1, 1);
    na = zeros(num_rectangles-1, 1);
    
    for i = 1:num_rectangles-1
        L = lengths(i+1);
        W = widths(i+1);
        
        % Count points in a rectangle centered at the origin
        indices  = (sample_x >= -L/2 & sample_x <= L/2) & (sample_y >= -W/2 & sample_y <= W/2);
        N_LW =sum(indices(:));
        
        % Area of the rectangle
        A_LW = L*W;
        np(i) = N_LW;
        na(i) =A_LW;
        
        % Calculate density for this rectangle
        densities(i) = N_LW / A_LW;
    end
    % Return the lower Beurling density (limit inferior)
    non_zeroDB =densities(densities>=1);
    D_B = min(non_zeroDB(:));
    delt_k = 1/sqrt(D_B);
end
