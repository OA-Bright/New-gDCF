function [DCF_time] = gDCF_extended_2D(k,mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Written by Oluyemi Aboyewa
% % Version 1.0:  March, 2025.
% %    - k space trajectory normalized between [-0.5,0.5]
% %    - k space dims = [Nreadout Nshots Ntimeframe]
% %    - mat = recon matrix size 
% %    Ref: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
Ns = size(k, 1);
Nproj = size(k, 2);
DCF_time = zeros(size(k));

UnitDistance = 1/ mat;
input_X = real(k(:,:,1));
input_Y = imag(k(:,:,1));

[Ang, ~] = cart2pol(real(k), imag(k));
Ang = mod(Ang, 2 * pi);

DCF = zeros(Ns, Nproj);

parfor Selected_Proj = 1:Nproj % use parallel computing or change parfor to for
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
DCF_time( find( DCF_time == 0 ) ) = eps;
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);

end