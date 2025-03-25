function [DCF_time] = gDCF_extended_3D(k,mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Written by Oluyemi Aboyewa
% % Version 1.0:  March, 2025.
% %    [Nsample,Nshot,Ndim,Nt] =size(k) %Ndim = 3 for x,y,z,Nt= timeframe
% %    - k space trajectory normalized between [-0.5,0.5,0.5]
% %    - mat = recon matrix size 
% %    Ref: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
Ns = size(k, 1);
Nproj = size(k, 2);
Nt = size(k,4);

for ii=1:Nt
    fprintf( 'TimeFrame = %3d\n', ii );
    UnitDistance = 1/mat;
    UnitArea =UnitDistance.^2;
    input_X = k(:,:,1,ii);
    input_Y = k(:,:,2,ii);
    input_Z = k(:,:,3,ii);

    DCF = zeros(Ns, Nproj);

    parfor Selected_Proj = 1:Nproj
        xsel = input_X(:, Selected_Proj);
        ysel = input_Y(:, Selected_Proj);
        zsel = input_Z(:, Selected_Proj);

        gDCF = zeros(Ns, 1);

        for Selected_Meas = 1:Ns
            Ref_X = xsel(Selected_Meas);
            Ref_Y = ysel(Selected_Meas);
            Ref_Z = zsel(Selected_Meas);

            % Logical conditions to reduce computation range
            tmp_a = (input_X >= (Ref_X-UnitDistance)) & (input_X <= (Ref_X+UnitDistance));
            tmp_b = (input_Y >= (Ref_Y-UnitDistance)) & (input_Y <= (Ref_Y+UnitDistance));
            tmp_c = (input_Z >= (Ref_Z-UnitDistance)) & (input_Z <= (Ref_Z+UnitDistance));
            tmp_abc = tmp_a & tmp_b & tmp_c;
            indx_1 = find( tmp_abc > 0 );

            % Calculate distance only for valid points
            if ~isempty(indx_1)
                calc_dist = (UnitArea - ((Ref_X-input_X(indx_1)).^2 + (Ref_Y-input_Y(indx_1)).^2 + (Ref_Z-input_Z(indx_1)).^2 ))/UnitArea;

                % Only keep positive distances
                calc_dist(calc_dist < 0) = 0;
                % Sum the contributions for this point
                gDCF(Selected_Meas) = 1./ sum(calc_dist); 
            else
                gDCF(Selected_Meas) = 0;  % Handle the case where no valid points were found
            end
        end

        % Store gDCF result for this projection
        DCF(:, Selected_Proj) = gDCF;
    end
  DCF_time(:,:,ii) = DCF;  
end
DCF_time( find( DCF_time == 0 ) ) = eps;
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);

end
