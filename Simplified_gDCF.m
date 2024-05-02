function [DCF_time] = Simplified_gDCF(k,mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Written by Oluyemi Aboyewa
% % Version 1.0:  January, 2024.
% %    - k space data normalized between [-0.5,0.5]
% %    - mat = recon matrix size 
% %    Ref: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
Ns =size(k,1);
Nproj =size(k,2);
DCF_time =zeros(size(k));
UnitDistance =1/mat;
input_X =real(k(:,:,1));
input_Y =imag(k(:,:,1));
    
[Ang, ~] = cart2pol(real(k), imag(k));
Ang = mod(Ang, 2*pi);

for Selected_Proj = 1 : size(k, 2)   %can replace with parfor
       xsel=input_X (:,Selected_Proj);
       ysel=input_Y (:,Selected_Proj);
       dist =[];
    for Selected_Meas = 1 : size(k,1) 
        indx_1 =[];
        Ref_X  = xsel (Selected_Meas);
        Ref_Y  = ysel (Selected_Meas);
        distance = -1*ones(size(input_X));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %To reduce computation time
        A= input_X >=Ref_X-UnitDistance & input_X<=Ref_X+UnitDistance; 
        B= input_Y >=Ref_Y-UnitDistance & input_Y<=Ref_Y+UnitDistance; 
        cond= A & B;
        indx_1= find(cond);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        calc_dist= (UnitDistance - sqrt( (Ref_X-input_X(indx_1)).^2 + (Ref_Y-input_Y(indx_1)).^2 ) )/UnitDistance;
        distance(indx_1) = calc_dist;
        dist(:,:,Selected_Meas) =distance;
    end
    ov = dist;
    [tmp_ind, ~] = find( ov(:) <= 0 );
    size(tmp_ind);
    ov(tmp_ind ) = 0;
    ov=reshape(ov,Ns*Nproj,[]);
    gDCF = 1./sum(ov,1);
    DCF(:,Selected_Proj)=gDCF';
end

dcf_1 =abs(DCF.*exp(1i*(Ang(:,:,2:end))));
DCF_time(:,:,1) = DCF;
DCF_time(:,:,2:end)=dcf_1;
elapsedTime = toc;
disp(['Elasped time: ' num2str(elapsedTime) 'second']);
end