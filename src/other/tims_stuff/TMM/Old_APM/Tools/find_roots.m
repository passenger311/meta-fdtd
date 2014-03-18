function modes = find_roots(contour,nrs)
% ----------------------- Find number of roots ---------------------------

Num_R = num_roots(contour);

if (Num_R<0.5)
    clc
    error = sprintf('\nNo roots found\n\n');
    fprintf(error);
    modes = [];
elseif isnan(Num_R)==1
    clc
    error = sprintf('\nToo many roots try a smaller contour\n\n');
    fprintf(error);
    modes = [];
else
    clc
    fprintf('\nNumber of roots inside contour = %.2f',Num_R);
    pause(1);
    % --------------------------- Isolate roots ------------------------------
    tic
    roots = [Num_R contour];
    roots_mat = isolate(roots);
    size(roots_mat)
    toc
   
    % -------------------------- Evaluate roots ------------------------------

    modes1 = evaluate(roots_mat,Num_R,nrs);

    % Remove duplicate roots

    z_mat = imag(modes1(:,1));
    z_mat = abs(z_mat/pi);
    ind = z_mat >1;
    modes1(ind,:)=[];
    modes = modes1;

end



