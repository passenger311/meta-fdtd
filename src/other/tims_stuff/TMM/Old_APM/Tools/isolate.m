
% -------------------------------------------------------------------------
% |                                                                       |
% |                           Function 'isolate'                          |
% |                                                                       |
% |  This function reqursively subdivides the original contour of         |
% |  integration using a quadtree type method evaluating each             |
% |  sub-rectangle for roots using the function 'Sub-divide'.  This       |
% |  process is repated until all roots have been isolated.               |
% |                                                                       |
% |                                                                       |
% |                   (re0,im1)---------------(re1,im1)                   |
% |                       |           |     |.    |                       |
% |   Im ^                |    .      |     |     |                       |
% |      |                |           |-----|-----|                       |
% |      |                |      .    |     |     |                       |
% |      |                |-----------|-----------|                       |
% |       ------>         |           |           |                       |
% |             Re        |           | .         |                       |
% |                       |    .      |           |                       |
% |                       |           |        .  |                       |
% |                   (re0,im0)---------------(re1,im0)                   |
% |                                                                       |
% |                                                                       |
% |  The function takes input 'roots' an array of the form:               |
% |                 (Num_R,re0,re1,im0,im1)                               |
% |                                                                       |
% |  Where Num_R is the number of roots inside the original contour.      |
% |  re0 and re1 are the min and max values of the real co-ordinates of   |
% |  the original contour.                                                |
% |  im0 and im1 are the min and max values of the imag co-ordinates of   |
% |  the original contour.                                                |
% -------------------------------------------------------------------------


function roots_mat = isolate(roots)

Num_R = roots(1);
roots_mat = [];

if Num_R == 1
    
    roots_mat = roots(2:end);
end


while max(roots(:,1)>1.5)

    hi = size(roots,1);

    for loop = 1:hi

        n  = roots(1,1);
        x0 = roots(1,2);
        x1 = roots(1,3);
        y0 = roots(1,4);
        y1 = roots(1,5);

        roots2  =   Sub_divide(x0,x1,y0,y1,n);
        ind     =   roots2(:,1)<1.1 & roots2(:,1)>0.9;
        ind2    =   roots2(:,1)>1.5;

        roots(1,:)  =   [];
        roots       =   [roots; roots2(ind2,:)];
        roots_mat   =   [roots_mat; roots2(ind,2:end)];

        clc
        txt = sprintf('%.0f out of %.0f',size(roots_mat,1),Num_R);
        fprintf('\nStage 1 Isolating roots  : %s',txt);
        fprintf('\nStage 2 Evaluating roots :');

    end

end

end



