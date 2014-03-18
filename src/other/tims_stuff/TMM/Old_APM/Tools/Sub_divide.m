
% -------------------------------------------------------------------------
% |                                                                       |
% |                           Function isolate                            |
% |                                                                       |
% |  This function subdivides a rectangluar contour into 4 subrectangles  |
% |  and then finds the number of roots in each subrectangle.  By         |
% |  recursively subdividing each rectangle until only one root remains   |
% |  in each, the function z_value can then be applied to find the exact  |
% |  position of the root.  The subdivision occurs as shown below:        |
% |                                                                       |
% |                               (rem,im1)                               |
% |                   (re0,im1)---------------(re1,im1)                   |
% |                       |           |           |                       |
% |                       |   Box 2   |   Box 3   |                       |
% |   Im ^                |           |           |                       |
% |      |                |           |           |                       |
% |      |      (re0,imm) |-----------|-----------| (re1,imm)             |
% |      |                |           |           |                       |
% |       ------>         |   Box 1   |   Box 4   |                       |
% |             Re        |           |           |                       |
% |                       |           |           |                       |
% |                   (re0,im0)---------------(re1,im0)                   |
% |                               (rem,im0)                               |
% |                                                                       |
% |  This function calls the function num_roots to calculate the number   |
% |  of roots contained in each subrectangle.  The function takes in the  |
% |  coordinates of the main rectangle along with the number of roots it  |
% |  contains.  A for loop integrates the 4 sub rectangles, however to    |
% |  improve efficiency if the number of roots found is equal to the      |
% |  expected number of roots before all subrectangles have been          |
% |  evaluated further integration is not performed. e.g if the main      |
% |  rectangle contains 10 roots and 10 roots are found in box 1 there is |
% |  no need to evaluate the other 3 boxes.  The output of this function  |
% |  has the form (number of roots, min real, max real, min imag,         |
% |  max imag).                                                           |
% -------------------------------------------------------------------------



function roots = Sub_divide(re0,re1,im0,im1,Num_R)

rem = ((re1-re0)/2)+re0;
imm = ((im1-im0)/2)+im0;

box1 = [re0, rem, im0, imm];
box2 = [re0, rem, imm, im1];
box3 = [rem, re1, imm, im1];
box4 = [rem, re1, im0, imm];

box = [box1; box2; box3; box4];
total = 0;
j=0;

for i = 1:3

    if total<Num_R

        x0 = box(i,1);
        x1 = box(i,2);
        y0 = box(i,3);
        y1 = box(i,4);

        num = num_roots([x0 x1 y0 y1]);

        if num>0.5

            j = j+1;
            total = total+num;
            roots(j,:) = [num x0 x1 y0 y1];

        end
        
    else
        break

    end

end

if total<Num_R
    x0 = box(4,1);
    x1 = box(4,2);
    y0 = box(4,3);
    y1 = box(4,4);
    num = Num_R-total;
    roots(j+1,:) = [num x0 x1 y0 y1];
end


end
