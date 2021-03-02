function fout = coeff_polynomial
%coefficients for cubic polynomial interpolation across layers with
%different level of refinement

fout.T_4 = @(y, y1,y2,y3,y4) (y-y1).*(y-y2).*(y-y3)./( (y4-y1).*(y4-y2).*(y4-y3));
  
  
fout.T_3 = @(y, y1,y2,y3,y4) (y-y1).*(y-y2).*(y-y4)./( (y3-y1).*(y3-y2).*(y3-y4));
  
fout.T_2 = @(y, y1,y2,y3,y4) (y-y1).*(y-y4).*(y-y3)./( (y2-y1).*(y2-y4).*(y2-y3));

fout.T_1 = @(y, y1,y2,y3,y4) (y-y4).*(y-y2).*(y-y3)./( (y1-y4).*(y1-y2).*(y1-y3));

fout.prime_T1 = @(y, y1,y2,y3,y4) (y-y4).*(y-y2)./( (y1-y4).*(y1-y2).*(y1-y3)) + ...
   (y-y4).*(y-y3)./( (y1-y4).*(y1-y2).*(y1-y3))+ ...
   (y-y2).*(y-y3)./( (y1-y4).*(y1-y2).*(y1-y3));

fout.prime_T2 = @(y, y1,y2,y3,y4) (y-y4).*(y-y3)./( (y2-y1).*(y2-y4).*(y2-y3)) +...
    (y-y1).*(y-y3)./( (y2-y1).*(y2-y4).*(y2-y3)) +...
    (y-y1).*(y-y4)./( (y2-y1).*(y2-y4).*(y2-y3));

fout.prime_T3 =  @(y, y1,y2,y3,y4) (y-y2).*(y-y4)./( (y3-y1).*(y3-y2).*(y3-y4))+...
    (y-y1).*(y-y4)./( (y3-y1).*(y3-y2).*(y3-y4)) +...
    (y-y1).*(y-y2)./( (y3-y1).*(y3-y2).*(y3-y4));

fout.prime_T4 =  @(y, y1,y2,y3,y4) (y-y2).*(y-y3)./( (y4-y1).*(y4-y2).*(y4-y3)) +...
    (y-y1).*(y-y3)./( (y4-y1).*(y4-y2).*(y4-y3)) + ...
    (y-y1).*(y-y2)./( (y4-y1).*(y4-y2).*(y4-y3));
