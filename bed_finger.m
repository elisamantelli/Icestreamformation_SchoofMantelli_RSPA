function fout = bed_finger(v_in, parameters)
%bed_finger.m returns bed elevatiom as a function of the downstream coordinate x

B=@(x) -(parameters. bed. b0+parameters. bed. b1.*x);
fout = B(v_in);
