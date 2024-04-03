data = load('matrices\regmatD1.mat', 'D');
D1 = data.D;
D1_singular = make_singular(D1,1);