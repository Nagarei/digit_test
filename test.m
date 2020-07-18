function test()

dim = 5;
testcase = 1000;
r1 = [zeros(1, testcase) ;random('Normal',0,1, [dim-1 testcase])];
r2 = [random('Normal',10,2, [dim-1 testcase]); zeros(1, testcase)];

[w, mu, S] = EM(2,[r1 r2 r1])

end
