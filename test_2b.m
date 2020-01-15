
h = [1/128, 1/64, 1/32, 1/16, 1/8, 1/4];
L2err = [2.5100e-05, 1.0040e-04, 4.0155e-04, 0.0016, 0.0064, 0.0255];

figure;
loglog(h,L2err);
title('log log plot of L2error against h');
xlabel('h');
ylabel('L2error');

log_h = log(h);
log_L2err = log(L2err);
slope = (log_L2err(1,6) - log_L2err(1,1))/(log_h(1,6) - log_h(1,1));
disp(slope);        % 1.9977