function B = ColdHot(m)
% A colormap for blue cold, white zero, Hot positives.

if nargin < 1, m = 256; end

n = fix(m/8);

% Create cold part:
A = [
    102 0 102;
    0 41 253;
    102 153 255;
     41 255 255;
     255 255 255]/255;
%A = ones(size(A)) - A;
 
v = [n-1 n n n];

cold = linspacev(A,v);

% Create hot part:
A = [
   255 255 255;
   255 255 0;
   255 102 41;
   255 0 0;
   102 41 0]/255;

v = [n n n n-1];
myhot = linspacev(A,v);


B = [cold; myhot];

%B = [B; flipud(hot(fix(m/2)))];


% Original cold part, 8/2/02:
A = [
    102 0 102;
     41 0 153;
     0 0 204;
     42 102 255;
     255 255 255]/255;