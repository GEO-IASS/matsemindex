clear all



%% Without keying
Q = fillwithindex( zeros([2,3,4]));
k = semindex(1,1);
i = semindex(2,3);
j = semindex(3,[1,2]);
j1 = semindex(3,4);

a = semmat(Q);
b = a(i,j) % 2x3x4 -> 2x1x2

c = a([i,j],k) % 2x3x4 -> 1x1x2 

d = a(j) % 2x3x4 -> 2x3x2

e = a(i,j1,k);


%% Without keying
Q = fillwithindex( zeros([2,3,4]));
a = semmat(Q,[20,30,40]);
k = semindex(20,1);
i = semindex(30,3);
j = semindex(40,[1,2]);
j1 = semindex(40,4);


b = a(i,j) % 2x3x4 -> 2x1x2

c = a([i,j],k) % 2x3x4 -> 1x1x2 

d = a(j) % 2x3x4 -> 2x3x2

e = a(i,j1,k) % scalar

%% Then some examples of the discrete manip - conditioning
Q = fillwithindex( zeros([2,3,4]));
a = semmat(Q,[20,30,40]);
i = semindex(30,2);
ac = a(i); % we conditioned a

%% this does automatically the memberwise
Q1 = fillwithindex( zeros([3,4]));
Q2 = fillwithindex( zeros([5,4]))*2;
a1 = semmat(Q1,[10,40]);
a2 = semmat(Q2,[20,40]);

% [3,4]*[4,5] = [3,5]
ac = binaryop(a1,a2,[],[],@(x,y) sum(reshape(x.*y,1,[])));

acX = a1.X*a2.X';

ac.X-acX

%% this does broadcast inversion
Q1 = fillwithindex( zeros([3,4,4]));
a1 = semmat(Q1,[10,20,30]);
ai = unaryop(a1,[20,30],[20,30],sizeofkey(a1,[20,30]),@(x) inv(x));

a1w = a1.X;
for I=1:size(a1w,1)
    a1w(I,:,:) = inv(squeeze(a1w(I,:,:)));
end

