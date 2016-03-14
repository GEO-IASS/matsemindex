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

%% marginalization 
Q = fillwithindex( zeros([2,3,4]));
a = semmat(Q,[20,30,40]);
ki = 30;
[ss,ke] = semmat.diffkeys(size(a),keys(a),ki);
ac = semmat(zeros(ss),ke);
% iterate over ac keys 
%  c = data(a(ii)); % this is sized sizeofkey(ki)
%  ac(ii) = sum(c(:));

%% marginalization except
Q = fillwithindex( zeros([2,3,4]));
a = semmat(Q,[20,30,40]);
ki = 30;
kidx = indexofkey(a,ki);
ac = semmat(zeros(size(a,kidx),1),ki); % 
for I=1:size(a,kidx)
    i = semindex(ki,I);
    c = data(a(i));
    %ac(i) = sum(c(:));
end

%% product
Q1 = fillwithindex( zeros([2,3,4]));
Q2 = fillwithindex( zeros([2,3,4]));
a1 = semmat(Q1,[20,30,40]);
a2 = semmat(Q2,[40]);

[s1,s2,s12] = semmat.splitkeys(keys(a1),keys(a2));
[ss,ke] = semmat.unionsizekey(size(a1),keys(a1),size(a2),keys(a2));

ac = semmat(zeros(size(ss)),ke);
% TODO reorder a1 and a2 with same order of ac, using permute

ai1 = semmat.iterateall(sizeofkey(a1,s1));
ai2 = semmat.iterateall(sizeofkey(a2,s2));
for I=1:size(ai1,1)
    % make index i1
    w1 = a1(i1); % now indexed by s12 => same order
    for J=1:size(ai2,2)
        % make index i2
        w2 = a2(i2); % now indexed by s12 => same order
        % ac(i1,i2) = w1 * w2
    end
end



