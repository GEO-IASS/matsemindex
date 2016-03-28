
%% marginalization (output size is all the )
Q = fillwithindex( zeros([2,3,4]));
a = semmat(Q,[20,30,40]);
ki = 30;
[ss,ke] = semmat.diffkeys(size(a),keys(a),ki);
ac = semmat(zeros(ss),ke);

% iterate over ac keys 
aic = iterateall(ac);
for I=1:size(aic,1)
    ic = multiindex(ac,aic(I,:));
    wc = data(a1(ic)); 
    ac(ic) = sum(wc(:));
end

%% marginalization except (output size is variable size)
Q = fillwithindex( zeros([2,3,4]));
a = semmat(Q,[20,30,40]);
ki = 30;
kidx = indexofkey(a,ki);
ac = semmat(zeros(size(a,kidx),1),ki); % 
for I=1:size(a,kidx)
    i = semindex(ki,I);
    c = data(a(i));
    ac(i) = sum(c(:));
end

%% product
Q1 = fillwithindex( zeros([2,3,4]));
Q2 = fillwithindex( zeros([2,3,4]));
a1 = semmat(Q1,[20,30,40]);
a2 = semmat(Q2,[40]);

%%
[s1,s2,s12] = semmat.splitkeys(keys(a1),keys(a2));
[ss,ke] = semmat.unionsizekey(size(a1),keys(a1),size(a2),keys(a2));

% build output 
ac = semmat(zeros(ss),ke);
% and reorder inputs to enforce same order in output
a1 = permute(a1,intersect(keys(ac),keys(a1),'stable'));
a2 = permute(a2,intersect(keys(ac),keys(a2),'stable'));
% build all the iterations
ai1 = semmat.iterateall(sizeofkey(a1,s1));
ai2 = semmat.iterateall(sizeofkey(a2,s2));

% one loop in a1 only
for I=1:size(ai1,1)
    i1 = semmat.multiindex(ai1(I,:),s1);
    w1 = a1(i1); % now indexed by s12 => same order
    % one loop in a2 only
    for J=1:size(ai2,2)
        i2 = semmat.multiindex(ai2(J,:),s2);
        w2 = a2(i2); % now indexed by s12 => same order
        ac(i1,i2) = double(w1) .* double(w2);        % because w1 and w2 have EXACT same layout
    end
end
