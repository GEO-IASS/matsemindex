function Q = fillwithindex(Q)

assert(max(size(Q)) < 10,'Test function only for dimensions < 10');
S = [];
S.type = '()';
S.subs = cell(1,ndims(Q));

Q = recfill(Q,S,1,0);

function Q = recfill(Q,S,d,b)
    if d == ndims(Q)
        for I=1:size(Q,d)
            S.subs{d} = I;
            Q = subsasgn(Q,S,10*b+I);
        end
    else
        for I=1:size(Q,d)
            S.subs{d} = I;
            Q = recfill(Q,S,d+1,10*b+I);
        end        
    end
