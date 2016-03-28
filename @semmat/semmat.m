classdef semmat
    properties
        X
        akeys
    end
    methods
        function self = semmat(X,keys)
            self.X = X;
            if nargin == 1
                self.akeys = [];
            else
                self.akeys = keys;
            end
        end
        function k = keys(self)
            k = self.akeys;
        end
        function s = size(self,index)
            if nargin == 1
                s = size(self.X);
            else
                s = size(self.X,index);
            end
        end
        function r = subsasgn(self,S,B)
            if strcmp(S.type,'()')
                oS = subsrefcommon(self,S);
                if isa(B,'semmat')
                    error('not implemented assign of semmat');
                else
                    self.X = subsasgn(self.X,oS,B); 
                end
                r = self;
            else
                error('unsupported');
            end
        end
        function oS = subsrefcommon(self,S)
                N = ndims(self.X);
                ii = [];
                for I=1:length(S.subs)
                    if ischar(S.subs{I})
                        error('Column indexing not needed in semmat');
                    else
                        w = S.subs{I};
                        if isempty(self.akeys)
                            for J=1:size(w,2)
                                ii(end+1) = w(J).dim;
                            end
                        else
                            % we need to match dim with the key, not with
                            % the effective index, this allows for
                            % generalized case
                            for J=1:size(w,2)
                                ii(end+1) = find(self.akeys == w(J).dim,1,'first');
                            end
                        end
                    end
                end
                [iio] = unique(ii,'sorted'); % iia is index in S.subs
                if length(iio) < length(ii)
                    error('Index used multiple times');
                end
                if iio(1) < 1 || iio(end) > N
                    error('Out of dimension rang');
                end
                oS = [];
                oS.type = '()';
                oS.subs = cell(1,N);
                for I=1:N
                    oS.subs{I} = ':'; 
                end
                k = 1;
                for I=1:length(S.subs)
                    w = S.subs{I};
                    for J=1:size(w,2)
                        oS.subs{ii(k)} = w(J).index; % resolved key
                        k = k+1;
                    end
                end            
        end
        function r = valuesubsref(self,dimindexes,indexes)  
                N = ndims(self.X);
                oS = [];
                oS.type = '()';
                oS.subs = cell(1,N);
                for I=1:N
                    oS.subs{I} = ':'; 
                end
                for I=1:length(dimindexes)
                    oS.subs{dimindexes(I)} = indexes(I);
                end
                b = subsref(self.X,oS);
                ss = size(b);
                if numel(b) == 1
                    r = b;
                else
                    % if there is any singleton dimension
                    if ~isempty(find(ss == 1, 1))
                        r = squeeze(b);
                    else
                        r = b;
                    end
                end
        end
        function r = valuesubsasgn(self,dimindexes,indexes,x)  
            N = ndims(self.X);
            oS = [];
            oS.type = '()';
            oS.subs = cell(1,N); % all ':'
            for I=1:N
                oS.subs{I} = ':'; 
            end
            for I=1:length(dimindexes)
                oS.subs{dimindexes(I)} = indexes(I); % 
            end
            self.X = subsasgn(self.X,oS,x);
            r = self;
        end
        function r = subsref(self,S)
            if strcmp(S.type,'()')
                oS = subsrefcommon(self,S);
                b = subsref(self.X,oS);
                ss = size(b);
                % subsref(...) eq scalar -> scalar -> no key
                % squeeze(subsref) eq 1x1x1xK -> Kx1
                if numel(b) == 1
                    self.akeys = [];
                    self.X = b;
                elseif isempty(self.akeys)
                    self.X = squeeze(b);
                else
                    % if there is any singleton dimension
                    if ~isempty(find(ss == 1, 1))
                        self.X = squeeze(b);
                        self.akeys = self.akeys(ss > 1);
                    else
                        % same keys (although reduction of size along a
                        % dimension makes it invalid)
                        self.X = b;
                    end
                end
                r = self;                
            else
                if strcmp(S.type,'.') && strcmp(S.subs,'X')
                    r = self.X;
                else
                    error('Unsupported subsref');
                end
            end
        end
        function it = aiterateall(self)
            it = semmat.iterateall(self,size(self.X));
        end
        function idxs = amultiindex(self,ii)
            idxs = semmat.multiindex(ii,self.keys);
        end
        function d = data(self)
            d = self.X;
        end
        function d = double(self)
            d = self.X;
        end
        function d = indexofkey(self,k)
            d = arrayfun(@(x) find(self.akeys==x,1,'first'),k);
        end
        function d = sizeofkey(self,k)
            d = arrayfun(@(x) size(self.X,find(self.akeys==x,1,'first')),k);
        end
                    
        % permute the content by keys
        function self = permute(self,neworderkeys)
            if isempty(self.akeys)
                self.X = permute(self.X,neworderkeys);
            else
                if(~all(neworderkeys == self.akeys))
                    self.X = permute(self.X,indexofkey(self,neworderkeys));
                    self.akeys = neworderkeys;
                end
            end
        end        

        % unary operation over submatrix specified by keys dimensions
        % the output is a matrix:
        %   R = allkeys - keys
        %   out = [R,newkeys]
        %
        % e.g. inverse of submatrix in given pair of keys
        function r = unaryop(self,keys,outkeys,outsizes,fun)

            dk = setdiff(self.keys,keys,'stable'); % broadcast keys
            self = permute(self,[dk,keys]); % keys order is good
            di = indexofkey(self,dk);
            ds = size(self.X,di);
            % build output
            rk = [dk,outkeys];
            rs = [ds,outsizes];
            r = semmat(zeros(rs),rk); % combine broadcast with specified output
            % and broadcast
            oi = 1:length(dk);
            ki = semmat.iterateall(ds); % extrat all indices from the broadcast - note preserved
            for I=1:size(ki,1)
                %qi = semmat.multiindex(ki(I,:),keys);
                %w = double(self(qi));
                %out(qi) = fun(w)
                r = valuesubsasgn(r,oi,ki(I,:),fun(valuesubsref(self,di,ki(I,:))));
            end
        end

        % applies a binary operation over the submatrices of self(a1) and a2 identified by the common keys
        % these submatricies can be of any dimension. 
        % 
        % Let's take the keys of self(a1) and a2 and tripartite as: s12,s1,s2 with s12 in common
        % The output is specified as: s1,outkeys,s2
        %
        % If the function is sum(reshape(x.*y,1,[])) then this generalizes the matrix product: all the same key dimension are memberwise multiplied and summed to scalar
        % Regular 2D case using notation (key,size): [(1,10),(2,20)] ** [(3,40),(2,20)] then gives [(1,10),(3,40)]
        %
        function r = binaryop(self,a2,outkeys,outsizes,fun)
            a1 = self;
            [k1,k2,k12] = semmat.splitkeys(keys(a1),keys(a2));
            s1 = sizeofkey(a1,k1);
            s2 = sizeofkey(a2,k2);
            s12_1 = sizeofkey(a1,k12);
            s12_2 = sizeofkey(a2,k12);
            assert(all(s12_1==s12_2),'same sizes in common')

            rk = [k1,outkeys,k2];
            rs = [s1,outsizes,s2];
            r = semmat(zeros(rs),rk);

            a1 = permute(a1,[k1,k12]); % adjust order of common
            a2 = permute(a2,[k2,k12]); % adjust order of common
            i1 = indexofkey(a1,k1);
            i2 = indexofkey(a2,k2);
            oi1 = indexofkey(r,k1);
            oi2 = indexofkey(r,k2);

            ai1 = semmat.iterateall(s1); % broadcast iteration
            ai2 = semmat.iterateall(s2); % broadcast iteration

            for I=1:size(ai1,1)
                %called outside classef:
                %i1 = semmat.multiindex(ai1(I,:),k1); 
                %a1(i1)
                w1 = valuesubsref(a1,i1,ai1(I,:)); % now indexed by s12 => same order
                % one loop in a2 only
                for J=1:size(ai2,1)
                    %i2 = semmat.multiindex(ai2(J,:),k2);
                    %w2 = a2.subsref(i2); % now indexed by s12 => same order
                    w2 = valuesubsref(a2,i2,ai2(J,:));
                    r = valuesubsasgn(r,[oi1,oi2],[ai1(I,:),ai2(J,:)],fun(w1,w2));
                end
            end
        end
    end    
    methods(Static)
        function [ss,ke] = unionsizekey(ss1,ke1,ss2,ke2)
            [ii,ii1,ii2] = intersect(ke1,ke2);
            dii2 = setdiff(1:length(ss2),ii2);
            assert(ss1(ii1) == ss2(ii2),'all common sizes should be same');
            ss = [ss1,ss2(dii2)];
            ke = [ke1,ke2(dii2)];
        end
        function [ss,ke] = diffkeys(ss1,ke1,ke2)
            % we use intersect to detect size mistakes, otherwise setdiff
            % would suffice
            [ii,ii1,ii2] = setdiff(ke1,ke2);
            ss = ss1(ii1);
            ke = ii;
        end
        function [ss,ke] = diffsizekeys(ss1,ke1,ss2,ke2)
            % we use intersect to detect size mistakes, otherwise setdiff
            % would suffice
            [ii,ii1,ii2] = intersect(ke1,ke2);
            dii1 = setdiff(1:length(ss1),ii1);
            assert(ss1(ii1) == ss2(ii2),'all common sizes should be same');
            ss = ss1(dii1);
            ke = ke1(dii1);
        end
        function [s1,s2,s12,i1,i2] = splitkeys(ke1,ke2)
            s12 = intersect(ke1,ke2);
            [s1,i1] = setdiff(ke1,s12,'stable');
            [s2,i2] = setdiff(ke2,s12,'stable');            
        end
        function a = iterateall(ss1)
            if isempty(ss1)
                a = 1;
            else
                % adapted from allcombi in FileExchange
                NC = length(ss1);
                q = arrayfun(@(x) 1:x,ss1,'UniformOutput',0);
                ii = NC:-1:1;
                [A{ii}] = ndgrid(q{ii});
                a = reshape(cat(NC+1,A{:}),[],NC);
            end
        end
        % given values i1...in and associated keys builds i1..in semindex
        function r = multiindex(indices,keys)
            % arrayfun does not work here
            r = semindex(keys(1),indices(1));
            for J=2:length(indices)
                r = [r; semindex(keys(J),indices(J))];
            end
            
        end
    end
end

