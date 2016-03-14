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
        
        function d = data(self)
            d = self.X;
        end
        function d = double(self)
            d = self.X;
        end
        function d = indexofkey(self,k)
            d = find(self.akeys==k,1,'first');
        end
        function d = sizeofkey(self,k)
            d = zeros(length(k),1);
            for I=1:length(k)
                d(I) = size(self.X,indexofkey(k(I)));
            end
        end
            
        
        % permute the content for some reason => means permute keys
        function self = permute(self,neworder)
            self.X = permute(self.X,neworder);
            if ~isempty(self.akeys)
                self.akeys = self.akeys(neworder);
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
        function [s1,s2,s12] = splitkeys(ke1,ke2)
            s12 = intersect(ke1,ke2);
            s1 = setdiff(ke1,s12);
            s2 = setdiff(ke2,s12);            
        end
        function a = iterateall(ss1)
            o = zeros(prod(ss1),length(ss1));
            k = 1;
            n = size(o,1);
            % ?
%            for I=1:size(o,2)
 %               o(k
  %              k = k * size(o,I);
   %         end
        end
        % given values i1...in and associated keys builds i1..in semindex
        function r = multiindex(indices,keys)
            r = semindex(keys(1),indices(1));
            for J=2:length(indices)
                r = [r; semindex(keys(J),indices(J))];
            end
            
        end
    end
end
end

