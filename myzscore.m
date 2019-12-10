function zscore = myzscore(x, xref, dim)
% returns pairwise zscores operated on columns independednlty
szx = size(x);
if ndims(x)>2
    x = reshape(x, [szx(1), prod(szx(2:end))]);
end

if nargin < 2
    xref = zeros(size(x));
    dim = 1;
end
if isempty(xref)
    xref = zeros(size(x));
end    
szxref = size(xref);
if ndims(xref)>2
    xref = reshape(xref, [szxref(1), prod(szxref(2:end))]);
end
if nargin == 3 && dim == 2
    zscore = myfun(x.', xref.');
    zscore = zscore.';
else
    zscore = myfun(x, xref);
end
zscore = reshape(zscore, szx);


end
    function zscore = myfun(x, xref)
        zscore = zeros(size(x));
        [n, p] = size(x);
        [nref, pref] = size(xref);
        muref = mean(xref);
        sref = std(xref);
        for i = 1:p
           zscore(:,i) = (x(:,i) - muref(i))/sref(i); 
        end
    end

    
