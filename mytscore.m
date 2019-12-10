function tscore = mytscore(x, xref, dim)
% returns pairwise t-scores operated on columns independednlty
if nargin < 2
    xref = zeros(size(x));
    dim = 1;
end
if isempty(xref)
    xref = zeros(size(x));
end    
if nargin == 3 && dim == 2
    tscore = myfun(x.', xref.');
    tscore = tscore.';
else
    tscore = myfun(x, xref);
end

end
    function tscore = myfun(x, xref)
        [n, p] = size(x);
        [nref, pref] = size(xref);
        muref = mean(xref);
        sref = std(xref);
        mu = mean(x);
        s = std(x);
        tscore = (mu - muref)./sqrt(s.*s/n + sref.*sref/nref); 
    end

    
