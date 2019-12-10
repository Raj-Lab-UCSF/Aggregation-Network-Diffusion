function outV = rescale_braindata(V, sig, method, whichcols)
% rescales brain data to convert t-stats or z-scores to more well-behaved
% numbers that refect atrophy or similar brain quantities
% Three methods: 'logistic', 'stdev' and 'explinear'
% sig is ratios wrt stdev of data
% Note there is a small center shift of 0.5*std to bias output towards
% above-mean values
% Works on both vectors and matrices (columnwise, unless last argument is 'allcols')

outV = V;
if nargin==4 && strcmp(whichcols, 'allcols')
    outV = V(:);
end

sz = size(outV);
for i = 1:sz(2)
    v = outV(:,i);
    stdv = std(v);
    if strcmp(method, 'logistic')
        center = 1*sig; %0.5*sig;
        if isempty(find(v>0))
                a0 = mean(v) + center*stdv;
        else
            a0 = mean(v(v>0)) + center*stdv;  
        end
        outv = 1./(1+exp(-(v-a0)/sig/stdv));
    elseif strcmp(method, 'stdev')
        a0 = mean(v(v>0));  
        outv = max(v, a0-sig*stdv);
        outv = min(outv, a0+sig*stdv);
    elseif strcmp(method,'explinear')
        center = 0; %0.5*sig;
        a0 = mean(v(v>0)) + center*stdv;
        outv = v;
        outv(v<a0) = a0*exp((v(v<a0)-a0)/sig/stdv);
    end
    outV(:,i) = outv;
end

if nargin==4 && strcmp(whichcols, 'allcols')
    outV = reshape(outV, size(V));
end

