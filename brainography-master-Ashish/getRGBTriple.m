function rgbVals = getRGBTriple(rawMap,minVal,maxVal,entriez)

% rgbVals = zeros(size(entriez,1),3);

valrange = maxVal-minVal;

entriez(find(entriez < minVal)) = minVal;
entriez(find(entriez > maxVal)) = maxVal;

if valrange > 0
    entriez(isnan(entriez)) = 0;
    mapIdx = 1 + round((size(rawMap,1)-1)*(entriez - minVal)/valrange);
    mapIdx(find(mapIdx < 0)) = 1;
    rgbVals = rawMap(mapIdx,:);
else
    rgbVals = repmat(rawMap(size(rawMap,1)/2,:),size(entriez,1),1);
end
end