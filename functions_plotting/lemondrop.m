function [] = lemondrop()

smallN = 256;
smallN2 = 258;
graymax = 240;
range = 0:1/(smallN+1):1;
muLD = range*(smallN + 1) / (smallN + 1);
lemonDrop = sqrt(muLD.*(1-muLD)) * 3.4;
sampleCount2 = (smallN + 2):-1:1;
colorRange = (255-graymax)*sampleCount2/(smallN + 2);
base = repmat(graymax, smallN + 2, 1);
col = (base + colorRange') / 255;
rgb = [col col col];
count2 = 1;
for ii = ceil(smallN2/2):smallN2-1
    ix = [ii ii+1 smallN2-ii smallN2-ii+1];
    fill(range(ix), lemonDrop(ix), rgb(count2, :),'edgecolor','none')
    fill(range(ix), -lemonDrop(ix), rgb(count2, :),'edgecolor','none')
    count2 = count2 + 2;
end
plot(muLD,lemonDrop,'k--');
plot(muLD,-lemonDrop,'k--');

end