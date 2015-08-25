function frameNo = getFrameNumbersFromMissing(firstFrame, lastFrame, missingFrames)

nGaps = size(missingFrames, 2);
n_missing = sum(missingFrames(2,:));

frameNo = firstFrame:(lastFrame-n_missing);
for g = 1:nGaps
    startIdx = find(frameNo == missingFrames(1,g));
    frameNo(startIdx:end) = frameNo(startIdx:end) + missingFrames(2,g);
end