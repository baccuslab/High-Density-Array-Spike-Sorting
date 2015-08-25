function AmplitudeSelectableAxesTest()
    r = [-10 10 50 100];
    ah = guiutil.AmplitudeSelectableAxes(@cb, 'clickableTicks', r, 'markerWidth', 10);
    plot(ah, r, ones(1, length(r)), 'kx');
    
    function b = cb(c, a)
        disp([c a])
        b = a<1;
    end
end