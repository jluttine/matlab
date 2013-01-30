function decades_equal(hAxes,xLimits,yLimits)

  if (nargin < 2) || isempty(xLimits)
    xLimits = get(hAxes,'XLim');
  end
  if (nargin < 3) || isempty(yLimits)
    yLimits = get(hAxes,'YLim');
  end

  logScale = diff(yLimits)/diff(xLimits);
  powerScale = diff(log10(yLimits))/diff(log10(xLimits));

  set(hAxes,'Xlim',xLimits,...
            'YLim',yLimits,...
            'DataAspectRatio',[1 logScale/powerScale 1]);

end
