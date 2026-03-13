function hfig = tightfig(hfig)
    if nargin == 0, hfig = gcf; end
    set(hfig, 'PaperPositionMode', 'auto');
end