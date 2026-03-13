function hfig = tightfig(hfig)
%tightfig strips excess whitespace from margins in figures
    if nargin == 0, hfig = gcf; end
    set(hfig, 'PaperPositionMode', 'auto');
end
