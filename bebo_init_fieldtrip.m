function bebo_init_fieldtrip(C)
% Add FieldTrip (grid version first, local copy for testing) and initialise it.
if exist('ft_defaults', 'file') ~= 2
    for k = 1:numel(C.fieldtrip)
        if exist(C.fieldtrip{k}, 'dir'), addpath(C.fieldtrip{k}); break, end
    end
end
ft_defaults
end
