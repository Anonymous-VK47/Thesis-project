function mkfld()

    directory = fullfile("C:\", 'Multibody Dynamics Generator');
    if ~exist(directory, 'dir')
        mkdir(directory);
    end
end