function seed_gen(num)
for n = 1:num
    a = clock;
    sd = a(6)*1e6;
    
    sd = num2str(sd);
    sd = fliplr(sd);
    sd = str2double(sd);
    rng(sd)
    sd = rand;
    sdo(n) = (round(sd*1e6));
    pause(rand)
end

filename = 'Seed.txt';
prev = dlmread(filename);
new = [prev;sdo'];
dlmwrite(filename,sdo', 'precision', '%.0f')
% fileID = fopen(filename,'w');
% fprintf(fileID,'%s',sd);

end