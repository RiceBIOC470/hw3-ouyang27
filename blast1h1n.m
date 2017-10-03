function zz = blast1h1n(Accession)
acc_data = topNblast(Accession);
N = length(acc_data);
human = '';
non_human = '';
for ii = 1:N
    acc = cell2mat(acc_data(ii));
    data = getgenbank(acc);
    species = data.Source;
    if contains(species, 'Homo sapiens') == 1
        human = acc;
        break
    end
end 
for ii = 1:N
    acc = cell2mat(acc_data(ii));
    data = getgenbank(acc);
    species = data.Source;
    if contains(species, 'Homo sapiens') == 0
        non_human = acc;
        break
    end
end
if length(human) == 0
    disp('None human match found');
end
zz = struct('human', human, 'non_human', non_human);

    