function zz = topNblast(Accession, N)
gb_data= getgenbank(Accession);
seq_begin = gb_data.Sequence(1:500);
[requestID,requestTime] = blastncbi(seq_begin,'blastn');
blast_data = getblast(requestID, 'WaitTime', requestTime);
if nargin == 2
    N = N;
else
    N = length(blast_data.Hits);
end
for ii = 1:N
    c = blast_data.Hits(ii);
    i = strfind(c.Name, '|');
    i = i(3);
    e = strfind(c.Name,'.');
    acc = c.Name(i+1:e-1);
    disp(acc);
    z(ii) = cellstr(acc);
end
zz = z;