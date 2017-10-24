% GB comments
1.	100
2a. 30 Code doesn’t read in data correctly and you didn’t calculate fraction of shared base pair 
2b. 30 No calculations 
2c. 30 No calculations
3a 100 
3b. 100
3c. 100  	
Overall: 70

%HW3

%% Problem 1 - Smith-Waterman alignment
% Consider two sequences 'GTAATCC' and 'GTATCCG'

% Construct the scoring matrix for this with the parameters:
% match value = 2, mismatch value = -1, and gap penalty = -1. Use your
% solution to get the optimal alignment. If you prefer, it is acceptable to do this with
% pencil and paper, you can then take a snapshot of your solution and
% include it in your repository. 

% Yu Ouyang answer:
%    -  G  T  A  A  T  C  C
%  - 0  0  0  0  0  0  0  0
%  G 0  2  1  0 -1  0  0  0 
%  T 0  1  4  3  2  4  3  2
%  A 0  0  3  6  7  6  5  4
%  T 0 -1  2  5  5  9  8  6
%  C 0 -2  1  4 -2  4 11 10
%  C 0 -3  0  3  3 -3  6 13
%  G 0 -1 -1  2 -4 -1 -4  5

% align:  -GTAATCC -
%          ||| |||
%         -GTA-TCCG-

%% Problem 2 - using the NCBI databases and sequence alignments

% Erk proteins are critical signal transducers of MAP kinase signaling.
% Accessions numbers for ERK1 (also called MAPK3) and ERK2 (also called MAPK1) human mRNA are NM_002746 and
% NM_002745, respectively. 

% Part 1. Perform an alignment of the coding DNA sequences of ERK1 and
% ERK2. What fraction of base pairs in ERK1 can align to ERK2? 
genbank_dat1 = genbankread('ERK1_sequence.gb.txt');
seq1 = genbank_dat1.Sequence;

genbank_dat2 = genbankread('ERK2_sequence.gb.txt');
seq2 = genbank_dat2.Sequence;
[score, align, start] = swalign(seq1,seq2,'Alphabet','nt','Showscore',true);

% Part2. Perform an alignment of the aminoacid sequences of ERK1 and ERK2.
% What fraction of amino acids align?

info1 = getgenbank('NM_002746');
protein1 = info1.CDS;
aa1 = protein1.translation;
info2 = getgenbank('NM_002745');
protein2 = info2.CDS;
aa2 = protein2.translation;
[score_a, align_a, start_a] = swalign(aa1,aa2,'Alphabet','aa','Showscore',true);


% Part 3.  Use the NCBI tools to get mRNA sequences for the mouse genes ERK1 and
% ERK2 and align both the coding DNA sequences and protein sequences to the
% human versions. How similar are they? 
%%
ERK1 = 'NM_002746';
ERK2 = 'NM_002745';
gb_data1= getgenbank(ERK1);
seq_begin1 = gb_data1.Sequence(1:500);
[requestID1,requestTime1] = blastncbi(seq_begin1,'blastn', 'Database', 'est_mouse');
blast_data1 = getblast(requestID1, 'WaitTime', requestTime1);

gb_data2= getgenbank(ERK2);
seq_begin2 = gb_data2.Sequence(1:500);
[requestID2,requestTime2] = blastncbi(seq_begin2,'blastn', 'Database', 'est_mouse');
blast_data2 = getblast(requestID2, 'WaitTime', requestTime2);

mseq1 = blast_data1.Hits(1).Name
mseq2 = blast_data2.Hits(1).Name
%%
% Get the accession for ERK1 mouse gene is 'AV498280' and ERK2 mouse gene
% 'AV498280'
mERK1 = getgenbank('AV498280');
mseq1 = mERK1.Sequence;
mERK2 = getgenbank('AV498280');
mseq2 = mERK2.Sequence;
[score_m, align_m, start_m] = swalign(mseq1,mseq2,'Alphabet','nt','Showscore',true);

%% Problem 3: using blast tools programatically

% Part 1. Write a function that takes an NCBI accession number and a number N as input and
% returns a cell array of the accession numbers for the top N blast hits. 

z1 = topNblast('AV498280', 3)


% Part 2. Write a function that takes an accession number as input, calls your function 
% from part 1, and returns two outputs - the closest scoring match in human DNA/RNA and the
% closest non-human match. Hint: see the "Source" or "SourceOrganism" field in the data
% returned by getgenbank. Make sure your function does something sensible
% if nothing human is found. 
%%
z2 = blast1h1n('AV498280')


% Part 3. Choose at least one gene from the human genome and one gene from
% another organism and run your code from part 2 on the relevant accession
% numbers. Comment on the results. 
%%
h = blast1h1n('NM_002746') % gene from the human genome
n = blast1h1n('AV498280') % gene from the other species

% Yu Ouyang comment: h.human = 'KJ534878', h.non-human = 'AB385073';
% n.human  = 'KJ534878', n.non_human = 'NM_011952'. The non_human gene is
% the homolog of 'NM_002746' from mouse genome. Both gene map the same
% human homolog, but different non-human homolog. 
