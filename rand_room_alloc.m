%%% Room 1: Rev Room
%%% Room 2: Anechoic
%%% Room 3: Classroom
%%% Room 4: Lect Hall

%%% Song 1: Staccato Piece EIL
%%% Song 2: Legato Piece ATR
%%% Song 3: Own Piece

%%% Subject order is order that I meet them in

num_subs = 5;
num_rms = 4;
num_songs = 3;

%%% Randomise Rooms
rooms = zeros(num_subs, num_rms);
songs_r1 = zeros(num_rms, num_songs);
songs_r2 = zeros(num_rms, num_songs);
songs_r3 = zeros(num_rms, num_songs);
songs_r4 = zeros(num_rms, num_songs);
songs_r5(rm,:) = randperm(s,num_songs);
s = RandStream('mt19937ar','Seed',0);

for subject = 1:num_subs
    rooms(subject,:) = randperm(s,num_rms);
end

for rm = 1:num_rms
    songs_r1(rm,:) = randperm(s,num_songs);
    songs_r2(rm,:) = randperm(s,num_songs);
    songs_r3(rm,:) = randperm(s,num_songs);
    songs_r4(rm,:) = randperm(s,num_songs);
    songs_r5(rm,:) = randperm(s,num_songs);
end


display('Room Order');
rooms

display('  ');
display('Songs Sub 1');
songs_r1

display('  ');
display('Songs Sub 2');
songs_r2

display('  ');
display('Songs Sub 3');
songs_r3

display('  ');
display('Songs Sub 4');
songs_r4

display('  ');
display('Songs Sub 5');
songs_r5