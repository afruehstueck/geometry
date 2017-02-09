
[V, F] = read_vertices_and_faces_from_obj_file('../data/teddy.obj');
%[V, F] = read_vertices_and_faces_from_obj_file('../data/cow.obj');
%[V, F] = read_vertices_and_faces_from_obj_file('../data/pumpkin.obj');
%[V, F] = read_vertices_and_faces_from_obj_file('../data/teapot.obj');
scr = get(0, 'ScreenSize');  

fig = figure('Name', '2D Viewer', 'NumberTitle', 'off', 'Position', [scr(3)/2 50 scr(3)/3 scr(3)/3]);
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3));%,'FaceColor',[0.26,0.33,1.0 ]);

light('Position', [0, 60.0, 30.0], 'Style', 'infinite');
light('Position', [0, -60.0, -30.0], 'Style', 'infinite');

lighting phong;