function colors = get_7_default_colors()

figure(3768)
a = plot([2 3]);
hold on
b = plot(2*[2 3]);
c = plot(3*[2 3]);
d = plot(4*[2 3]);
e = plot(5*[2 3]);
f = plot(6*[2 3]);
g = plot(7*[2 3]);

% colors = [a.Color; b.Color; c.Color; d.Color; e.Color; f.Color; g.Color];
colors = {a.Color; b.Color; c.Color; d.Color; e.Color; f.Color; g.Color};
% colors = [a.Color; b.Color; c.Color; d.Color; e.Color; [0 0 0]; g.Color];
% colors = kron(colors,[1; 1]);
% colors = kron(colors,[1; 1; 1; 1]);

close

end