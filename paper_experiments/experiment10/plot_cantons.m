function plot_cantons(S,hidd_cant,clstr,n)
    O = size(S,1);
    N = 26; 
    R_earth = 6371;
    cantons = 1:N;
    cantons(hidd_cant) = [];
    obs_idx = cantons;
    %hacer que funcione para distinto número de cantones
    
    col_full_names = {'Zurich','Bern','Luzern','Uri','Schwyz','Obwalden','Nidwalden','Glarus','Zug','Fribourg','Solothurn','Basel-Stadt','Basel-Land','Schaffhausen','App Aus','App Inn','St. Gallen','Graubunden','Argau','Thurgau','Ticino','Vaud','Valais','Neuchatel','Geneve','Jura'};
    Longitudes = [8.656389, 7.616667,  8.116667,  8.616667,  8.75,       8.033333,  8.4,  9.066667,  8.548,   7.083333,  7.633333, 7.6,       7.755833,  8.566667,   9.3128,  9.4,   9.166667,  9.5,  8.11,   9.066667, 8.816667, 6.55, 7.6, 6.783333, 6.166, 7.15];
    Latitudes = [47.413056, 46.833333, 47.083333, 46.783333, 47.066667, 46.866667, 46.95, 46.983333, 47.1648, 46.716667, 47.15,    47.566667, 47.463056, 47.716667,  47.3784, 47.3, 47.333333, 46.75, 47.418, 47.583333, 46.316667, 46.616667, 46.066667, 46.983333, 46.218, 47.366667];
    
    cfn = col_full_names(obs_idx);
    lng = Longitudes(obs_idx);
    lat = Latitudes(obs_idx);

    Coords = zeros(O,2);
    Coords(:,1) = lng*pi/180;
    Coords(:,2) = lat*pi/180;
    
    Coords_km = zeros(O,2);
    Coords_km(:, 1) = R_earth*Coords(:, 1).*cos(Coords(:, 2));
    Coords_km(:, 2) = R_earth*Coords(:, 2);
    Coords_km = round(Coords_km);
        
    S_bin = mbinarize(S,2);%Adj binarizada 
    S_hat = S.*S_bin; % Adj binarizada con pesos

    if n == 1
        S_plot = S_bin;
    else
        S_plot = S_hat;
    end
    G = graph(S_plot, cfn);
    AP = plot(G,'XData',Coords_km(:,1),'YData',Coords_km(:,2),'LineWidth',G.Edges.Weight/max(G.Edges.Weight)*5); 
    %hacer la comparación de clusters con el groundtruth
    highlight(AP,find(clstr(:,n)==1),'NodeColor','r');
    highlight(AP,find(clstr(:,n)==2),'NodeColor','blue');
    highlight(AP,find(clstr(:,n)==3),'NodeColor','yellow');
    AP.MarkerSize = 10;
    
%     subplot(1,2,2)
%     G = graph(S_hat, cfn);
%     AP = plot(G,'XData',Coords_km(:,1),'YData',Coords_km(:,2),'LineWidth',G.Edges.Weight*10); 
%     %hacer la comparación de clusters con el groundtruth
%     highlight(AP,find(clstr(:,2)==1),'NodeColor','r');
%     highlight(AP,find(clstr(:,2)==2),'NodeColor','blue');
%     highlight(AP,find(clstr(:,2)==3),'NodeColor','yellow');
%     AP.MarkerSize = 10;
end