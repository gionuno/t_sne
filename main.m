DIR = dir('train_set');
X = [];
C = [];
for i = 3:length(DIR)
    dirname = strcat('train_set/',DIR(i).name);
    disp(dirname);
    AUX_DIR = dir(dirname);
    for j = 3:length(AUX_DIR)
        auxfilename = strcat(dirname,'/',AUX_DIR(j).name);
        aux = double(imread(auxfilename))/255.0;
        d = size(aux);
        X = [X ; reshape(aux,[1,numel(aux)])];
        C = [C ; i-3];
    end
end

mX = mean(X,1);
Y = X - repmat(mX,[size(X,1),1]);
sX = sqrt(mean(Y.^2,1));
Y = Y ./ repmat(sX,[size(X,1),1]);
[Z,S,P,E] = tSNE(Y,3,40.0,500);

figure;
plot(log(E+1e-8));

figure;
for i = 3:length(DIR)
    disp(DIR(i).name);
    idx = (C == i-3);
    scatter3(Z(idx,1),Z(idx,2),Z(idx,3),7,C(idx),'filled','MarkerEdgeColor','k','DisplayName',DIR(i).name);
    if i == 3
        hold on;
    end
end
colormap(colorcube);
legend(gca,'show');
hold off;